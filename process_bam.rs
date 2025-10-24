use anyhow::{Context, Result};
use cr_types::reference::ReferenceInfo;
use rust_htslib::bam::{self, Read, Record};
use std::path::Path;
use tx_annotation::transcript::{AnnotationParams, TranscriptAnnotator, AnnotationRegion};
use fastq_set::WhichEnd;
use cr_types::ReqStrand;

fn main() -> Result<()> {
    // 命令行参数解析
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("Usage: {} <reference_path> <input_bam> <output_bam>", args[0]);
        std::process::exit(1);
    }

    let reference_path = Path::new(&args[1]);
    let input_bam = &args[2];
    let output_bam = &args[3];

    // 创建转录组注释器
    let params = AnnotationParams {
        chemistry_strandedness: ReqStrand::Forward, // 可以根据实际情况调整
        chemistry_endedness: WhichEnd::FivePrime,   // 可以根据实际情况调整
        intergenic_trim_bases: 5,
        intronic_trim_bases: 5,
        junction_trim_bases: 5,
        region_min_overlap: 0.5,
        include_exons: true,
        include_introns: true,
    };

    println!("Initializing transcript annotator from {}", reference_path.display());
    let annotator = TranscriptAnnotator::new(reference_path, params)?;

    // 打开输入BAM文件
    let mut bam = bam::Reader::from_path(input_bam)
        .with_context(|| format!("Failed to open input BAM file: {}", input_bam))?;
    
    // 获取头部信息
    let header = bam::Header::from_template(bam.header());
    
    // 创建输出BAM文件
    let mut out = bam::Writer::from_path(output_bam, &header, bam::Format::Bam)
        .with_context(|| format!("Failed to create output BAM file: {}", output_bam))?;

    println!("Processing BAM file...");
    let mut total_reads = 0;
    let mut exonic_reads = 0;
    let mut intronic_reads = 0;
    let mut intergenic_reads = 0;

    // 处理每条记录
    for result in bam.records() {
        let mut record = result?;
        total_reads += 1;

        if record.is_unmapped() {
            // 直接写入未比对的记录
            out.write(&record)?;
            continue;
        }

        // 注释记录
        let annotation = annotator.annotate_alignment(&record);
        
        // 添加区域标签 (RE tag)
        if let Some(re_tag) = annotation.make_re_tag() {
            record.push_aux("RE", re_tag)?;
        }

        // 添加基因标签 (GX, GN tags)
        if let Some((gx_tag, gn_tag)) = annotation.make_gx_gn_tags() {
            record.push_aux("GX", &gx_tag)?;
            record.push_aux("GN", &gn_tag)?;
        }

        // 添加转录本标签 (TX tag)
        if let Some(tx_tag) = annotation.make_tx_tag() {
            record.push_aux("TX", &tx_tag)?;
        }

        // 添加反义转录本标签 (AN tag)
        if let Some(an_tag) = annotation.make_an_tag() {
            record.push_aux("AN", &an_tag)?;
        }

        // 统计不同区域的读数
        match annotation.region {
            AnnotationRegion::Exonic => exonic_reads += 1,
            AnnotationRegion::Intronic => intronic_reads += 1,
            AnnotationRegion::Intergenic => intergenic_reads += 1,
        }

        // 写入处理后的记录
        out.write(&record)?;

        // 每处理100万条记录打印一次进度
        if total_reads % 1_000_000 == 0 {
            println!("Processed {} million reads", total_reads / 1_000_000);
        }
    }

    println!("Processing complete!");
    println!("Total reads: {}", total_reads);
    println!("Exonic reads: {} ({:.2}%)", exonic_reads, 100.0 * exonic_reads as f64 / total_reads as f64);
    println!("Intronic reads: {} ({:.2}%)", intronic_reads, 100.0 * intronic_reads as f64 / total_reads as f64);
    println!("Intergenic reads: {} ({:.2}%)", intergenic_reads, 100.0 * intergenic_reads as f64 / total_reads as f64);

    Ok(())
}