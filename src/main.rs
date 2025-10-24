use anyhow::{Context, Result};
use clap::Parser;
use fastq_set::WhichEnd;
use rust_htslib::bam::{self, Read, Record};
use std::path::PathBuf;
use tx_annotation::transcript::{AnnotationParams, TranscriptAnnotator};

/// Tool to annotate BAM files with exon/intron information
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Path to the reference directory containing STAR index files
    #[arg(short, long)]
    reference: PathBuf,

    /// Input BAM file
    #[arg(short, long)]
    input: PathBuf,

    /// Output BAM file
    #[arg(short, long)]
    output: PathBuf,

    /// Chemistry strandedness (forward or reverse)
    #[arg(long, default_value = "forward")]
    strandedness: String,

    /// Chemistry endedness (five_prime or three_prime)
    #[arg(long, default_value = "three_prime")]
    endedness: String,

    /// Minimum overlap fraction to consider a read exonic/intronic
    #[arg(long, default_value_t = 0.5)]
    region_min_overlap: f64,

    /// Whether to include exons in annotation
    #[arg(long, default_value_t = true)]
    include_exons: bool,

    /// Whether to include introns in annotation
    #[arg(long, default_value_t = false)]
    include_introns: bool,

    /// Bases to trim from intergenic regions
    #[arg(long, default_value_t = 0)]
    intergenic_trim_bases: i64,

    /// Bases to trim from intronic regions
    #[arg(long, default_value_t = 0)]
    intronic_trim_bases: i64,

    /// Bases to trim from junction regions
    #[arg(long, default_value_t = 0)]
    junction_trim_bases: i64,
    
    /// Add gene ID (GX) and gene name (GN) tags
    #[arg(long, default_value_t = false)]
    add_gene_tags: bool,
    
    /// Add transcript (TX) tag
    #[arg(long, default_value_t = false)]
    add_tx_tag: bool,
    
    /// Add antisense transcript (AN) tag
    #[arg(long, default_value_t = false)]
    add_an_tag: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Parse chemistry strandedness
    let chemistry_strandedness = match args.strandedness.to_lowercase().as_str() {
        "forward" => cr_types::ReqStrand::Forward,
        "reverse" => cr_types::ReqStrand::Reverse,
        _ => anyhow::bail!("Invalid strandedness: {}", args.strandedness),
    };

    // Parse chemistry endedness
    let chemistry_endedness = match args.endedness.to_lowercase().as_str() {
        "five_prime" => WhichEnd::FivePrime,
        "three_prime" => WhichEnd::ThreePrime,
        _ => anyhow::bail!("Invalid endedness: {}", args.endedness),
    };

    // Create annotation parameters
    let params = AnnotationParams {
        chemistry_strandedness,
        chemistry_endedness,
        intergenic_trim_bases: args.intergenic_trim_bases,
        intronic_trim_bases: args.intronic_trim_bases,
        junction_trim_bases: args.junction_trim_bases,
        region_min_overlap: args.region_min_overlap,
        include_exons: args.include_exons,
        include_introns: args.include_introns,
    };

    println!("Initializing transcript annotator...");
    let annotator = TranscriptAnnotator::new(&args.reference, params)?;

    // Open input BAM file
    println!("Opening input BAM file: {}", args.input.display());
    let mut bam = rust_htslib::bam::Reader::from_path(&args.input)
        .with_context(|| format!("Failed to open BAM file: {}", args.input.display()))?;

    // Create output BAM file with same header
    println!("Creating output BAM file: {}", args.output.display());
    let header = bam::Header::from_template(bam.header());
    let mut out = bam::Writer::from_path(
        &args.output,
        &header,
        bam::Format::Bam,
    )
    .with_context(|| format!("Failed to create output BAM file: {}", args.output.display()))?;

    // Process each read
    println!("Processing reads...");
    let mut count = 0;
    let mut exonic_count = 0;
    let mut intronic_count = 0;
    let mut intergenic_count = 0;

    for result in bam.records() {
        let mut record = result.with_context(|| "Failed to read BAM record")?;

        // Skip unmapped reads
        if record.is_unmapped() {
            out.write(&record)?;
            continue;
        }

        // Annotate the read
        let annotation = annotator.annotate_alignment(&record);

        // Add region tag (RE) - always add this tag
        if let Some(tag) = annotation.make_re_tag() {
            record.push_aux(b"RE", rust_htslib::bam::record::Aux::Char(tag.as_bytes()[0]))?;

            // Count by region type
            match tag {
                "E" => exonic_count += 1,
                "N" => intronic_count += 1,
                "I" => intergenic_count += 1,
                _ => {}
            }
        }

        // Optionally add gene tags (GX, GN)
        if args.add_gene_tags {
            if let Some((gx_tag, gn_tag)) = annotation.make_gx_gn_tags() {
                record.push_aux(b"GX", rust_htslib::bam::record::Aux::String(&gx_tag))?;
                record.push_aux(b"GN", rust_htslib::bam::record::Aux::String(&gn_tag))?;
            }
        }

        // Optionally add transcript tag (TX)
        if args.add_tx_tag {
            if let Some(tx_tag) = annotation.make_tx_tag() {
                record.push_aux(b"TX", rust_htslib::bam::record::Aux::String(&tx_tag))?;
            }
        }

        // Optionally add antisense transcript tag (AN)
        if args.add_an_tag {
            if let Some(an_tag) = annotation.make_an_tag() {
                record.push_aux(b"AN", rust_htslib::bam::record::Aux::String(&an_tag))?;
            }
        }

        // Write the annotated record
        out.write(&record)?;
        
        count += 1;
        if count % 1_000_000 == 0 {
            println!("Processed {} million reads", count / 1_000_000);
        }
    }

    println!("Annotation complete!");
    println!("Total reads processed: {}", count);
    println!("Exonic reads: {} ({:.2}%)", exonic_count, 100.0 * exonic_count as f64 / count as f64);
    println!("Intronic reads: {} ({:.2}%)", intronic_count, 100.0 * intronic_count as f64 / count as f64);
    println!("Intergenic reads: {} ({:.2}%)", intergenic_count, 100.0 * intergenic_count as f64 / count as f64);

    Ok(())
}