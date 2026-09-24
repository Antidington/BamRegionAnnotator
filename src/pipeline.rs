use crate::{
    bam_io,
    cli::Args,
    reference,
    stats::Stats,
    tags::{self, TagOptions},
};
use anyhow::{Context, Result};
use rust_htslib::bam::{self, Read};
use tx_annotation::transcript::TranscriptAnnotator;

pub(crate) fn run(args: Args) -> Result<()> {
    args.validate()?;
    bam_io::ensure_output_absent(&args.output)?;
    let params = args.annotation_params()?;
    println!("Initializing transcript annotator...");
    let annotator = TranscriptAnnotator::new(&args.reference, params)?;
    println!("Opening input BAM file: {}", args.input.display());
    let mut reader = bam::Reader::from_path(&args.input)
        .with_context(|| format!("Failed to open BAM file: {}", args.input.display()))?;
    bam_io::check_eof(&mut reader)?;
    reference::validate_reference(&args.reference, reader.header())?;

    println!("Creating output BAM file: {}", args.output.display());
    let header = bam::Header::from_template(reader.header());
    let temporary = bam_io::temporary_output(&args.output)?;
    let mut writer = bam::Writer::from_path(temporary.path(), &header, bam::Format::Bam)
        .with_context(|| {
            format!(
                "Failed to create output BAM file: {}",
                args.output.display()
            )
        })?;
    println!("Processing reads...");
    let stats = process_records(&mut reader, &mut writer, &annotator, args.tag_options())?;
    drop(writer);
    bam_io::publish(temporary, &args.output, &header, stats.total)?;
    stats.report();
    Ok(())
}

fn process_records(
    reader: &mut bam::Reader,
    writer: &mut bam::Writer,
    annotator: &TranscriptAnnotator,
    options: TagOptions,
) -> Result<Stats> {
    let header = bam::HeaderView::from_header(&bam::Header::from_template(reader.header()));
    let mut stats = Stats::default();
    for result in reader.records() {
        let mut record = result.context("Failed to read BAM record")?;
        if record.is_unmapped() {
            writer.write(&record)?;
            stats.observe_unmapped();
            continue;
        }
        reference::validate_alignment(&record, &header)?;
        let region = tags::annotate(&mut record, annotator, options)?;
        writer.write(&record)?;
        stats.observe_region(region);
        if stats.mapped % 1_000_000 == 0 {
            println!("Processed {} million reads", stats.mapped / 1_000_000);
        }
    }
    Ok(stats)
}
