use anyhow::{bail, ensure, Context, Result};
use clap::Parser;
use fastq_set::WhichEnd;
use rust_htslib::bam::{self, Read, Record};
use std::fs;
use std::io::ErrorKind;
use std::path::{Path, PathBuf};
use tx_annotation::transcript::{AnnotationParams, TranscriptAnnotator};

/// Tool to annotate BAM files with exon/intron information
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// Cell Ranger reference root (star/, genes/genes.gtf[.gz], reference.json)
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

    /// Minimum exonic overlap fraction (finite number from 0 to 1)
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
    validate_args(&args)?;
    ensure_output_absent(&args.output)?;

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

    check_eof(&mut bam)?;
    validate_reference(&args.reference, bam.header())?;

    // Create output BAM file with same header
    println!("Creating output BAM file: {}", args.output.display());
    let header = bam::Header::from_template(bam.header());
    let parent = args
        .output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let temporary = tempfile::Builder::new()
        .prefix(".bam-region-annotator-")
        .suffix(".bam")
        .tempfile_in(parent)?;
    let mut out = bam::Writer::from_path(temporary.path(), &header, bam::Format::Bam)
        .with_context(|| {
            format!(
                "Failed to create output BAM file: {}",
                args.output.display()
            )
        })?;

    // Process each read
    println!("Processing reads...");
    let mut count = 0_u64;
    let mut total_count = 0_u64;
    let mut unmapped_count = 0_u64;
    let mut exonic_count = 0_u64;
    let mut intronic_count = 0_u64;
    let mut intergenic_count = 0_u64;

    let input_header = bam::HeaderView::from_header(&header);
    for result in bam.records() {
        let mut record = result.with_context(|| "Failed to read BAM record")?;

        total_count += 1;

        // Skip unmapped reads
        if record.is_unmapped() {
            unmapped_count += 1;
            out.write(&record)?;
            continue;
        }

        validate_alignment(&record, &input_header)?;
        check_target_tags(&record, &args)?;

        // Annotate the read
        let annotation = annotator.annotate_alignment(&record);

        // Add region tag (RE) - always add this tag
        if let Some(tag) = annotation.make_re_tag() {
            record.push_aux(
                b"RE",
                rust_htslib::bam::record::Aux::Char(tag.as_bytes()[0]),
            )?;

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

    // Writer::drop does not expose hts_close errors. Reopen and fully decode
    // the completed temporary file before publishing it.
    drop(out);
    temporary
        .as_file()
        .sync_all()
        .context("Failed to sync output BAM")?;
    verify_output(temporary.path(), &header, total_count)?;
    temporary
        .persist_noclobber(&args.output)
        .map_err(|error| {
            let tempfile::PersistError { error, file } = error;
            drop(file);
            error
        })
        .with_context(|| {
            format!(
                "Failed to publish output without overwriting: {}",
                args.output.display()
            )
        })?;

    println!("Annotation complete!");
    println!("Total alignment records: {}", total_count);
    println!("Mapped alignment records: {}", count);
    println!("Unmapped alignment records: {}", unmapped_count);
    println!("Region percentages use mapped alignment records as denominator.");
    println!(
        "Exonic records: {} ({:.2}%)",
        exonic_count,
        percentage(exonic_count, count)
    );
    println!(
        "Intronic records: {} ({:.2}%)",
        intronic_count,
        percentage(intronic_count, count)
    );
    println!(
        "Intergenic records: {} ({:.2}%)",
        intergenic_count,
        percentage(intergenic_count, count)
    );
    Ok(())
}

fn validate_args(args: &Args) -> Result<()> {
    ensure!(
        args.region_min_overlap.is_finite() && (0.0..=1.0).contains(&args.region_min_overlap),
        "--region-min-overlap must be finite and between 0 and 1"
    );
    for (name, value) in [
        ("intergenic", args.intergenic_trim_bases),
        ("intronic", args.intronic_trim_bases),
        ("junction", args.junction_trim_bases),
    ] {
        ensure!(value >= 0, "--{name}-trim-bases must be nonnegative");
    }
    Ok(())
}

fn ensure_output_absent(output: &Path) -> Result<()> {
    // symlink_metadata also rejects dangling symlinks. Existing hardlinks and
    // symlinks to the input are covered by this no-overwrite policy.
    match fs::symlink_metadata(output) {
        Ok(_) => bail!(
            "Output already exists (refusing to overwrite): {}",
            output.display()
        ),
        Err(error) if error.kind() == ErrorKind::NotFound => Ok(()),
        Err(error) => Err(error).context("Cannot inspect output path"),
    }
}

fn check_eof(reader: &mut bam::Reader) -> Result<()> {
    // SAFETY: Reader owns a live htsFile throughout this call; mutable borrowing
    // prevents concurrent use. HTSlib performs a seek and restores the position.
    let status = unsafe { rust_htslib::htslib::hts_check_EOF(reader.htsfile()) };
    ensure!(
        status == 1 || status == 3,
        "BAM EOF integrity check failed (status {status}); input may be truncated or unseekable"
    );
    Ok(())
}

fn read_lines(path: &Path) -> Result<Vec<String>> {
    Ok(fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?
        .lines()
        .map(str::to_owned)
        .collect())
}

fn validate_reference(reference: &Path, header: &bam::HeaderView) -> Result<()> {
    let names = read_lines(&reference.join("star/chrName.txt"))?;
    let lengths = read_lines(&reference.join("star/chrLength.txt"))?;
    let starts = read_lines(&reference.join("star/chrStart.txt"))?;
    ensure!(
        !names.is_empty()
            && names.len() == lengths.len()
            && names.len() == header.target_count() as usize,
        "BAM/reference contig count mismatch or empty reference"
    );
    ensure!(
        starts.len() == names.len() + 1,
        "Invalid STAR chrStart.txt: expected contig starts and terminal offset"
    );
    let starts: Vec<i64> = starts
        .iter()
        .map(|x| x.parse())
        .collect::<std::result::Result<_, _>>()
        .context("Invalid STAR chromosome offset")?;
    ensure!(starts[0] == 0, "STAR chromosome offsets must start at zero");
    let mut seen = std::collections::HashSet::new();
    for (tid, (name, length)) in names.iter().zip(lengths).enumerate() {
        let length: u64 = length.parse().context("Invalid STAR chromosome length")?;
        ensure!(
            seen.insert(name) && !name.is_empty(),
            "Duplicate or empty reference contig name"
        );
        ensure!(
            length > 0 && length <= i64::MAX as u64,
            "Invalid reference contig length: {name}"
        );
        ensure!(header.tid2name(tid as u32) == name.as_bytes() && header.target_len(tid as u32) == Some(length),
            "BAM/reference contig mismatch at tid {tid}: expected {name}, length {length}; names, order and lengths must match");
        ensure!(
            starts[tid] >= 0
                && starts[tid]
                    .checked_add(length as i64)
                    .is_some_and(|end| end <= starts[tid + 1]),
            "Invalid STAR chromosome offsets for {name}"
        );
    }
    Ok(())
}

fn validate_alignment(record: &Record, header: &bam::HeaderView) -> Result<()> {
    let tid = record.tid();
    let name = String::from_utf8_lossy(record.qname());
    ensure!(
        tid >= 0 && (tid as u32) < header.target_count(),
        "Invalid reference id for record {name}"
    );
    let end = record.cigar().end_pos();
    ensure!(
        record.pos() >= 0
            && end > record.pos()
            && end as u64 <= header.target_len(tid as u32).unwrap(),
        "Alignment outside reference bounds for record {name}"
    );
    Ok(())
}

fn check_target_tags(record: &Record, args: &Args) -> Result<()> {
    for entry in record.aux_iter() {
        let (tag, _) = entry.context("Invalid auxiliary tag")?;
        let target = tag == b"RE"
            || (args.add_gene_tags && (tag == b"GX" || tag == b"GN"))
            || (args.add_tx_tag && tag == b"TX")
            || (args.add_an_tag && tag == b"AN");
        ensure!(
            !target,
            "Record {} already has target tag {}",
            String::from_utf8_lossy(record.qname()),
            String::from_utf8_lossy(tag)
        );
    }
    Ok(())
}

fn verify_output(path: &Path, header: &bam::Header, expected: u64) -> Result<()> {
    let mut reader = bam::Reader::from_path(path).context("Cannot reopen temporary BAM")?;
    check_eof(&mut reader)?;
    ensure!(
        bam::Header::from_template(reader.header()).to_bytes() == header.to_bytes(),
        "Output BAM header changed"
    );
    let mut count = 0_u64;
    for result in reader.records() {
        result.context("Output BAM failed full decode")?;
        count += 1;
    }
    ensure!(
        count == expected,
        "Output BAM record count mismatch: expected {expected}, found {count}"
    );
    Ok(())
}

fn percentage(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        100.0 * part as f64 / total as f64
    }
}
