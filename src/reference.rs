use anyhow::{ensure, Context, Result};
use rust_htslib::bam::{self, Record};
use std::{fs, path::Path};

fn read_lines(path: &Path) -> Result<Vec<String>> {
    Ok(fs::read_to_string(path)
        .with_context(|| format!("Cannot read {}", path.display()))?
        .lines()
        .map(str::to_owned)
        .collect())
}

pub(crate) fn validate_reference(reference: &Path, header: &bam::HeaderView) -> Result<()> {
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

pub(crate) fn validate_alignment(record: &Record, header: &bam::HeaderView) -> Result<()> {
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
