use crate::tags::TagOptions;
use anyhow::{ensure, Result};
use clap::Parser;
use fastq_set::WhichEnd;
use std::path::PathBuf;
use tx_annotation::transcript::AnnotationParams;

/// Tool to annotate BAM files with exon/intron information
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub(crate) struct Args {
    /// Cell Ranger reference root (star/, genes/genes.gtf[.gz], reference.json)
    #[arg(short, long)]
    pub(crate) reference: PathBuf,

    /// Input BAM file
    #[arg(short, long)]
    pub(crate) input: PathBuf,

    /// Output BAM file
    #[arg(short, long)]
    pub(crate) output: PathBuf,

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

impl Args {
    pub(crate) fn validate(&self) -> Result<()> {
        ensure!(
            self.region_min_overlap.is_finite() && (0.0..=1.0).contains(&self.region_min_overlap),
            "--region-min-overlap must be finite and between 0 and 1"
        );
        for (name, value) in [
            ("intergenic", self.intergenic_trim_bases),
            ("intronic", self.intronic_trim_bases),
            ("junction", self.junction_trim_bases),
        ] {
            ensure!(value >= 0, "--{name}-trim-bases must be nonnegative");
        }
        Ok(())
    }

    pub(crate) fn annotation_params(&self) -> Result<AnnotationParams> {
        // Parse chemistry strandedness
        let chemistry_strandedness = match self.strandedness.to_lowercase().as_str() {
            "forward" => cr_types::ReqStrand::Forward,
            "reverse" => cr_types::ReqStrand::Reverse,
            _ => anyhow::bail!("Invalid strandedness: {}", self.strandedness),
        };

        // Parse chemistry endedness
        let chemistry_endedness = match self.endedness.to_lowercase().as_str() {
            "five_prime" => WhichEnd::FivePrime,
            "three_prime" => WhichEnd::ThreePrime,
            _ => anyhow::bail!("Invalid endedness: {}", self.endedness),
        };

        // Create annotation parameters
        Ok(AnnotationParams {
            chemistry_strandedness,
            chemistry_endedness,
            intergenic_trim_bases: self.intergenic_trim_bases,
            intronic_trim_bases: self.intronic_trim_bases,
            junction_trim_bases: self.junction_trim_bases,
            region_min_overlap: self.region_min_overlap,
            include_exons: self.include_exons,
            include_introns: self.include_introns,
        })
    }

    pub(crate) fn tag_options(&self) -> TagOptions {
        TagOptions {
            genes: self.add_gene_tags,
            transcripts: self.add_tx_tag,
            antisense: self.add_an_tag,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn legacy_defaults_and_case_insensitive_chemistry_are_preserved() {
        let base = ["annotator", "-r", "ref", "-i", "input", "-o", "output"];
        let args = Args::try_parse_from(base).unwrap();
        let params = args.annotation_params().unwrap();
        assert_eq!(params.chemistry_strandedness, cr_types::ReqStrand::Forward);
        assert!(matches!(params.chemistry_endedness, WhichEnd::ThreePrime));
        assert!(params.include_exons);
        assert!(!params.include_introns);
        assert_eq!(params.region_min_overlap, 0.5);
        assert_eq!(
            (
                params.intergenic_trim_bases,
                params.intronic_trim_bases,
                params.junction_trim_bases
            ),
            (0, 0, 0)
        );
        let options = args.tag_options();
        assert!(!options.genes && !options.transcripts && !options.antisense);
        let args = Args::try_parse_from(base.into_iter().chain([
            "--strandedness",
            "REVERSE",
            "--endedness",
            "FIVE_PRIME",
            "--include-exons",
        ]))
        .unwrap();
        let params = args.annotation_params().unwrap();
        assert_eq!(params.chemistry_strandedness, cr_types::ReqStrand::Reverse);
        assert!(matches!(params.chemistry_endedness, WhichEnd::FivePrime));
        assert!(params.include_exons);
    }
}
