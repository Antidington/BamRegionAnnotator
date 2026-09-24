# BamRegionAnnotator
* Add 10X standard `RE` attribute for BAM files

## Build

Use Rust 1.90.0 (pinned in `rust-toolchain.toml`), a C/C++ toolchain,
CMake, pkg-config, libclang, and Python 3. Dependencies and recovery provenance
are documented in [DEPENDENCIES.md](DEPENDENCIES.md).

```bash
cargo build --locked --release
```

## Usage
```bash
bam_region_annotator -h

Tool to annotate BAM files with exon/intron information

Usage: bam_region_annotator [OPTIONS] --reference <REFERENCE> --input <INPUT> --output <OUTPUT>

Options:
  -r, --reference <REFERENCE>
          Path to the reference directory containing STAR index files
  -i, --input <INPUT>
          Input BAM file
  -o, --output <OUTPUT>
          Output BAM file
      --strandedness <STRANDEDNESS>
          Chemistry strandedness (forward or reverse) [default: forward]
      --endedness <ENDEDNESS>
          Chemistry endedness (five_prime or three_prime) [default: three_prime]
      --region-min-overlap <REGION_MIN_OVERLAP>
          Minimum overlap fraction to consider a read exonic/intronic [default: 0.5]
      --include-exons
          Whether to include exons in annotation
      --include-introns
          Whether to include introns in annotation
      --intergenic-trim-bases <INTERGENIC_TRIM_BASES>
          Bases to trim from intergenic regions [default: 0]
      --intronic-trim-bases <INTRONIC_TRIM_BASES>
          Bases to trim from intronic regions [default: 0]
      --junction-trim-bases <JUNCTION_TRIM_BASES>
          Bases to trim from junction regions [default: 0]
      --add-gene-tags
          Add gene ID (GX) and gene name (GN) tags
      --add-tx-tag
          Add transcript (TX) tag
      --add-an-tag
          Add antisense transcript (AN) tag
  -h, --help
          Print help
  -V, --version
          Print version          
```
