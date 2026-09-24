# BamRegionAnnotator

A Rust CLI that annotates mapped BAM alignment records as **exonic**, **intronic**
or **intergenic** using the vendored Cell Ranger transcript annotator. It adds
`RE` region tags and can also write gene (`GX`, `GN`), transcript (`TX`) and
antisense (`AN`) tags.

## Requirements

- Rust **1.90.0**, pinned in [rust-toolchain.toml](rust-toolchain.toml), with a
  rustup-managed Cargo toolchain.
- A C/C++ toolchain, CMake, pkg-config, libclang and Python 3 for native builds.
- A BAM file and the matching Cell Ranger-style reference described below.

On macOS, Xcode Command Line Tools and Homebrew `cmake`, `pkg-config`, `llvm` and
`python@3.11` provide the native build tools. If discovery fails, set
`LIBCLANG_PATH` to the LLVM library directory and `PYO3_PYTHON` to the Python
executable. A standalone Homebrew Cargo does not select the pinned Rust version.

## Installation

From the repository root:

```bash
cargo build --locked --release
```

The executable is `target/release/bam_region_annotator`. Use the committed
`Cargo.lock`; `cargo update` is not an installation step.

## Configuration: reference and input

Pass a Cell Ranger-style reference root to `--reference`:

```text
reference.json
star/chrName.txt
star/chrLength.txt
star/chrStart.txt
star/transcriptInfo.tab
star/exonInfo.tab
genes/genes.gtf             # genes.gtf.gz is also accepted
```

A standalone STAR index is insufficient. BAM `@SQ` names, order and lengths must
exactly match the reference; subsets and reordered headers are rejected. The
STAR tables and GTF must belong to the same reference build. Dictionary matching
does not establish reference sequence identity.

Input BAM must be seekable and have a valid EOF marker. Choose an **output path
that does not exist**, in an existing directory with space for the output BAM.

## Quick start

Replace the reference and input paths with your own:

```bash
target/release/bam_region_annotator \
  --reference /path/to/reference \
  --input input.bam --output annotated.bam
```

To include intronic gene/transcript assignments and all optional tags:

```bash
target/release/bam_region_annotator \
  -r /path/to/reference -i input.bam -o annotated_with_genes.bam \
  --include-introns --add-gene-tags --add-tx-tag --add-an-tag
```

Full CLI help:

```bash
target/release/bam_region_annotator --help
```

## Annotation behavior

| Region | BAM tag |
| --- | --- |
| Exonic | `RE:A:E` |
| Intronic | `RE:A:N` |
| Intergenic | `RE:A:I` |

Every alignment record is retained in input order. Mapped secondary and
supplementary alignments are annotated; unmapped records and their tags are
preserved unchanged. There is no deduplication or mate-level aggregation.

| Option | Default and behavior |
| --- | --- |
| `--strandedness` | `forward`; also accepts `reverse`, case insensitive |
| `--endedness` | `three_prime`; also accepts `five_prime`, case insensitive |
| `--region-min-overlap` | `0.5`; finite value in `[0,1]` for exonic classification |
| `--include-exons` | Enabled; compatibility flag with no disabling option |
| `--include-introns` | Disabled; enables intronic gene/transcript assignments, not intronic RE classification |
| `--intergenic-trim-bases`, `--intronic-trim-bases`, `--junction-trim-bases` | `0`; nonnegative integers |
| `--add-gene-tags` | Disabled; add `GX` and `GN` when an assignment exists |
| `--add-tx-tag` | Disabled; add `TX` when an assignment exists |
| `--add-an-tag` | Disabled; add `AN` when an antisense assignment exists |

Region classification and gene/transcript assignment are distinct. Intronic
classification requires full overlap of the alignment span with a transcript;
`--region-min-overlap` is not an intronic threshold. With intronic assignments,
TX/AN can contain a gene ID and strand instead of transcript coordinates.

For mapped records, an existing `RE` is an error. Enabling a tag-writing option
also makes any existing corresponding tag an error, even if the new annotation
would have no value. Optional tags whose options are disabled are preserved.
There is no tag-overwrite mode.

## Output and statistics

Output is written to a temporary file in the destination directory. After the
writer closes, the file is synced, its EOF and header checked, and every record
decoded and counted before publication. This adds a complete output read pass.
Publication refuses to replace existing files, including symlinks, hardlinks and
files created concurrently. No BAM index is generated.

Reported errors remove the temporary file. A forced process kill or machine
failure can leave a `.bam-region-annotator-*.bam` file. The success message is
printed only after the completed output is published.

Statistics count **alignment records**, not unique reads, fragments or molecules.
Total includes mapped and unmapped records; region percentages use mapped
records as the denominator and report `0.00%` when there are none.

## Documentation

- [Validation results and limits](VALIDATION.md): macOS ARM and Linux x86_64
  container checks; no production-cohort or large-BAM performance acceptance.
- [Dependency provenance and license text](DEPENDENCIES.md): pinned recovery
  source and the limits of historical source equivalence.
- [Regression fixtures](tests/fixtures/README.md): synthetic inputs and original
  wrapper snapshots, not independent validation against a Cell Ranger release.

## Development and testing

```bash
cargo test --locked
cargo fmt --check
cargo clippy --locked --all-targets --no-deps -- -D warnings
```

The [pipeline](src/pipeline.rs) connects CLI configuration, reference checks,
tag writing and statistics. [BAM I/O](src/bam_io.rs) handles output validation
and publication. The upstream annotation algorithm remains under `lib/rust`.

### Linux x86_64 validation

[Dockerfile.validation](Dockerfile.validation) pins the Rust base image; Debian
build packages are not pinned to a dated package snapshot. From the repository
root, with Docker available:

```bash
docker build --platform linux/amd64 -t bam-region-annotator-validation:rust1.90 - < Dockerfile.validation
mkdir -p target/linux-amd64
docker run --rm --platform linux/amd64 \
  -e RUSTUP_TOOLCHAIN=1.90.0 -e CARGO_TARGET_DIR=/work/target/linux-amd64 \
  -v "$PWD:/work:ro" \
  -v "$PWD/target/linux-amd64:/work/target/linux-amd64" \
  bam-region-annotator-validation:rust1.90 cargo test --locked
```
