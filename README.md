# BamRegionAnnotator

Annotate alignment records with Cell Ranger-style region tags (`RE:A:E`,
`RE:A:N`, `RE:A:I`) and optional gene/transcript tags. The existing vendored
transcript annotation algorithm is retained.

## Build and test

Rust 1.90.0 is pinned in `rust-toolchain.toml`. Install a C/C++ toolchain,
CMake, pkg-config, libclang and Python 3. On macOS, Xcode Command Line Tools
and Homebrew `cmake`, `pkg-config`, `llvm` and `python@3.11` provide these tools.
If necessary, set `LIBCLANG_PATH` to the LLVM library directory and
`PYO3_PYTHON` to the Python executable.

```bash
cargo build --locked --release
cargo test --locked
cargo fmt --check
cargo clippy --locked --all-targets --no-deps -- -D warnings
```

Dependency versions are locked; do not run `cargo update` as an installation
step. The recovered dependency source, local compatibility patch and license
text are described in [DEPENDENCIES.md](DEPENDENCIES.md).

Linux x86_64 is the formal target; macOS ARM is the development target. The
local container recipe pins the Rust base image and installs native build
packages (the Debian package repository is not snapshotted):

```bash
docker build --platform linux/amd64 -t bam-region-annotator-validation:rust1.90 - < Dockerfile.validation
mkdir -p target/linux-amd64
docker run --rm --platform linux/amd64 \
  -e RUSTUP_TOOLCHAIN=1.90.0 -e CARGO_TARGET_DIR=/work/target/linux-amd64 \
  -v "$PWD:/work:ro" \
  -v "$PWD/target/linux-amd64:/work/target/linux-amd64" \
  bam-region-annotator-validation:rust1.90 cargo test --locked
```

## Reference and input contract

`--reference` is a Cell Ranger-style reference root, containing:

```text
reference.json
star/chrName.txt
star/chrLength.txt
star/chrStart.txt
star/transcriptInfo.tab
star/exonInfo.tab
genes/genes.gtf             # genes.gtf.gz is also accepted
```

A standalone STAR index is insufficient. BAM `@SQ` names, order and lengths
must exactly match this reference; subsets and reordered headers are rejected,
not remapped. This dictionary check does not prove reference sequence identity.
The STAR tables and GTF must come from the same reference build.

Input BAM must be seekable and have a valid EOF marker. Every alignment record
is retained in order. Secondary and supplementary alignments are annotated;
there is no deduplication, primary-only filtering or mate-level aggregation.
Unmapped records, including their existing tags, are written unchanged.

## Usage

```bash
target/release/bam_region_annotator \
  --reference /path/to/reference \
  --input input.bam --output annotated.bam

# Also generate gene, transcript and antisense tags; include intronic assignments.
target/release/bam_region_annotator \
  -r /path/to/reference -i input.bam -o annotated_with_genes.bam \
  --include-introns --add-gene-tags --add-tx-tag --add-an-tag
```

| Option | Default and behavior |
| --- | --- |
| `--strandedness` | `forward`; also accepts `reverse`, case insensitive |
| `--endedness` | `three_prime`; also accepts `five_prime`, case insensitive |
| `--region-min-overlap` | `0.5`; finite value in `[0,1]`, used for exonic classification |
| `--include-exons` | Enabled by default; retained as a compatibility flag, with no disabling option |
| `--include-introns` | Disabled by default; enables intronic gene/transcript assignments, not intronic RE classification |
| `--intergenic-trim-bases`, `--intronic-trim-bases`, `--junction-trim-bases` | `0`; nonnegative integers |
| `--add-gene-tags` | Add `GX` and `GN` when an assignment exists |
| `--add-tx-tag` | Add `TX` when an assignment exists |
| `--add-an-tag` | Add `AN` when an antisense assignment exists |

Region classification and gene/transcript assignment are distinct. Intronic
classification uses the existing algorithm's full-transcript-span overlap
rule; `--region-min-overlap` is not an intronic threshold. In intronic assignment
mode, TX/AN can contain gene ID and strand rather than transcript coordinates.

## Conflicts, output and statistics

For mapped records, an existing `RE` tag is an error. If a tag-writing option
is enabled, any existing corresponding tag (`GX`/`GN`, `TX`, or `AN`) is also an
error, even when the new annotation would have no value. There is no overwrite
mode. Optional tags whose options are disabled are preserved.

The output path must not exist, including symlinks and hardlinks. This also
prevents overwriting the input. Output is first written to a temporary file in
the destination directory. After the writer closes, the file is synced, its EOF
and header checked, and every record decoded and counted before publication.
Publication refuses to replace a file created concurrently. Reported errors
remove the temporary file; a forced process kill or machine failure can leave a
`.bam-region-annotator-*.bam` temporary file. No BAM index is generated.

Validation adds a complete output read pass and requires space for the output
in its destination directory. A completed output is published only after this
pass succeeds.

Statistics count **alignment records**, not unique reads, fragments or molecules.
Total includes mapped and unmapped records. Region percentages use mapped
records as their denominator and report `0.00%` when there are none. The success
message is emitted only after publication.

## Validation scope

The [synthetic fixtures](tests/fixtures/README.md) capture the original wrapper's
results using the recovered dependency. Tests compare complete decoded records
and headers across four parameter configurations, and exercise invalid inputs,
tag conflicts, output preservation and cleanup. They do not establish identity
with the unavailable historical dependency source or with a particular
Cell Ranger release.

## Code layout

`src/main.rs` calls the pipeline. `cli.rs` owns argument parsing and annotation
parameters; `reference.rs` checks reference dictionaries and alignment bounds;
`tags.rs` owns tag conflicts and annotation writes; `pipeline.rs` streams records;
`stats.rs` defines counters and reporting; `bam_io.rs` validates and publishes
temporary output. Upstream annotation algorithms remain under `lib/rust`.

The obsolete, uncompiled `process_bam.rs` has been removed; use the Cargo binary
and documented CLI as the single application entry point.
