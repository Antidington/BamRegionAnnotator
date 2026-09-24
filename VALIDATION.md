# Refactor validation

Last verified: 2026-09-24.
Application source commit: `a45760de409eb2c6339d7c74e6e645c767c0b1f9`.

## Results

| Check | Result |
| --- | --- |
| macOS ARM, Rust 1.90.0, `cargo test --locked --offline` | Passed: 3 unit tests, 1 baseline test with 4 configurations, 10 reliability tests |
| Application `cargo clippy --locked --offline --all-targets --no-deps -- -D warnings` | Passed |
| `cargo fmt --check` | Passed |
| README and fixture documentation local links | Passed |
| Recovered fastq_set comparison against pinned upstream revision | Only the documented redundant-parentheses patch differs |
| Linux x86_64 container, Rust 1.90.0, `cargo test --locked` | Passed: the same 14 tests, including all 4 original-result configurations |
| macOS ARM, `cargo build --locked --offline --release` | Passed; optimized executable built |
| macOS release executable end-to-end smoke | Passed: all optional tags, samtools 1.21 quickcheck and full decoded SAM comparison with original all-tags snapshot |

The existing Homebrew rustc (1.87.0) failed to load its LLVM library. Validation
used the already installed Rust 1.90.0 toolchain binaries directly; no system
toolchain or user configuration was changed. Dependency cache was isolated at
`/private/tmp/bra-cargo-home`. Build products are ignored under `target/`.

Linux validation uses the pinned base in `Dockerfile.validation`, running
`linux/amd64` on the local ARM Docker engine. It is actual x86_64 container
execution under emulation, not validation on a remote physical x86_64 host.
Source is mounted read-only and Linux build products use `target/linux-amd64`.

macOS release artifact: `target/release/bam_region_annotator` (Mach-O arm64).
SHA-256: `9952343fd2556f80330e1dcc996a3261f30679367db33635c5c953ed8c16d5dc`.

## What the checks establish

- Complete decoded record and header equality with snapshots from the original
  application using the recovered dependency (see [fixture provenance](tests/fixtures/README.md)).
- Existing output preservation, input hardlink/symlink protection, concurrent
  destination creation protection and temporary-output cleanup on errors.
- Target-tag conflict handling and preservation of unselected tags and unmapped
  records; secondary/supplementary records remain included.
- Contig mismatch rejection, parameter validation, alignment bounds, empty-input
  statistics, truncated/corrupt BAM rejection and validation before publication.

The existing upstream `barcode`/PyO3 macro emits a `non_local_definitions`
warning. It is outside the application code and was not suppressed or rewritten.
The deliberately malformed BAM tests can emit HTSlib diagnostics while passing.

## Limits

- Historical fastq_set source identity remains unknown; see [dependency provenance](DEPENDENCIES.md).
- These synthetic cases are not independent validation against a released
  Cell Ranger executable, nor acceptance on a production cohort.
- No throughput or large-BAM memory benchmark was performed. Publication now
  includes a full decode pass, adding output I/O by design.
- Vendored annotation algorithms and their disabled upstream test targets were
  retained; application tests supply the regression protection for this change.
- The base image and Rust dependencies are pinned; Debian build-package versions
  are not pinned to a dated package snapshot.
