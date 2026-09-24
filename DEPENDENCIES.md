# Dependency provenance

The existing `lib/rust` source tree is retained from repository commit
`2bbbaf75c06e4a921d0d06a54968af9f383e4607`. Its original upstream revision is not recorded.

## Recovered fastq_set

- Repository: https://github.com/10XGenomics/fastq_set
- Revision: `8896646805641d824890362e762bd68652f0271c`
- Package version: `0.5.3`, matching the original Cargo.lock entry.
- Files: Cargo.toml, src/, benches/, LICENSE.txt, README.md, CHANGELOG.md.
- Local compatibility patch: remove redundant closure parentheses in
  src/filenames/bcl2fastq.rs; upstream denies warnings and Rust 1.90 rejects them.
  No algorithm change.
- License text: retained verbatim in lib/rust/fastq_set/LICENSE.txt. Although
  the upstream manifest says MIT, the bundled text additionally restricts use
  to a 10x Genomics product or data generated with such a product (condition 2).
  Do not describe this source as unrestricted MIT.

The original local fastq_set source was absent and its lockfile entry has no
revision or checksum. This recovery is pinned but is not proof of historical
source identity. Regression comparisons use the original application with this
recovered dependency; they do not establish equivalence with an unavailable
historical executable. Other external dependencies are fixed by Cargo.lock;
use `--locked` for builds and tests.
