use anyhow::{bail, ensure, Context, Result};
use rust_htslib::bam::{self, Read};
use std::{fs, io::ErrorKind, path::Path};
use tempfile::NamedTempFile;

pub(crate) fn ensure_output_absent(output: &Path) -> Result<()> {
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

pub(crate) fn check_eof(reader: &mut bam::Reader) -> Result<()> {
    // SAFETY: Reader owns a live htsFile throughout this call; mutable borrowing
    // prevents concurrent use. HTSlib performs a seek and restores the position.
    let status = unsafe { rust_htslib::htslib::hts_check_EOF(reader.htsfile()) };
    ensure!(
        status == 1 || status == 3,
        "BAM EOF integrity check failed (status {status}); input may be truncated or unseekable"
    );
    Ok(())
}

pub(crate) fn temporary_output(destination: &Path) -> Result<NamedTempFile> {
    let parent = destination
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    tempfile::Builder::new()
        .prefix(".bam-region-annotator-")
        .suffix(".bam")
        .tempfile_in(parent)
        .context("Cannot create temporary output BAM")
}

// Caller must drop the BAM writer before publication. rust-htslib's Writer
// destructor discards hts_close errors, so verify the complete closed file.
pub(crate) fn publish(
    temporary: NamedTempFile,
    destination: &Path,
    header: &bam::Header,
    expected: u64,
) -> Result<()> {
    temporary
        .as_file()
        .sync_all()
        .context("Failed to sync output BAM")?;
    verify_output(temporary.path(), header, expected)?;
    temporary
        .persist_noclobber(destination)
        .map_err(|error| {
            let tempfile::PersistError { error, file } = error;
            drop(file);
            error
        })
        .with_context(|| {
            format!(
                "Failed to publish output without overwriting: {}",
                destination.display()
            )
        })?;
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

#[cfg(test)]
mod tests {
    use super::*;

    fn complete_empty_bam(temporary: &NamedTempFile, header: &bam::Header) {
        let writer = bam::Writer::from_path(temporary.path(), header, bam::Format::Bam).unwrap();
        drop(writer);
    }

    #[test]
    fn publication_race_preserves_competing_file_and_removes_temporary() {
        let dir = tempfile::tempdir().unwrap();
        let destination = dir.path().join("output.bam");
        ensure_output_absent(&destination).unwrap();
        let temporary = temporary_output(&destination).unwrap();
        let path = temporary.path().to_path_buf();
        let header = bam::Header::new();
        complete_empty_bam(&temporary, &header);
        fs::write(&destination, b"created by another process").unwrap();
        assert!(publish(temporary, &destination, &header, 0).is_err());
        assert_eq!(
            fs::read(destination).unwrap(),
            b"created by another process"
        );
        assert!(!path.exists());
    }

    #[test]
    fn failed_output_validation_prevents_publication_and_cleans_up() {
        let dir = tempfile::tempdir().unwrap();
        let header = bam::Header::new();
        for failure in ["truncated", "count", "header", "corrupt"] {
            let destination = dir.path().join(format!("{failure}.bam"));
            let temporary = temporary_output(&destination).unwrap();
            let path = temporary.path().to_path_buf();
            complete_empty_bam(&temporary, &header);
            let mut expected_header = header.clone();
            let mut expected_count = 0;
            match failure {
                "truncated" => {
                    let size = temporary.as_file().metadata().unwrap().len();
                    temporary.as_file().set_len(size - 28).unwrap();
                }
                "count" => expected_count = 1,
                "header" => {
                    expected_header.push_comment(b"unexpected");
                }
                "corrupt" => fs::write(&path, b"not a BAM").unwrap(),
                _ => unreachable!(),
            }
            assert!(
                publish(temporary, &destination, &expected_header, expected_count).is_err(),
                "{failure}"
            );
            assert!(!destination.exists());
            assert!(!path.exists());
        }
    }
}
