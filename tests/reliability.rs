mod common;

use common::*;
use rust_htslib::bam::{self, record::Aux};
use std::fs;
use std::path::Path;

fn assert_failure(result: &std::process::Output, message: &str) {
    let stderr = String::from_utf8_lossy(&result.stderr);
    assert!(!result.status.success(), "unexpected success");
    assert!(stderr.contains(message), "expected {message:?}: {stderr}");
}

fn assert_no_temporary_files(dir: &Path) {
    for entry in fs::read_dir(dir).unwrap() {
        assert!(!entry
            .unwrap()
            .file_name()
            .to_string_lossy()
            .starts_with(".bam-region-annotator-"));
    }
}

#[test]
fn existing_output_and_input_aliases_are_never_overwritten() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let original = fs::read(&input).unwrap();
    let output = dir.path().join("output.bam");
    fs::write(&output, b"keep me").unwrap();
    assert_failure(&run(&input, &output, &[]), "Output already exists");
    assert_eq!(fs::read(&output).unwrap(), b"keep me");
    assert_failure(&run(&input, &input, &[]), "Output already exists");
    let hardlink = dir.path().join("hardlink.bam");
    fs::hard_link(&input, &hardlink).unwrap();
    assert_failure(&run(&input, &hardlink, &[]), "Output already exists");
    #[cfg(unix)]
    for (name, target) in [
        ("symlink.bam", input.clone()),
        ("dangling.bam", dir.path().join("missing")),
    ] {
        let link = dir.path().join(name);
        std::os::unix::fs::symlink(target, &link).unwrap();
        assert_failure(&run(&input, &link, &[]), "Output already exists");
    }
    assert_eq!(fs::read(input).unwrap(), original);
    assert_no_temporary_files(dir.path());
}

#[test]
fn contig_order_names_lengths_and_subsets_are_checked() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    let cases = [
        vec![("chr2", 1000), ("chr1", 1000)],
        vec![("chrX", 1000), ("chr2", 1000)],
        vec![("chr1", 999), ("chr2", 1000)],
        vec![("chr1", 1000)],
        vec![("chr1", 1000), ("chr2", 1000), ("chr3", 1000)],
    ];
    for contigs in cases {
        let mut header = bam::Header::new();
        for (name, length) in contigs {
            header.push_record(
                bam::header::HeaderRecord::new(b"SQ")
                    .push_tag(b"SN", name)
                    .push_tag(b"LN", length),
            );
        }
        write_bam(&input, &header, &[]);
        assert_failure(&run(&input, &output, &[]), "BAM/reference contig");
        assert!(!output.exists());
        assert_no_temporary_files(dir.path());
    }
}

#[test]
fn all_selected_tag_conflicts_abort_and_remove_partial_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    input_bam(&input);
    let (header, records) = read_bam(&input);
    for (tag, flags) in [
        (b"RE", vec![]),
        (b"GX", vec!["--add-gene-tags"]),
        (b"GN", vec!["--add-gene-tags"]),
        (b"TX", vec!["--add-tx-tag"]),
        (b"AN", vec!["--add-an-tag"]),
    ] {
        let mut changed = records.clone();
        // Place the conflict late so earlier records have already been written.
        changed[9].push_aux(tag, Aux::String("existing")).unwrap();
        write_bam(&input, &header, &changed);
        assert_failure(&run(&input, &output, &flags), "already has target tag");
        assert!(!output.exists());
        assert_no_temporary_files(dir.path());
    }
    // Also reject selected tags even if annotation would not produce a new value.
    let mut changed = records.clone();
    changed[2].push_aux(b"GX", Aux::String("stale")).unwrap();
    write_bam(&input, &header, &changed);
    assert_failure(
        &run(&input, &output, &["--add-gene-tags"]),
        "already has target tag GX",
    );
    assert!(!output.exists());
    assert_no_temporary_files(dir.path());
}

#[test]
fn unselected_tags_and_unmapped_records_are_preserved() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    input_bam(&input);
    let (header, mut records) = read_bam(&input);
    for tag in [b"GX", b"GN", b"TX", b"AN"] {
        records[0].push_aux(tag, Aux::String("keep")).unwrap();
        records[10].push_aux(tag, Aux::String("keep")).unwrap();
    }
    records[10].push_aux(b"RE", Aux::Char(b'N')).unwrap();
    write_bam(&input, &header, &records);
    assert_success(&run(&input, &output, &[]));
    let (_, actual) = read_bam(&output);
    for tag in [b"GX", b"GN", b"TX", b"AN"] {
        assert_eq!(actual[0].aux(tag).unwrap(), Aux::String("keep"));
    }
    assert_eq!(actual[10], records[10]);
    assert_no_temporary_files(dir.path());
}

#[test]
fn empty_and_unmapped_only_inputs_have_finite_zero_percentages() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let (header, records) = read_bam(&input);
    for (name, subset) in [("empty", &records[0..0]), ("unmapped", &records[10..])] {
        write_bam(&input, &header, subset);
        let output = dir.path().join(format!("{name}.bam"));
        let result = run(
            &input,
            &output,
            &["--add-gene-tags", "--add-tx-tag", "--add-an-tag"],
        );
        assert_success(&result);
        let stdout = String::from_utf8(result.stdout).unwrap();
        assert!(stdout.contains("Mapped alignment records: 0"));
        assert!(stdout.contains("Exonic records: 0 (0.00%)"));
        assert!(!stdout.contains("NaN"));
        assert_eq!(read_bam(&output).1, subset);
    }
}

#[test]
fn statistics_count_alignment_records_including_secondary_and_supplementary() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let result = run(&input, &dir.path().join("output.bam"), &[]);
    assert_success(&result);
    let stdout = String::from_utf8(result.stdout).unwrap();
    for expected in [
        "Total alignment records: 11",
        "Mapped alignment records: 10",
        "Unmapped alignment records: 1",
        "Exonic records: 8 (80.00%)",
        "Intronic records: 1 (10.00%)",
        "Intergenic records: 1 (10.00%)",
    ] {
        assert!(stdout.contains(expected), "{stdout}");
    }
}

#[test]
fn invalid_parameters_do_not_create_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let output = dir.path().join("output.bam");
    for args in [
        vec!["--region-min-overlap=NaN"],
        vec!["--region-min-overlap=inf"],
        vec!["--region-min-overlap=-0.1"],
        vec!["--region-min-overlap=1.1"],
        vec!["--intergenic-trim-bases=-1"],
        vec!["--intronic-trim-bases=-1"],
        vec!["--junction-trim-bases=-1"],
        vec!["--strandedness=invalid"],
        vec!["--endedness=invalid"],
    ] {
        assert!(!run(&input, &output, &args).status.success(), "{args:?}");
        assert!(!output.exists());
    }
    assert_no_temporary_files(dir.path());
}

#[test]
fn missing_eof_and_corrupt_input_never_publish_a_result() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let original = fs::read(&input).unwrap();
    let output = dir.path().join("output.bam");
    fs::write(&input, &original[..original.len() - 28]).unwrap();
    assert_failure(&run(&input, &output, &[]), "EOF integrity check failed");
    assert!(!output.exists());
    let mut corrupt = original;
    // Damage compressed content while retaining the final EOF block.
    let offset = corrupt.len() / 2;
    corrupt[offset] ^= 0xff;
    fs::write(&input, corrupt).unwrap();
    assert!(!run(&input, &output, &[]).status.success());
    assert!(!output.exists());
    assert_no_temporary_files(dir.path());
}

#[test]
fn out_of_bounds_alignment_fails_without_output() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let (header, mut records) = read_bam(&input);
    records[9].set_pos(999);
    write_bam(&input, &header, &records);
    let output = dir.path().join("output.bam");
    assert_failure(
        &run(&input, &output, &[]),
        "Alignment outside reference bounds",
    );
    assert!(!output.exists());
    assert_no_temporary_files(dir.path());
}

#[test]
fn missing_input_or_output_directory_fails_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    let output = dir.path().join("output.bam");
    assert!(!run(&input, &output, &[]).status.success());
    assert!(!output.exists());
    input_bam(&input);
    assert!(!run(&input, &dir.path().join("missing/output.bam"), &[])
        .status
        .success());
    assert_no_temporary_files(dir.path());
}
