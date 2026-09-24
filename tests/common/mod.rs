#![allow(dead_code)]

use rust_htslib::bam::{self, Read};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

pub fn fixture(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(path)
}

pub fn input_bam(path: &Path) {
    let mut reader = bam::Reader::from_path(fixture("reads.sam")).unwrap();
    let header = bam::Header::from_template(reader.header());
    let records: Vec<_> = reader.records().map(Result::unwrap).collect();
    write_bam(path, &header, &records);
}

pub fn write_bam(path: &Path, header: &bam::Header, records: &[bam::Record]) {
    let mut writer = bam::Writer::from_path(path, header, bam::Format::Bam).unwrap();
    for record in records {
        writer.write(record).unwrap();
    }
}

pub fn read_bam(path: &Path) -> (bam::Header, Vec<bam::Record>) {
    let mut reader = bam::Reader::from_path(path).unwrap();
    let header = bam::Header::from_template(reader.header());
    (header, reader.records().map(Result::unwrap).collect())
}

pub fn run(input: &Path, output: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_bam_region_annotator"))
        .arg("--reference")
        .arg(fixture("reference"))
        .arg("--input")
        .arg(input)
        .arg("--output")
        .arg(output)
        .args(args)
        .output()
        .unwrap()
}

pub fn assert_success(result: &Output) {
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}
