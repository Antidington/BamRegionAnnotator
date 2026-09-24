mod common;

use common::*;
use rust_htslib::bam::{HeaderView, Record};

#[test]
fn original_annotation_results_are_preserved() {
    let cases: &[(&str, &[&str])] = &[
        ("default", &[]),
        (
            "all_tags",
            &[
                "--include-introns",
                "--add-gene-tags",
                "--add-tx-tag",
                "--add-an-tag",
            ],
        ),
        (
            "reverse_five_prime",
            &[
                "--strandedness",
                "reverse",
                "--endedness",
                "five_prime",
                "--include-introns",
                "--add-gene-tags",
                "--add-tx-tag",
                "--add-an-tag",
            ],
        ),
        (
            "trimmed",
            &[
                "--intergenic-trim-bases",
                "5",
                "--intronic-trim-bases",
                "5",
                "--junction-trim-bases",
                "5",
                "--region-min-overlap",
                "0.75",
                "--include-introns",
                "--add-gene-tags",
                "--add-tx-tag",
                "--add-an-tag",
            ],
        ),
    ];
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.bam");
    input_bam(&input);
    let (input_header, _) = read_bam(&input);
    for (name, args) in cases {
        let output = dir.path().join(format!("{name}.bam"));
        assert_success(&run(&input, &output, args));
        let (header, records) = read_bam(&output);
        assert_eq!(header.to_bytes(), input_header.to_bytes());
        let view = HeaderView::from_header(&header);
        let sam = std::fs::read_to_string(fixture(&format!("expected/{name}.sam"))).unwrap();
        let expected: Vec<_> = sam
            .lines()
            .map(|line| Record::from_sam(&view, line.as_bytes()).unwrap())
            .collect();
        assert_eq!(records.len(), expected.len(), "{name}");
        for (actual, expected) in records.iter().zip(&expected) {
            assert_eq!(
                actual,
                expected,
                "{name}: {}",
                String::from_utf8_lossy(expected.qname())
            );
        }
    }
}
