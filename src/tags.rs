use anyhow::{ensure, Context, Result};
use rust_htslib::bam::{record::Aux, Record};
use tx_annotation::transcript::{AnnotationRegion, TranscriptAnnotator};

#[derive(Clone, Copy)]
pub(crate) struct TagOptions {
    pub(crate) genes: bool,
    pub(crate) transcripts: bool,
    pub(crate) antisense: bool,
}

pub(crate) fn annotate(
    record: &mut Record,
    annotator: &TranscriptAnnotator,
    options: TagOptions,
) -> Result<AnnotationRegion> {
    check_target_tags(record, options)?;
    let annotation = annotator.annotate_alignment(record);
    if let Some(tag) = annotation.make_re_tag() {
        record.push_aux(b"RE", Aux::Char(tag.as_bytes()[0]))?;
    }
    // Optionally add gene tags (GX, GN)
    if options.genes {
        if let Some((gx_tag, gn_tag)) = annotation.make_gx_gn_tags() {
            record.push_aux(b"GX", Aux::String(&gx_tag))?;
            record.push_aux(b"GN", Aux::String(&gn_tag))?;
        }
    }

    // Optionally add transcript tag (TX)
    if options.transcripts {
        if let Some(tx_tag) = annotation.make_tx_tag() {
            record.push_aux(b"TX", Aux::String(&tx_tag))?;
        }
    }

    // Optionally add antisense transcript tag (AN)
    if options.antisense {
        if let Some(an_tag) = annotation.make_an_tag() {
            record.push_aux(b"AN", Aux::String(&an_tag))?;
        }
    }

    Ok(annotation.region)
}

fn check_target_tags(record: &Record, options: TagOptions) -> Result<()> {
    for entry in record.aux_iter() {
        let (tag, _) = entry.context("Invalid auxiliary tag")?;
        let target = tag == b"RE"
            || (options.genes && (tag == b"GX" || tag == b"GN"))
            || (options.transcripts && tag == b"TX")
            || (options.antisense && tag == b"AN");
        ensure!(
            !target,
            "Record {} already has target tag {}",
            String::from_utf8_lossy(record.qname()),
            String::from_utf8_lossy(tag)
        );
    }
    Ok(())
}
