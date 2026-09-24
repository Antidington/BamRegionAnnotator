use tx_annotation::transcript::AnnotationRegion;

#[derive(Default)]
pub(crate) struct Stats {
    pub(crate) total: u64,
    pub(crate) mapped: u64,
    unmapped: u64,
    exonic: u64,
    intronic: u64,
    intergenic: u64,
}

impl Stats {
    pub(crate) fn observe_unmapped(&mut self) {
        self.total += 1;
        self.unmapped += 1;
    }

    pub(crate) fn observe_region(&mut self, region: AnnotationRegion) {
        self.total += 1;
        self.mapped += 1;
        match region {
            AnnotationRegion::Exonic => self.exonic += 1,
            AnnotationRegion::Intronic => self.intronic += 1,
            AnnotationRegion::Intergenic => self.intergenic += 1,
        }
    }

    pub(crate) fn report(&self) {
        println!("Annotation complete!");
        println!("Total alignment records: {}", self.total);
        println!("Mapped alignment records: {}", self.mapped);
        println!("Unmapped alignment records: {}", self.unmapped);
        println!("Region percentages use mapped alignment records as denominator.");
        for (label, count) in [
            ("Exonic", self.exonic),
            ("Intronic", self.intronic),
            ("Intergenic", self.intergenic),
        ] {
            println!(
                "{label} records: {count} ({:.2}%)",
                percentage(count, self.mapped)
            );
        }
    }
}

fn percentage(part: u64, total: u64) -> f64 {
    if total == 0 {
        0.0
    } else {
        100.0 * part as f64 / total as f64
    }
}
