use super::bamlift::*;
use super::extract::*;
use rust_htslib::bam::record;
use rust_htslib::{bam, bam::ext::BamRecordExtensions, bam::HeaderView};

pub struct CenteredFiberData {
    fiber: FiberseqData,
    record: record::Record,
    offset: i64,
    forward: bool,
    reference_position: i64,
}

impl CenteredFiberData {
    pub fn new(
        fiber: FiberseqData,
        record: bam::Record,
        reference_position: i64,
        forward: bool,
    ) -> Self {
        let offset = CenteredFiberData::find_offset(&record, reference_position);
        Self {
            fiber,
            record,
            offset,
            forward,
            reference_position,
        }
    }

    pub fn get_centering_position(&self) -> i64 {
        self.reference_position
    }

    // TODO
    fn find_offset(record: &bam::Record, reference_position: i64) -> i64 {
        let z = get_exact_query_positions(record, &[reference_position]);
        if z.is_empty() {
            0
        } else {
            z[0] - reference_position
        }
    }

    // TODO
    fn apply_offset(&self, positions: &[i64]) -> Vec<i64> {
        if self.forward {
            positions.iter().map(|&p| p + self.offset).collect()
        } else {
            positions.iter().map(|&p| p + self.offset).collect()
        }
    }

    pub fn m6a_positions(&self) -> Vec<i64> {
        self.apply_offset(&self.fiber.base_mods.m6a_positions(false))
    }

    // TODO
    pub fn write(&self, head_view: &HeaderView) -> String {
        let ct = std::str::from_utf8(head_view.tid2name(self.record.tid() as u32)).unwrap();
        format!(
            "{ct}\t{}\t{}",
            self.record.reference_start(),
            self.record.reference_end()
        )
    }
}

pub fn center(
    records: Vec<bam::Record>,
    reference_position: i64,
    forward: bool,
    head_view: &HeaderView,
) {
    let fiber_data = FiberseqData::from_records(&records);
    let iter = fiber_data.into_iter().zip(records.into_iter());
    for (fiber, record) in iter {
        let centered_fiber = CenteredFiberData::new(fiber, record, reference_position, forward);
        centered_fiber.write(head_view);
    }
}
