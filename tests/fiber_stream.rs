use fibertools_rs::fiber::FiberseqRecords;
use fibertools_rs::utils::input_bam::FiberFilters;
use rust_htslib::bam::{self, Read};

/// A chunk whose records are all filtered out must not terminate the
/// stream. Regression test: BamChunk::next used to return None (meaning
/// "exhausted") whenever one pull of chunk_size records filtered to
/// empty, silently dropping everything after e.g. a long run of
/// secondary alignments under -F 256.
#[test]
fn stream_survives_a_fully_filtered_chunk() {
    // chunk_size caps at 2,500; a longer run guarantees at least one
    // all-filtered pull regardless of thread count.
    const N_SECONDARY: usize = 3_000;
    const N_PRIMARY: usize = 5;

    let tmp = tempfile::NamedTempFile::with_suffix(".bam").unwrap();
    let header = bam::Header::new();
    {
        let mut writer = bam::Writer::from_path(tmp.path(), &header, bam::Format::Bam).unwrap();
        let mut rec = bam::Record::new();
        rec.set(b"read", None, b"ACGTACGT", &[255u8; 8]);
        rec.set_tid(-1);
        rec.set_pos(-1);
        for i in 0..(N_SECONDARY + N_PRIMARY) {
            let flags = if i < N_SECONDARY { 0x100 } else { 0x4 };
            rec.set_flags(flags);
            writer.write(&rec).unwrap();
        }
    }

    let mut bam = bam::Reader::from_path(tmp.path()).unwrap();
    let filters = FiberFilters {
        bit_flag: Some(0x100),
        ..FiberFilters::default()
    };
    let n = FiberseqRecords::new(&mut bam, filters).count();
    assert_eq!(
        n, N_PRIMARY,
        "reads after a fully-filtered chunk must still stream"
    );
}
