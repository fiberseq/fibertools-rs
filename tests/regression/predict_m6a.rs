use super::common::{fixture, ft};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{Read, Reader};
use std::collections::BTreeSet;
use std::process::Command;
use tempfile::NamedTempFile;

/// Regression guard for the MM/ML strand convention emitted by the m6a
/// producer path (`predict-m6a`). The command writes a BAM, so we read it back
/// and project to the *distinct* MM group headers — skip-base + strand + mod
/// code, e.g. `A+a`, `T-a`, `C+m`.
///
/// This projection is intentionally coarse: it depends only on which base is
/// modified, not on the model's exact per-base scores, so it is deterministic
/// across ML backends (candle vs libtorch differ on counts/positions — see
/// `tests/m6a_prediction_test.rs` — but never on the base→strand mapping).
/// It is exactly the axis a strand-encoding bug corrupts: minus-strand m6a
/// must serialize as `T-a`, never `T+a`. Snapshotting ML values or call counts
/// would be flaky; the header set is not.
#[test]
fn mm_strand_groups() {
    let out_bam = NamedTempFile::new().unwrap();
    let status = Command::new(ft())
        .args([
            "predict-m6a",
            "-t",
            "1",
            fixture("revio.bam").to_str().unwrap(),
            out_bam.path().to_str().unwrap(),
        ])
        .status()
        .expect("failed to spawn ft predict-m6a");
    assert!(status.success(), "ft predict-m6a exited {status}");

    let mut reader = Reader::from_path(out_bam.path()).expect("open predicted bam");
    let mut groups: BTreeSet<String> = BTreeSet::new();
    for rec in reader.records() {
        let record = rec.expect("read record");
        if let Ok(Aux::String(mm)) = record.aux(b"MM") {
            for group in mm.split(';').filter(|g| !g.is_empty()) {
                // The group header is everything before the first delta count,
                // e.g. `A+a` in `A+a,5,3,...`.
                let header = group.split(',').next().unwrap_or("");
                groups.insert(header.to_string());
            }
        }
    }

    let projection = groups.into_iter().collect::<Vec<_>>().join("\n");
    insta::assert_snapshot!(projection);
}
