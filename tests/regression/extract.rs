use super::common::{fixture, run, select_bed12_cols, select_tsv_cols};
use tempfile::NamedTempFile;

// bed12 columns worth snapshotting: locator + per-record feature data.
// score/thick_*/item_rgb are constants or duplicates for these outputs.
const BED12_COLS: &[&str] = &[
    "chrom",
    "start",
    "end",
    "name",
    "strand",
    "block_count",
    "block_sizes",
    "block_starts",
];

#[test]
fn extract_m6a() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--m6a",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

#[test]
fn extract_nuc() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--nuc",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

#[test]
fn extract_msp() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("all.bam").to_str().unwrap(),
        "--msp",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

// NAPA.bam has FIRE annotations (aq tag) — exercises FIRE in the --all tabular output
// Select MA-relevant columns so new output columns don't count as regressions.
#[test]
fn extract_all_napa() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("NAPA.bam").to_str().unwrap(),
        "--all",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#ct",
            "st",
            "en",
            "fiber",
            "nuc_starts",
            "nuc_lengths",
            "msp_starts",
            "msp_lengths",
            "fire_starts",
            "fire_lengths",
            "fire_qual",
            "m6a",
        ]
    ));
}

// `ma_spelled.bam` is a committed fixture carrying only the canonical
// Ma/Aq/An tag spellings (samtools/hts-specs#862), generated once with
// `ft convert-tags` from msp_nuc.bam. It pins the on-disk format
// independently of the current writer: a same-build write+read round trip
// can be self-consistently wrong, this fixture cannot.
#[test]
fn extract_reads_ma_spelled_fixture() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "extract",
        fixture("ma_spelled.bam").to_str().unwrap(),
        "--nuc",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_bed12_cols(&out, BED12_COLS));
}

// The reader drops annotations that do not fit SEQ (hard-clipped supplementary
// reads keep the full-length read's tags, #136). Before this, extract printed
// misplaced coordinates for the forward read and u32-wrapped ones for the
// reverse read. Supplementaries must report `.` for nuc, msp and m6a; the
// primaries must not.
fn assert_supplementaries_untagged(bam: &str, n_primary: usize, n_supp: usize) {
    let out = run(&["extract", "--all", "-", fixture(bam).to_str().unwrap()]);
    let mut lines = out.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (flag, nuc, msp, m6a) = (
        col("sam_flag"),
        col("nuc_starts"),
        col("msp_starts"),
        col("m6a"),
    );
    let (mut seen_primary, mut seen_supp) = (0, 0);
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let supplementary = f[flag].parse::<u16>().unwrap() & 2048 != 0;
        for c in [nuc, msp, m6a] {
            assert_eq!(f[c] == ".", supplementary, "{bam} {}: {}", header[c], f[c]);
        }
        if supplementary {
            seen_supp += 1;
        } else {
            seen_primary += 1;
        }
    }
    assert_eq!((seen_primary, seen_supp), (n_primary, n_supp), "{bam}");
}

// Legacy ns/nl/as/al tags, no MM/ML (dorado aligner strips them): caught by
// the hard-clip rule for legacy tags.
#[test]
fn extract_drops_legacy_annotations_on_hard_clipped_reads() {
    assert_supplementaries_untagged("ont_hardclip_supplementary.bam", 2, 2);
}

// MM/ML/MN copied verbatim onto a 2376H hard-clipped supplementary (the
// plain minimap2 shape): caught by MN != SEQ length. Until 0.14 this record
// was silently removed from every output instead.
#[test]
fn extract_drops_mm_ml_on_hard_clipped_reads() {
    assert_supplementaries_untagged("ont_hardclip_mmml.bam", 1, 1);
}
