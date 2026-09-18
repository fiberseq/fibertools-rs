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
// reverse read.
#[test]
fn extract_drops_annotations_that_exceed_the_sequence() {
    let out = run(&[
        "extract",
        "--all",
        "-",
        fixture("ont_hardclip_supplementary.bam").to_str().unwrap(),
    ]);
    let mut lines = out.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (flag, nuc, msp, m6a) = (
        col("sam_flag"),
        col("nuc_starts"),
        col("msp_starts"),
        col("m6a"),
    );
    let mut n_primary = 0;
    let mut n_supp = 0;
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        let supplementary = f[flag].parse::<u16>().unwrap() & 2048 != 0;
        for c in [nuc, msp, m6a] {
            assert_eq!(f[c] == ".", supplementary, "{}: {}", header[c], f[c]);
        }
        if supplementary {
            n_supp += 1;
        } else {
            n_primary += 1;
        }
    }
    assert_eq!((n_primary, n_supp), (2, 2));
}
