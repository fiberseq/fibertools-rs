use super::common::{fixture, run, run_capture, select_bed12_cols, select_tsv_cols};
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

/// A hard-clipped supplementary whose tags describe the full-length read
/// (#136): nuc/msp are kept in the full read's frame and lifted through the
/// hard clip, m6A is dropped. Expected values were computed with pysam
/// get_aligned_pairs, independently of ft.
struct FullFrame {
    /// qname prefix
    qname: &'static str,
    flag: u16,
    /// SEQ length
    /// The annotation frame: on a full-read-frame record the whole read,
    /// not SEQ, so it matches the molecular nuc/msp columns and the
    /// molecular-mode BED12 end.
    fiber_length: i64,
    /// every nucleosome of the full read stays in the tag
    n_nuc: usize,
    /// nuc_starts[0]: BAM orientation, full-read frame
    first_nuc_start: i64,
    /// ref_nuc_starts entries != -1, in output order (a prefix when shorter
    /// than n_lifted)
    lifted: &'static [i64],
    n_lifted: usize,
    same_nucs_as_primary: bool,
}

fn assert_full_frame_supplementaries(bam: &str, n_primary: usize, expected: &[FullFrame]) {
    let out = run(&["extract", "--all", "-", fixture(bam).to_str().unwrap()]);
    let mut lines = out.lines();
    let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let (fiber, flag, len, nuc, ref_nuc, m6a) = (
        col("fiber"),
        col("sam_flag"),
        col("fiber_length"),
        col("nuc_starts"),
        col("ref_nuc_starts"),
        col("m6a"),
    );
    let ints = |s: &str| -> Vec<i64> {
        s.trim_end_matches(',')
            .split(',')
            .map(|x| x.parse().unwrap())
            .collect()
    };
    let mut primary_nucs = std::collections::HashMap::new();
    let mut supp_rows = Vec::new();
    for line in lines {
        let f: Vec<String> = line.split('\t').map(str::to_string).collect();
        if f[flag].parse::<u16>().unwrap() & 2048 == 0 {
            assert_ne!(f[m6a], ".", "{bam} {}: primary lost m6A", f[fiber]);
            primary_nucs.insert(f[fiber].clone(), f[nuc].clone());
        } else {
            supp_rows.push(f);
        }
    }
    assert_eq!(
        (primary_nucs.len(), supp_rows.len()),
        (n_primary, expected.len()),
        "{bam}"
    );
    for f in supp_rows {
        let sam_flag: u16 = f[flag].parse().unwrap();
        let e = expected
            .iter()
            .find(|e| f[fiber].starts_with(e.qname) && e.flag == sam_flag)
            .unwrap_or_else(|| {
                panic!(
                    "{bam}: unexpected supplementary {} flag {sam_flag}",
                    f[fiber]
                )
            });
        assert_eq!(
            f[len].parse::<i64>().unwrap(),
            e.fiber_length,
            "{bam} {}",
            e.qname
        );
        assert_eq!(
            f[m6a], ".",
            "{bam} {} {}: m6A must be dropped on a full-read frame",
            e.qname, e.flag
        );
        let nucs = ints(&f[nuc]);
        assert_eq!(
            nucs.len(),
            e.n_nuc,
            "{bam} {} {}: nuc_starts",
            e.qname,
            e.flag
        );
        assert_eq!(
            nucs[0], e.first_nuc_start,
            "{bam} {} {}: nuc_starts[0]",
            e.qname, e.flag
        );
        assert!(
            nucs.windows(2).all(|w| w[0] < w[1]),
            "{bam}: nuc_starts not ascending"
        );
        if e.same_nucs_as_primary {
            assert_eq!(
                f[nuc], primary_nucs[&f[fiber]],
                "{bam} {}: molecular nuc_starts must be unchanged",
                e.qname
            );
        }
        let refs = ints(&f[ref_nuc]);
        assert_eq!(refs.len(), e.n_nuc);
        let lifted: Vec<i64> = refs.into_iter().filter(|r| *r != -1).collect();
        assert_eq!(
            lifted.len(),
            e.n_lifted,
            "{bam} {} {}: lifted nucleosomes",
            e.qname,
            e.flag
        );
        assert_eq!(
            &lifted[..e.lifted.len()],
            e.lifted,
            "{bam} {} {}: ref_nuc_starts",
            e.qname,
            e.flag
        );
    }
}

#[test]
fn extract_lifts_full_frame_ma_records() {
    assert_full_frame_supplementaries(
        "ont_hardclip_full_frame.bam",
        2,
        &[
            FullFrame {
                qname: "f2009f4d",
                flag: 2048,
                fiber_length: 5376,
                n_nuc: 31,
                first_nuc_start: 64,
                n_lifted: 13,
                lifted: &[
                    2171, 2392, 2557, 2698, 2862, 3006, 3209, 3438, 3626, 3837, 3978, 4218, 4380,
                ],
                same_nucs_as_primary: true,
            },
            FullFrame {
                qname: "f2009f4d",
                flag: 2064,
                fiber_length: 5376,
                n_nuc: 31,
                first_nuc_start: 73,
                n_lifted: 13,
                lifted: &[
                    2207, 2333, 2575, 2721, 2925, 3093, 3231, 3498, 3666, 3806, 3982, 4166, 4359,
                ],
                same_nucs_as_primary: false,
            },
            FullFrame {
                qname: "7b40cfd0",
                flag: 2048,
                fiber_length: 9693,
                n_nuc: 44,
                first_nuc_start: 84,
                n_lifted: 6,
                lifted: &[102454, 102616, 102810, 102955, 103409, 103554],
                same_nucs_as_primary: true,
            },
        ],
    );
}

// Legacy ns/nl/as/al on hard clips (dorado aligner shape): full-read frame
// with read_length = SEQ + H_lead + H_trail (29940 = 7440+22492+8;
// 33088 = 9910+0+23178). The 8ac3be13 nucleosome that straddles the leading
// clip snaps to the first aligned base (3834036), like across a soft clip.
#[test]
fn extract_lifts_full_frame_legacy_records() {
    assert_full_frame_supplementaries(
        "ont_hardclip_supplementary.bam",
        2,
        &[
            FullFrame {
                qname: "8ac3be13",
                flag: 2048,
                fiber_length: 29940,
                n_nuc: 150,
                first_nuc_start: 172,
                n_lifted: 34,
                lifted: &[
                    3834036, 3834112, 3834330, 3834492, 3834646, 3834860, 3834979, 3835452,
                    3835701, 3835898, 3836264, 3836417, 3836689, 3836805,
                ],
                same_nucs_as_primary: false,
            },
            FullFrame {
                qname: "4bd15181",
                flag: 2064,
                fiber_length: 33088,
                n_nuc: 162,
                first_nuc_start: 217,
                n_lifted: 48,
                lifted: &[
                    10475654, 10475848, 10476012, 10476197, 10476540, 10476719, 10476912, 10477059,
                    10477293, 10477492, 10477716, 10477893, 10478079,
                ],
                same_nucs_as_primary: false,
            },
        ],
    );
}

// One WARN per run that m6A was dropped, with the realignment remedy; no
// stale-frame noise.
#[test]
fn full_frame_m6a_dropped_warns_once() {
    let (_, err) = run_capture(&[
        "extract",
        "--all",
        "-",
        fixture("ont_hardclip_full_frame.bam").to_str().unwrap(),
    ]);
    assert_eq!(
        err.matches("their m6A (MM/ML) describes bases this record does not carry and was dropped")
            .count(),
        1,
        "{err}"
    );
    assert!(err.contains("minimap2 -Y -y"), "remedy missing: {err}");
    assert!(
        !err.contains("dropping annotations for"),
        "full-frame records are not stale: {err}"
    );
}

// MM/ML/MN copied verbatim onto a 2376H hard-clipped supplementary (the
// plain minimap2 shape): caught by MN != SEQ length. Until 0.14 this record
// was deleted from every output with a per-record warning instead.
#[test]
fn extract_drops_mm_ml_on_hard_clipped_reads() {
    assert_supplementaries_untagged("ont_hardclip_mmml.bam", 1, 1);
}
