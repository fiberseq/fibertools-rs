use super::common::{fixture, run, select_tsv_cols};

// Subset to 3 motifs and cap features to ±200 bp of center to keep the
// snapshot small. The center math is uniform across distance, so this
// still exercises liftover + strand flipping for m6a/nuc/msp types.
#[test]
fn center_default() {
    let out = run(&[
        "center",
        fixture("center.bam").to_str().unwrap(),
        "--bed",
        fixture("center.small.bed").to_str().unwrap(),
        "--dist",
        "200",
    ]);
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "chrom",
            "centering_position",
            "strand",
            "query_name",
            "centered_position_type",
            "centered_start",
            "centered_end",
        ]
    ));
}

// Centering on ref 2171 (the first retained nucleosome of the forward
// full-frame record). Rows of the two 2400 bp full-frame records
// (query_length 2400; the primary is 5376) must use the clip offset:
// molecular mode reproduces the primary's rows for the forward record and
// the reverse record's flipped ones; reference mode lists only the
// nucleosomes that lift. Expected values from pysam (#136).
#[test]
fn center_full_frame_records_use_the_clip_offset() {
    let bam = fixture("ont_hardclip_full_frame.bam");
    let bed = fixture("ont_hardclip_full_frame.center.bed");
    let rows = |extra: &[&str], query_length: &str| -> Vec<(i64, i64)> {
        let mut args = vec![
            "center",
            bam.to_str().unwrap(),
            "--bed",
            bed.to_str().unwrap(),
            "--dist",
            "200",
        ];
        args.extend_from_slice(extra);
        let out = run(&args);
        let mut lines = out.lines();
        let header: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let col = |n: &str| header.iter().position(|h| *h == n).unwrap();
        let (strand, qlen, ty, st, en, cqs, cqe) = (
            col("strand"),
            col("query_length"),
            col("centered_position_type"),
            col("centered_start"),
            col("centered_end"),
            col("centered_query_start"),
            col("centered_query_end"),
        );
        assert!(
            lines
                .clone()
                .any(|l| l.split('\t').nth(strand) == Some("-")),
            "minus-strand centering must produce rows"
        );
        let mut v: Vec<(i64, i64)> = lines
            .map(|l| l.split('\t').collect::<Vec<_>>())
            .filter(|f| f[strand] == "+" && f[ty] == "nuc" && f[qlen] == query_length)
            .inspect(|f| {
                // molecular mode: leading columns are SEQ-relative (anchor at
                // SEQ position 37, SEQ length 2400)
                if query_length == "2400" && extra.is_empty() {
                    assert_eq!(
                        (f[cqs], f[cqe]),
                        ("-37", "2363"),
                        "SEQ frame leading columns"
                    );
                }
            })
            .map(|f| (f[st].parse().unwrap(), f[en].parse().unwrap()))
            .collect();
        v.sort();
        v
    };
    // the primary, for reference (unchanged behaviour)
    assert_eq!(rows(&[], "5376"), vec![(-163, -62), (0, 138)]);
    assert_eq!(rows(&["--reference"], "5376"), vec![(-164, -62), (0, 137)]);
    // forward full-frame record: same rows as the primary; reverse record: (36,116),(163,280)
    assert_eq!(
        rows(&[], "2400"),
        vec![(-163, -62), (0, 138), (36, 116), (163, 280)]
    );
    // reference mode: forward keeps only its lifted nucleosome, reverse its two
    assert_eq!(
        rows(&["--reference"], "2400"),
        vec![(0, 137), (36, 116), (162, 278)]
    );
}
