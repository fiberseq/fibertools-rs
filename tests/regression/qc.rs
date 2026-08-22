use super::common::{fixture, run, select_tsv_cols, tagged_bam};
use tempfile::NamedTempFile;

#[test]
fn qc_default() {
    let out = run(&["qc", fixture("all.bam").to_str().unwrap()]);
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &["statistic", "value", "count", "count_filtered"]
    ));
}

// FIRE quals live on the `fire` annotation type (MA spec), not the MSPs;
// the m6a_per_msp_size statistic must overlay them or is_fire is always
// false for fire-scored BAMs.
#[test]
fn qc_m6a_per_msp_sees_fire_elements() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("all.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let out = run(&["qc", "-m", scored.path().to_str().unwrap()]);
    let fire_rows: Vec<&str> = out
        .lines()
        .filter(|l| {
            l.starts_with("m6a_per_msp_size") && l.split('\t').nth(1).unwrap().ends_with(",true")
        })
        .collect();
    assert!(
        !fire_rows.is_empty(),
        "no m6a_per_msp_size rows with is_fire=true; fire quals were not overlaid onto MSPs"
    );
    insta::assert_snapshot!(fire_rows.join("\n"));
}

/// Parse "stat\tvalue\tcount\tcount_filtered" rows for one statistic.
fn rows_for<'a>(out: &'a str, stat: &str) -> Vec<(f64, i64, i64)> {
    out.lines()
        .filter(|l| l.split('\t').next() == Some(stat))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[1].parse().unwrap_or(f64::NAN),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
            )
        })
        .collect()
}

/// On a BAM with the tag on disk, the filtered column must obey the span
/// invariant: the surviving nuc+MSP calls tile the span exactly, so
/// sum(nuc_length x filtered) + sum(msp_length x filtered) equals the
/// filtered phased_bp, which equals sum(fiber_length keys x filtered).
/// This is the strongest available proof the bp denominator is right.
#[test]
fn qc_filtered_columns_on_tagged_bam() {
    // add-nucleosomes writes the tag on disk (not just read-time backfill).
    // Note it drops FIRE_TYPE, so no is_fire assertions on this BAM.
    let tagged = tagged_bam("all.bam");
    let out = run(&["qc", tagged.path().to_str().unwrap()]);

    let fl = rows_for(&out, "fiber_length");
    assert_eq!(fl.len(), 44, "22 read-length + 22 span-length union keys");
    let max_count_key = fl
        .iter()
        .filter(|r| r.1 > 0)
        .map(|r| r.0)
        .fold(f64::MIN, f64::max);
    let max_filt_key = fl
        .iter()
        .filter(|r| r.2 > 0)
        .map(|r| r.0)
        .fold(f64::MIN, f64::max);
    assert!(
        max_filt_key <= max_count_key,
        "span keys cannot exceed read keys"
    );

    let weighted = |stat: &str, col: fn(&(f64, i64, i64)) -> i64| -> f64 {
        rows_for(&out, stat)
            .iter()
            .map(|r| r.0 * col(r) as f64)
            .sum()
    };
    let filt_bp = rows_for(&out, "phased_bp")[0].2 as f64;
    let all_bp = rows_for(&out, "phased_bp")[0].1 as f64;
    assert_eq!(
        weighted("nuc_length", |r| r.2) + weighted("msp_length", |r| r.2),
        filt_bp,
        "surviving calls must tile the filtered bp exactly"
    );
    assert_eq!(weighted("fiber_length", |r| r.2), filt_bp);
    assert!(filt_bp < all_bp, "span bp strictly below read bp");
}

/// Mixed callability: three_two.bam has 20 NotCallable reads under the
/// default minimums, zero-nuc reads (inf key), and degenerate-span reads.
#[test]
fn qc_mixed_callability() {
    let tagged = tagged_bam("three_two.bam");
    let out = run(&["qc", tagged.path().to_str().unwrap()]);

    // NotCallable reads are excluded from every filtered tally: the read
    // count in the filtered column is below the unfiltered one.
    let pr = rows_for(&out, "phased_reads");
    let (all, filt): (i64, i64) = pr.iter().fold((0, 0), |a, r| (a.0 + r.1, a.1 + r.2));
    assert_eq!(all, 92);
    assert_eq!(filt, 72, "72 of 92 reads clear the callability minimums");

    // Zero-nuc reads keep their inf key in the unfiltered column.
    let inf_rows: Vec<&str> = out
        .lines()
        .filter(|l| l.starts_with("read_length_per_nuc\tinf"))
        .collect();
    assert!(
        !inf_rows.is_empty(),
        "inf key preserved in unfiltered column"
    );
    for r in &inf_rows {
        assert_eq!(r.split('\t').nth(3), Some("0"), "inf never enters filtered");
    }
    // And no NaN key in either column.
    assert!(!out.contains("NaN"), "no NaN keys anywhere");
}

/// --filtered-min-callable-length layers on top of the callable state:
/// it shrinks count_filtered while leaving every count byte-identical.
#[test]
fn qc_filtered_min_callable_length() {
    let tagged = tagged_bam("all.bam");
    let base = run(&["qc", tagged.path().to_str().unwrap()]);
    // all.bam spans run ~9k-49k; 15000 splits them.
    let cut = run(&[
        "qc",
        "--filtered-min-callable-length",
        "15000",
        tagged.path().to_str().unwrap(),
    ]);
    // Keyed, not positional: a read failing the threshold no longer inserts
    // its span-key row at all, so the row sets differ while every surviving
    // (statistic, value) keeps its count.
    let counts = |out: &str| -> std::collections::HashMap<(String, String), String> {
        out.lines()
            .skip(1)
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                ((f[0].to_string(), f[1].to_string()), f[2].to_string())
            })
            .collect()
    };
    let (b_counts, c_counts) = (counts(&base), counts(&cut));
    for k in b_counts.keys().chain(c_counts.keys()) {
        let get = |m: &std::collections::HashMap<(String, String), String>| {
            m.get(k).map(String::as_str).unwrap_or("0").to_string()
        };
        assert_eq!(get(&b_counts), get(&c_counts), "count changed for {k:?}");
    }
    let filt_total = |out: &str| -> i64 { rows_for(out, "phased_reads").iter().map(|r| r.2).sum() };
    let (b, c) = (filt_total(&base), filt_total(&cut));
    assert!(c < b, "threshold must shrink filtered reads ({c} !< {b})");
    assert!(c > 0, "threshold must not empty the filtered set");
}

/// The callability minimums select the filtered column in ft qc and never
/// drop reads: the unfiltered totals stay complete.
#[test]
fn qc_minimums_never_shrink_unfiltered() {
    let out = run(&[
        "qc",
        "--min-msp",
        "50",
        fixture("all.bam").to_str().unwrap(),
    ]);
    let pr = rows_for(&out, "phased_reads");
    let (all, filt): (i64, i64) = pr.iter().fold((0, 0), |a, r| (a.0 + r.1, a.1 + r.2));
    assert_eq!(all, 22, "unfiltered total must stay complete");
    assert!(filt < 22, "--min-msp 50 selects the filtered column");
}

/// Every line matches the header's field count, m6a_acf rows included
/// (select_tsv_cols silently tolerates ragged trailing fields).
#[test]
fn qc_output_is_never_ragged() {
    let tagged = tagged_bam("all.bam");
    let out = run(&["qc", "--acf", tagged.path().to_str().unwrap()]);
    let n = out.lines().next().unwrap().split('\t').count();
    for l in out.lines() {
        assert_eq!(l.split('\t').count(), n, "ragged line: {l}");
    }
}

/// Dual ACF: subsample first, filter within that sample. Deterministic here
/// because the huge sample rate keeps exactly the first qualifying reads.
#[test]
fn qc_dual_acf() {
    let tagged = tagged_bam("three_two.bam");
    let args = [
        "qc",
        "--acf",
        "--acf-max-reads",
        "5",
        "--acf-sample-rate",
        "1000000000",
        tagged.path().to_str().unwrap(),
    ];
    let out = run(&args);
    // Reproducible without a seed at p(sample) ~ 1e-9.
    assert_eq!(out, run(&args));

    let acf_reads: Vec<&str> = out.lines().filter(|l| l.starts_with("acf_reads")).collect();
    assert_eq!(acf_reads.len(), 1);
    let f: Vec<&str> = acf_reads[0].split('\t').collect();
    assert_eq!(f[2], "5", "reservoir holds the first 5 qualifying reads");
    let passing: i64 = f[3].parse().unwrap();
    assert!(passing <= 5);

    // Every m6a_acf row has a real filtered value (4 numeric fields).
    let acf_rows: Vec<&str> = out.lines().filter(|l| l.starts_with("m6a_acf")).collect();
    assert!(!acf_rows.is_empty());
    for r in &acf_rows {
        let f: Vec<&str> = r.split('\t').collect();
        assert_eq!(f.len(), 4);
        assert!(
            f[3].parse::<f64>().is_ok(),
            "filtered ACF must be numeric: {r}"
        );
    }
}

/// A filter no sampled read satisfies: filtered ACF is NA, exit 0.
#[test]
fn qc_dual_acf_empty_filtered_set() {
    let tagged = tagged_bam("three_two.bam");
    let out = run(&[
        "qc",
        "--acf",
        "--acf-max-reads",
        "5",
        "--acf-sample-rate",
        "1000000000",
        "--filtered-min-callable-length",
        "100000000",
        tagged.path().to_str().unwrap(),
    ]);
    let acf_rows: Vec<&str> = out.lines().filter(|l| l.starts_with("m6a_acf")).collect();
    assert!(!acf_rows.is_empty());
    for r in &acf_rows {
        assert_eq!(
            r.split('\t').nth(3),
            Some("NA"),
            "empty filtered set is NA: {r}"
        );
    }
}

/// Pre-existing crash: an ACF threshold no read clears used to panic on the
/// empty reservoir (x.len() - 1 underflow). Must exit 0 with no acf rows.
#[test]
fn qc_acf_empty_reservoir_exits_cleanly() {
    let out = run(&[
        "qc",
        "--acf",
        "--acf-min-m6a",
        "1000000",
        fixture("all.bam").to_str().unwrap(),
    ]);
    assert!(
        !out.lines().any(|l| l.starts_with("m6a_acf")),
        "no m6a_acf rows on an empty reservoir"
    );
}

/// Seeded reservoir replacement is reproducible; the passing subset is a
/// subset of the sample.
#[test]
fn qc_acf_seeded_replacement() {
    let tagged = tagged_bam("three_two.bam");
    let args = [
        "qc",
        "--acf",
        "--acf-max-reads",
        "3",
        "--acf-sample-rate",
        "1",
        "--seed",
        "42",
        tagged.path().to_str().unwrap(),
    ];
    let out = run(&args);
    assert_eq!(out, run(&args));
    let f: Vec<&str> = out
        .lines()
        .find(|l| l.starts_with("acf_reads"))
        .unwrap()
        .split('\t')
        .collect();
    assert_eq!(f[2], "3");
    assert!(f[3].parse::<i64>().unwrap() <= 3);
}

/// Under custom minimums, the fiberseq_callable state rows and the
/// filtered column must agree: the loop re-ensures with the real filters
/// because the stream runs on zeroed ones.
#[test]
fn qc_custom_minimums_reach_state_rows() {
    let tagged = tagged_bam("all.bam");
    let out = run(&["qc", "--min-msp", "100000", tagged.path().to_str().unwrap()]);
    let fc: Vec<(String, i64, i64)> = out
        .lines()
        .filter(|l| l.starts_with("fiberseq_callable"))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (
                f[1].to_string(),
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
            )
        })
        .collect();
    let get = |k: &str| fc.iter().find(|r| r.0 == k).unwrap();
    assert_eq!(get("Callable").1, 0, "no read clears minimum 100000");
    assert_eq!(get("NotCallable").1, 22, "state rows use the user minimums");
    let filt: i64 = rows_for(&out, "phased_reads").iter().map(|r| r.2).sum();
    assert_eq!(filt, 0, "filtered column agrees with the statet rows");
}
