use super::common::{ft, run};
use tempfile::{NamedTempFile, TempPath};

/// s3 has two overlapping intervals (20000-20200, 20100-20400) so the tests can pin that
/// one file only ever contributes 1 to a peak's support count.
fn fixture() -> [TempPath; 3] {
    let beds = [
        "chr1\t1000\t1200\nchr1\t5000\t5300\nchr1\t9000\t9100\nchr2\t100\t400\n",
        "chr1\t1020\t1220\nchr1\t5050\t5350\nchr2\t150\t420\n",
        "chr1\t1050\t1250\nchr1\t20000\t20200\nchr1\t20100\t20400\nchr2\t120\t380\n",
    ];
    beds.map(|contents| {
        let bed = NamedTempFile::with_suffix(".bed").unwrap();
        std::fs::write(bed.path(), contents).unwrap();
        bed.into_temp_path()
    })
}

fn union_peaks(beds: &[TempPath], extra: &[&str]) -> String {
    let mut args = vec!["union-peaks"];
    args.extend(beds.iter().map(|b| b.to_str().unwrap()));
    args.extend(["--names", "s1,s2,s3"]);
    args.extend(extra);
    run(&args)
}

// Peak boundaries are the peak caller's median boundaries, not the outer span of the
// supporting intervals, and support is counted per input file. Peak 4 is the collapse
// case: s3's two overlapping intervals must yield one peak with support 1, not two peaks
// or a support of 2.
#[test]
fn union_peaks_support_counts_and_boundaries() {
    let beds = fixture();
    assert_eq!(
        union_peaks(&beds, &[]),
        "#chrom\tstart\tend\tname\tn_support\tfrac_support\tsupport\tunion_start\tunion_end\tpeak_summit\n\
         chr1\t1020\t1220\tunion_peak_1\t3\t1.0000\ts1,s2,s3\t1000\t1250\t1125\n\
         chr1\t5050\t5350\tunion_peak_2\t2\t0.6667\ts1,s2\t5000\t5350\t5175\n\
         chr1\t9000\t9100\tunion_peak_3\t1\t0.3333\ts1\t9000\t9100\t9050\n\
         chr1\t20000\t20400\tunion_peak_4\t1\t0.3333\ts3\t20000\t20400\t20200\n\
         chr2\t120\t400\tunion_peak_5\t3\t1.0000\ts1,s2,s3\t100\t420\t265\n"
    );
}

// --min-support is an output filter, so it must drop exactly the low support rows and
// leave the surviving rows byte for byte the same as the unfiltered run.
#[test]
fn union_peaks_min_support_filters_only() {
    let beds = fixture();
    let all = union_peaks(&beds, &[]);
    let filtered = union_peaks(&beds, &["-n", "2"]);
    let kept: Vec<&str> = filtered.lines().skip(1).collect();
    assert_eq!(kept.len(), 3, "got: {filtered}");
    for (line, expected) in kept.iter().zip(["1020\t1220", "5050\t5350", "120\t400"]) {
        assert!(line.contains(expected), "{line} lacks {expected}");
        // same peak, only the sequential name changes
        let coords = line.split('\t').take(3).collect::<Vec<_>>().join("\t");
        assert!(all.contains(&coords), "{coords} not in unfiltered output");
    }
}

// mock-fire writes in read-name order, so a hand rolled `mock-fire | samtools index |
// call-peaks` pipeline fails with "unsorted positions" whenever sample names sort against
// their positions. union-peaks never writes a BAM, so it has to work here.
#[test]
fn union_peaks_names_ordered_against_positions() {
    let beds: [TempPath; 2] = [
        "chr1\t100\t200\nchr10\t120\t220\nchr2\t150\t250\n",
        "chr1\t5000\t5100\nchr10\t100\t200\nchr2\t100\t200\n",
    ]
    .map(|contents| {
        let bed = NamedTempFile::with_suffix(".bed").unwrap();
        std::fs::write(bed.path(), contents).unwrap();
        bed.into_temp_path()
    });
    let out = run(&[
        "union-peaks",
        beds[0].to_str().unwrap(),
        beds[1].to_str().unwrap(),
        "--names",
        "zz,aa",
    ]);
    let chroms: Vec<&str> = out
        .lines()
        .skip(1)
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    assert_eq!(chroms, ["chr1", "chr1", "chr10", "chr2"], "got: {out}");
    assert_eq!(out.matches("zz,aa").count(), 2, "got: {out}");
}

// An empty BED is a peak caller that found nothing, not an error, but it still counts in
// the frac_support denominator.
#[test]
fn union_peaks_empty_input_is_a_sample_with_no_peaks() {
    let beds = fixture();
    let empty = NamedTempFile::with_suffix(".bed").unwrap();
    let out = run(&[
        "union-peaks",
        beds[0].to_str().unwrap(),
        empty.path().to_str().unwrap(),
    ]);
    assert_eq!(out.lines().count(), 5, "got: {out}");
    assert!(
        out.lines().skip(1).all(|l| l.contains("\t0.5000\t")),
        "{out}"
    );
}

#[test]
fn union_peaks_rejects_bad_input() {
    let beds = fixture();
    let bad = NamedTempFile::with_suffix(".bed").unwrap();
    std::fs::write(bad.path(), "chr1\t500\t100\n").unwrap();
    let cases: Vec<(Vec<&str>, &str)> = vec![
        (
            vec![beds[0].to_str().unwrap(), "--names", "only-one,and,three"],
            "--names has 3 names",
        ),
        (vec![bad.path().to_str().unwrap()], "invalid interval"),
        (
            vec![beds[0].to_str().unwrap(), "--window-size", "1"],
            "--window-size must be at least 2",
        ),
    ];
    for (args, expected) in cases {
        let out = std::process::Command::new(ft())
            .arg("union-peaks")
            .args(&args)
            .output()
            .unwrap();
        assert!(!out.status.success(), "{args:?} unexpectedly succeeded");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(stderr.contains(expected), "{args:?} stderr: {stderr}");
    }
}
