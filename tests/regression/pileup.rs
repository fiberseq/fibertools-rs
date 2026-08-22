use super::common::{fixture, run, select_tsv_cols, tagged_bam};
use tempfile::NamedTempFile;

#[test]
fn pileup_default() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
            "msp_coverage",
        ]
    ));
}

#[test]
fn pileup_m6a() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--m6a",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
            "msp_coverage",
            "m6a_coverage",
        ]
    ));
}

#[test]
fn pileup_no_msp() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--no-msp",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "nuc_coverage",
        ]
    ));
}

#[test]
fn pileup_no_nuc() {
    let tmp = NamedTempFile::new().unwrap();
    run(&[
        "pileup",
        fixture("all.bam").to_str().unwrap(),
        "--no-nuc",
        "-o",
        tmp.path().to_str().unwrap(),
    ]);
    let out = std::fs::read_to_string(tmp.path()).unwrap();
    insta::assert_snapshot!(select_tsv_cols(
        &out,
        &[
            "#chrom",
            "start",
            "end",
            "coverage",
            "fire_coverage",
            "score",
            "msp_coverage",
        ]
    ));
}

fn index(path: &std::path::Path) {
    rust_htslib::bam::index::build(path, None, rust_htslib::bam::index::Type::Bai, 1)
        .expect("index tagged bam");
}

/// --callable-fibers off must be inert: pileup on a tagged BAM is
/// byte-identical to the untagged input (an accidental default flip would
/// silently shrink every cached bigWig).
#[test]
fn pileup_callable_fibers_default_is_inert() {
    let tagged = tagged_bam("all.bam");
    index(tagged.path());
    let pile = |bam: &str| -> String {
        let tmp = NamedTempFile::new().unwrap();
        run(&["pileup", bam, "-o", tmp.path().to_str().unwrap()]);
        std::fs::read_to_string(tmp.path()).unwrap()
    };
    // Build a copy of the tagged BAM with ONLY the callable section
    // removed, so the two inputs differ in nothing else. With the flag
    // off, pileup must not read the tag: outputs are byte-identical.
    use fibertools_rs::utils::ma_io::{read_record, write_record, FIBERSEQ_CALLABLE_TYPE};
    use rust_htslib::bam::Read;
    let untagged = NamedTempFile::with_suffix(".bam").unwrap();
    {
        let mut reader = rust_htslib::bam::Reader::from_path(tagged.path()).unwrap();
        let header = rust_htslib::bam::Header::from_template(reader.header());
        let mut writer = rust_htslib::bam::Writer::from_path(
            untagged.path(),
            &header,
            rust_htslib::bam::Format::Bam,
        )
        .unwrap();
        for rec in reader.records() {
            let mut record = rec.unwrap();
            let mut annot = read_record(&record).unwrap();
            annot
                .annotation_types
                .retain(|t| t.name != FIBERSEQ_CALLABLE_TYPE);
            write_record(&mut record, &annot);
            writer.write(&record).unwrap();
        }
    }
    index(untagged.path());
    let with_tag = pile(tagged.path().to_str().unwrap());
    let without_tag = pile(untagged.path().to_str().unwrap());
    assert!(!with_tag.is_empty());
    assert_eq!(with_tag, without_tag, "flag off must never read the tag");
}

/// With the flag on, coverage shrinks; and combined with --fire-filter the
/// span can only narrow further (intersection semantics, never widen).
#[test]
fn pileup_callable_fibers_shrinks_and_intersects() {
    let tagged = tagged_bam("three_two.bam");
    index(tagged.path());
    let total_cov = |extra: &[&str]| -> i64 {
        let tmp = NamedTempFile::new().unwrap();
        let mut args = vec!["pileup", tagged.path().to_str().unwrap()];
        args.extend_from_slice(extra);
        args.extend_from_slice(&["-o", tmp.path().to_str().unwrap()]);
        run(&args);
        let out = std::fs::read_to_string(tmp.path()).unwrap();
        let mut lines = out.lines();
        let hdr: Vec<&str> = lines.next().unwrap().split('\t').collect();
        let (s, e, c) = (
            hdr.iter().position(|h| *h == "start").unwrap(),
            hdr.iter().position(|h| *h == "end").unwrap(),
            hdr.iter().position(|h| *h == "coverage").unwrap(),
        );
        lines
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                let (st, en): (i64, i64) = (f[s].parse().unwrap(), f[e].parse().unwrap());
                let cov: i64 = f[c].parse().unwrap();
                (en - st) * cov
            })
            .sum()
    };
    let plain = total_cov(&[]);
    let callable = total_cov(&["--callable-fibers"]);
    let fiber = total_cov(&["--fiber-coverage"]);
    let fire = total_cov(&["--fire-coverage"]);
    let fire_filter = total_cov(&["--fire-filter"]);
    assert!(callable < plain, "the flag must shrink coverage");
    assert_eq!(callable, fiber, "--fiber-coverage is the same flag");
    assert_eq!(callable, fire, "--fire-coverage is the same flag");
    assert_eq!(callable, fire_filter, "--fire-filter is the same flag");
}
