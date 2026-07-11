use super::common::{fixture, ft};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, Read};
use std::process::Command;
use tempfile::NamedTempFile;

const LEGACY_TAGS: &[&[u8]] = &[b"ns", b"nl", b"as", b"al", b"aq"];

fn records(path: &std::path::Path) -> Vec<bam::Record> {
    let mut reader = bam::Reader::from_path(path).unwrap();
    reader.records().map(|r| r.unwrap()).collect()
}

fn convert(input: &std::path::Path, output: &std::path::Path) {
    let status = Command::new(ft())
        .args([
            "convert-tags",
            input.to_str().unwrap(),
            output.to_str().unwrap(),
        ])
        .status()
        .expect("spawn ft convert-tags");
    assert!(status.success(), "convert-tags exited {status}");
}

fn mm(rec: &bam::Record) -> Option<String> {
    match rec.aux(b"MM") {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

fn ml(rec: &bam::Record) -> Option<Vec<u8>> {
    match rec.aux(b"ML") {
        Ok(Aux::ArrayU8(arr)) => Some(arr.iter().collect()),
        _ => None,
    }
}

fn ma(rec: &bam::Record) -> Option<String> {
    match rec.aux(b"MA") {
        Ok(Aux::String(s)) => Some(s.to_string()),
        _ => None,
    }
}

// `all.bam` carries legacy nuc/msp tags (ns/nl/as/al/aq) and MM/ML base mods.
// convert-tags should: drop the legacy tags, emit an MA tag, and leave MM/ML
// byte-identical (it never touches base modifications).
#[test]
fn convert_tags_rewrites_legacy_as_ma() {
    let input = fixture("all.bam");
    let out = NamedTempFile::new().unwrap();
    convert(&input, out.path());

    let before = records(&input);
    let after = records(out.path());
    assert_eq!(before.len(), after.len(), "record count changed");

    for (b, a) in before.iter().zip(&after) {
        for tag in LEGACY_TAGS {
            assert!(
                a.aux(tag).is_err(),
                "legacy tag {} survived conversion",
                String::from_utf8_lossy(tag)
            );
        }
        assert!(ma(a).is_some(), "MA tag not written");
        assert_eq!(mm(b), mm(a), "MM changed");
        assert_eq!(ml(b), ml(a), "ML changed");
    }
}

// Converting an already-MA BAM is a no-op on the annotation tags: the MA tag is
// unchanged and no legacy tags are (re)introduced.
#[test]
fn convert_tags_is_idempotent() {
    let input = fixture("all.bam");
    let once = NamedTempFile::new().unwrap();
    let twice = NamedTempFile::new().unwrap();
    convert(&input, once.path());
    convert(once.path(), twice.path());

    let first = records(once.path());
    let second = records(twice.path());
    assert_eq!(first.len(), second.len());

    for (a, b) in first.iter().zip(&second) {
        assert_eq!(ma(a), ma(b), "MA tag changed on re-convert");
        for tag in LEGACY_TAGS {
            assert!(
                b.aux(tag).is_err(),
                "legacy tag {} reintroduced on re-convert",
                String::from_utf8_lossy(tag)
            );
        }
    }
}
