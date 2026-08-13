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

// `msp_nuc.bam` carries legacy nuc/msp tags (ns/nl/as/al/aq) and MM/ML base
// mods. convert-tags should: drop the legacy tags, emit an MA tag, and leave
// MM/ML byte-identical (it never touches base modifications).
#[test]
fn convert_tags_rewrites_legacy_as_ma() {
    let input = fixture("msp_nuc.bam");
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
    let input = fixture("msp_nuc.bam");
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

// Legacy fibertig fs/fl/fa tags are consumed by the reader (into the MA
// fibertig type) and must be stripped on conversion just like ns/nl/as/al/aq
// — otherwise the source tags survive as duplicate, potentially-stale copies.
#[test]
fn convert_tags_strips_consumed_fibertig_tags() {
    use rust_htslib::bam::Header;

    // synthesize a legacy fibertig BAM from the first msp_nuc.bam record
    let input = fixture("msp_nuc.bam");
    let mut reader = bam::Reader::from_path(&input).unwrap();
    let header = Header::from_template(reader.header());
    let synth = NamedTempFile::with_suffix(".bam").unwrap();
    {
        let mut writer = bam::Writer::from_path(synth.path(), &header, bam::Format::Bam).unwrap();
        let mut rec = reader.records().next().unwrap().unwrap();
        rec.push_aux(b"fs", Aux::ArrayU32((&vec![100u32, 500]).into()))
            .unwrap();
        rec.push_aux(b"fl", Aux::ArrayU32((&vec![50u32, 60]).into()))
            .unwrap();
        rec.push_aux(b"fa", Aux::String("gene_a|")).unwrap();
        writer.write(&rec).unwrap();
    }

    let out = NamedTempFile::with_suffix(".bam").unwrap();
    convert(synth.path(), out.path());
    let rec = &records(out.path())[0];
    for tag in [b"fs", b"fl", b"fa"] {
        assert!(
            rec.aux(tag).is_err(),
            "fibertig tag {} survived conversion",
            String::from_utf8_lossy(tag)
        );
    }
    let ma_tag = ma(rec).expect("MA tag not written");
    assert!(
        ma_tag.contains("fibertig"),
        "fibertig annotations not carried into MA tag: {ma_tag}"
    );
    let an = match rec.aux(b"AN") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => String::new(),
    };
    assert!(
        an.contains("gene_a"),
        "fibertig name not carried into AN tag: {an}"
    );
}
