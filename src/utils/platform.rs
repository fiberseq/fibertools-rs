//! Detect the sequencing platform (PacBio vs ONT) per read, so downstream
//! commands can auto-select models/parameters. Aux tags are checked first;
//! when they've been stripped, the read name is used as a fallback.

use lazy_static::lazy_static;
use regex::Regex;
use rust_htslib::bam::Record;

/// PacBio-only per-read aux tags: `zm` ZMW, `np` passes, `rq` read quality,
/// `sn` SNR. Never emitted by ONT or re-used by fibertools.
/// See <https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html#use-of-read-tags-for-per-read-information>.
const PACBIO_READ_TAGS: [&[u8; 2]; 4] = [b"zm", b"np", b"rq", b"sn"];

/// ONT-only per-read aux tags: `mx` mux, `ch` channel, `st` start time,
/// `du` duration, `ts` trimmed samples. Excludes `ns`/`nl`/`fn`/`rn`, which
/// fibertools re-uses for its own annotation arrays.
/// See <https://nanoporetech.github.io/ont-output-specifications/latest/read_formats/bam/#read-tags>.
const ONT_READ_TAGS: [&[u8; 2]; 5] = [b"mx", b"ch", b"st", b"du", b"ts"];

lazy_static! {
    /// PacBio read names are `<movie>/<zmw>/<type>`, e.g.
    /// `m84081_231207_200206_s1/240456095/ccs`.
    /// See <https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html#qname-convention>.
    static ref PACBIO_QNAME: Regex = Regex::new(r"^m\w+/\d+/").unwrap();
    /// ONT read names are the per-read UUID (`read_id`), e.g.
    /// `d78c9b31-cec7-4381-bc1c-5f1513d070b9`.
    /// See <https://nanoporetech.github.io/ont-output-specifications/latest/read_formats/bam/>.
    static ref ONT_QNAME: Regex = Regex::new(
        r"^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}-[0-9a-fA-F]{12}$"
    ).unwrap();
}

/// The sequencing platform a single read originated from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeqPlatform {
    PacBio,
    Ont,
    /// No usable signal on this read.
    Unknown,
}

/// Detect the platform of a single read: aux tags first, then the read name as
/// a fallback for reads whose basecaller tags have been stripped.
pub fn platform_from_record(rec: &Record) -> SeqPlatform {
    match platform_from_aux(rec) {
        SeqPlatform::Unknown => platform_from_qname(rec.qname()),
        known => known,
    }
}

/// Detect from per-read aux tag signatures. Presence is a positive signal only;
/// absence means nothing (many tags are optional/flag-gated).
pub fn platform_from_aux(rec: &Record) -> SeqPlatform {
    let has_any = |tags: &[&[u8; 2]]| tags.iter().any(|t| rec.aux(&t[..]).is_ok());
    let pacbio = has_any(&PACBIO_READ_TAGS);
    let ont = has_any(&ONT_READ_TAGS);
    match (pacbio, ont) {
        (true, false) => SeqPlatform::PacBio,
        (false, true) => SeqPlatform::Ont,
        // Neither, or contradictory: don't guess.
        _ => SeqPlatform::Unknown,
    }
}

/// Detect from the read-name convention: PacBio `<movie>/<zmw>/...` vs an ONT
/// UUID. Survives aux-tag stripping.
pub fn platform_from_qname(qname: &[u8]) -> SeqPlatform {
    let name = String::from_utf8_lossy(qname);
    if PACBIO_QNAME.is_match(&name) {
        SeqPlatform::PacBio
    } else if ONT_QNAME.is_match(&name) {
        SeqPlatform::Ont
    } else {
        SeqPlatform::Unknown
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{self, Read};

    fn first_record(path: &str) -> Record {
        let mut reader = bam::Reader::from_path(path).unwrap();
        let mut rec = Record::new();
        reader
            .read(&mut rec)
            .expect("fixture has at least one record")
            .expect("read ok");
        rec
    }

    const PACBIO_FIXTURES: [&str; 5] = [
        "tests/data/revio.bam",
        "tests/data/two_two.bam",
        "tests/data/three_two.bam",
        "tests/data/all.bam",
        "tests/data/ctcf.bam",
    ];
    const ONT_FIXTURE: &str = "tests/data/ONT.NAPA.bam";

    #[test]
    fn aux_detects_pacbio() {
        for f in PACBIO_FIXTURES {
            assert_eq!(
                platform_from_aux(&first_record(f)),
                SeqPlatform::PacBio,
                "{f}"
            );
        }
    }

    #[test]
    fn aux_detects_ont() {
        assert_eq!(
            platform_from_aux(&first_record(ONT_FIXTURE)),
            SeqPlatform::Ont
        );
    }

    #[test]
    fn qname_detects_pacbio() {
        for f in PACBIO_FIXTURES {
            assert_eq!(
                platform_from_qname(first_record(f).qname()),
                SeqPlatform::PacBio,
                "{f}"
            );
        }
    }

    #[test]
    fn qname_detects_ont() {
        assert_eq!(
            platform_from_qname(first_record(ONT_FIXTURE).qname()),
            SeqPlatform::Ont
        );
    }

    #[test]
    fn qname_recovers_platform_after_aux_stripping() {
        // Strip the discriminator aux tags; the read name must still classify.
        let mut rec = first_record(ONT_FIXTURE);
        for t in ONT_READ_TAGS {
            rec.remove_aux(t).ok();
        }
        assert_eq!(platform_from_aux(&rec), SeqPlatform::Unknown);
        assert_eq!(platform_from_record(&rec), SeqPlatform::Ont);

        let mut rec = first_record("tests/data/revio.bam");
        for t in PACBIO_READ_TAGS {
            rec.remove_aux(t).ok();
        }
        assert_eq!(platform_from_aux(&rec), SeqPlatform::Unknown);
        assert_eq!(platform_from_record(&rec), SeqPlatform::PacBio);
    }

    #[test]
    fn qname_patterns_and_non_matches() {
        assert_eq!(
            platform_from_qname(b"m84081_231207_200206_s1/240456095/ccs"),
            SeqPlatform::PacBio
        );
        assert_eq!(
            platform_from_qname(b"m54329U_210810_004956/120523046/0_5000"),
            SeqPlatform::PacBio
        );
        assert_eq!(
            platform_from_qname(b"d78c9b31-cec7-4381-bc1c-5f1513d070b9"),
            SeqPlatform::Ont
        );
        assert_eq!(platform_from_qname(b"read_12345"), SeqPlatform::Unknown);
    }
}
