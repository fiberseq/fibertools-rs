//! Repro: fibertig peaks on a *split* record vanish once the record is
//! flipped to the reverse strand (what `vg surject` does when a fibertig
//! lands on the reverse strand of the target).
use anyhow::Result;
use fibertools_rs::utils::fibertig::{FiberTig, FIBERTIG_TYPE};
use fibertools_rs::utils::ma_io;
use molecular_annotation::{Encoding, MolecularAnnotations, QualitySpec, Strand};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::HeaderView;
use std::io::Write;
use tempfile::NamedTempFile;

/// Build split records exactly as `ft pg-inject --split-size` does.
fn split_records(contig_len: usize, split_size: i64) -> Result<Vec<rust_htslib::bam::Record>> {
    let mut fasta_file = NamedTempFile::new()?;
    writeln!(fasta_file, ">chr1")?;
    writeln!(fasta_file, "{}", "ACGT".repeat(contig_len / 4))?;
    fasta_file.flush()?;
    let sequences = FiberTig::read_fasta_into_vec(fasta_file.path().to_str().unwrap())?;
    let header = FiberTig::create_mock_bam_header_from_sequences(&sequences);
    let header_view = HeaderView::from_header(&header);

    // Peaks: one per window, as read_bed_annotations would build them.
    let mut annot = MolecularAnnotations::new(contig_len as u32);
    {
        let t = annot.add_annotation_type(FIBERTIG_TYPE, QualitySpec::none(), Encoding::Ma);
        t.add(100, 50, Strand::Forward, vec![], Some("p1".into()));
        t.add(1_100, 50, Strand::Forward, vec![], Some("p2".into()));
        t.add(2_100, 50, Strand::Forward, vec![], Some("p3".into()));
    }
    let mut splits = FiberTig::approximately_divide_annotations_by_window_size(
        contig_len as i64,
        split_size,
        &annot,
    );
    FiberTig::create_annotated_records_from_splits(
        "chr1",
        &sequences[0].1,
        &mut splits,
        &header_view,
    )
}

/// Simulate vg surject landing this record on the reverse strand: set flag
/// 0x10, revcomp SEQ. CIGAR stays a single `=` so aligned blocks are identity.
/// Aux tags (including Ma/An) are kept verbatim, as vg does.
fn flip_to_reverse(rec: &rust_htslib::bam::Record) -> rust_htslib::bam::Record {
    let mut out = rec.clone();
    let seq = bio::alphabets::dna::revcomp(rec.seq().as_bytes());
    let qual = rec.qual().to_vec();
    let cigar = rec.cigar().take();
    out.set(rec.qname(), Some(&cigar), &seq, &qual);
    out.set_reverse();
    out
}

#[test]
fn split_record_ma_read_length_matches_seq_len() -> Result<()> {
    let records = split_records(4_000, 1_000)?;
    assert!(records.len() > 1, "expected splitting to happen");
    for rec in &records {
        let Ok(Aux::String(ma)) = rec.aux(b"Ma") else {
            panic!("no Ma tag")
        };
        let read_len: usize = ma.split(';').next().unwrap().parse()?;
        eprintln!("pos={} seq_len={} Ma={ma}", rec.pos(), rec.seq_len());
        assert_eq!(
            read_len,
            rec.seq_len(),
            "Ma read_length ({read_len}) != SEQ length ({}) for record at pos {}",
            rec.seq_len(),
            rec.pos()
        );
    }
    Ok(())
}

#[test]
fn split_record_reverse_strand_keeps_peaks() -> Result<()> {
    let records = split_records(4_000, 1_000)?;
    assert!(records.len() > 1, "expected splitting to happen");
    for rec in &records {
        let fwd: Vec<_> = ma_io::read_record(rec)?
            .iter_type(FIBERTIG_TYPE)
            .into_iter()
            .flatten()
            .map(|i| (i.ref_start, i.ref_end, i.name.map(String::from)))
            .collect();
        let rev_rec = flip_to_reverse(rec);
        let rev: Vec<_> = ma_io::read_record(&rev_rec)?
            .iter_type(FIBERTIG_TYPE)
            .into_iter()
            .flatten()
            .map(|i| {
                (
                    i.query_start,
                    i.query_end,
                    i.ref_start,
                    i.ref_end,
                    i.name.map(String::from),
                )
            })
            .collect();
        eprintln!(
            "pos={} seq_len={} fwd={fwd:?} rev={rev:?}",
            rec.pos(),
            rec.seq_len()
        );
        if fwd.is_empty() {
            continue;
        }
        // Identity alignment: forward ref [s, e) must mirror to
        // [seq_len - e, seq_len - s) on the reverse strand, same name.
        // Ref coords include the record's pos, so mirror within the record.
        let pos = rec.pos() as u32;
        let n = rec.seq_len() as u32;
        let expected: Vec<_> = fwd
            .iter()
            .map(|(rs, re, name)| {
                (
                    Some(pos + n - (re.unwrap() - pos)),
                    Some(pos + n - (rs.unwrap() - pos)),
                    name.clone(),
                )
            })
            .collect();
        let got: Vec<_> = rev
            .iter()
            .map(|(_, _, rs, re, name)| (*rs, *re, name.clone()))
            .collect();
        assert_eq!(
            got,
            expected,
            "reverse-strand peaks wrong on record at pos {}",
            rec.pos()
        );
    }
    Ok(())
}

#[test]
fn unsplit_record_reverse_strand_keeps_peaks() -> Result<()> {
    // Control: no splitting (split_size >= contig) -> read_length is correct.
    let records = split_records(4_000, 1_000_000)?;
    assert_eq!(records.len(), 1);
    let rev_rec = flip_to_reverse(&records[0]);
    let rev: Vec<_> = ma_io::read_record(&rev_rec)?
        .iter_type(FIBERTIG_TYPE)
        .into_iter()
        .flatten()
        .map(|i| (i.ref_start, i.ref_end))
        .collect();
    eprintln!("unsplit rev={rev:?}");
    assert!(rev.iter().all(|(rs, re)| rs.is_some() && re.is_some()));
    Ok(())
}
