use crate::cli::DddaToM6aOptions;
use crate::utils::basemods;
use crate::utils::ma_io;
use crate::*;
use anyhow::Error;
use bio::alphabets::dna::revcomp;
use molecular_annotation::{MolecularAnnotations, QualitySpec, Strand};
use rayon::iter::ParallelIterator;
use rayon::prelude::*;
use rust_htslib::bam::Record;

pub fn ddda_to_m6a_record(record: &mut Record, _opts: &DddaToM6aOptions) {
    let was_reverse = record.is_reverse();
    // clear the flag of any reverse or positive strand, by resetting the 4th bit, aka 16, to zero
    // let mut flag = record.flags() & !(1 << 4);
    // save some values we will need to modify
    let cigar = bam::record::CigarString(record.cigar().to_vec());
    let q_name = record.qname().to_vec();
    let qual = record.qual().to_vec();
    let mut forward_seq = record.seq().as_bytes();
    if was_reverse {
        forward_seq = revcomp(&forward_seq);
    }

    // Positions to mark as m6a, paired with the canonical post-Y/R base
    // they land on (T for Y, A for R) in molecular orientation. We carry
    // the base alongside the position so the write path knows the right
    // MM group header without re-indexing the (possibly revcomp'd) seq.
    let mut modified_bases_forward: Vec<(u32, u8)> = vec![];
    let mut y_count = 0;
    let mut r_count = 0;
    let mut new_forward_seq = vec![];
    for (idx, bp) in forward_seq.iter().enumerate() {
        if bp == &b'Y' {
            modified_bases_forward.push((idx as u32, b'T'));
            y_count += 1;
            new_forward_seq.push(b'T');
        } else if bp == &b'R' {
            modified_bases_forward.push((idx as u32, b'A'));
            r_count += 1;
            new_forward_seq.push(b'A');
        } else {
            new_forward_seq.push(*bp);
        }
    }

    // validate
    if !(y_count == 0 || r_count == 0) {
        panic!("Y and R cannot be in the same sequence");
    } else if y_count == 0 && r_count == 0 {
        return;
    }

    if was_reverse {
        // must reverse complement the new sequence if we are on the bottom strand
        new_forward_seq = revcomp(&new_forward_seq);
    }

    // set values for the new modified record
    record.set(&q_name, Some(&cigar), &new_forward_seq, &qual);

    // Load existing MM/ML into a MolecularAnnotations, drop any
    // pre-existing m6a (we're replacing it with the Y/R-derived calls),
    // then append the synthesized m6a and write back. Each call gets
    // its canonical MM group header (A+a or T-a) tagged via
    // `canonical_header` so the library's MM/ML serializer emits the
    // right groups.
    let mut annot = ma_io::read_record(record).unwrap_or_else(|e| {
        log::warn!("read_record failed for {:?}: {e}", String::from_utf8_lossy(record.qname()));
        MolecularAnnotations::from_record(record)
    });
    annot
        .annotation_types
        .retain(|t| t.name != basemods::M6A_TYPE);
    let qspec = "Q".parse::<QualitySpec>().expect("Q parses");
    let t = annot.add_annotation_type(basemods::M6A_TYPE, qspec);
    for (pos, base) in modified_bases_forward {
        let header = basemods::canonical_header(basemods::M6A_TYPE, base)
            .expect("ddda_to_m6a only places calls on A/T bases");
        t.add(pos, 1, Strand::Forward, vec![255], Some(header.to_string()));
    }
    ma_io::ensure_basemod_encoding(&mut annot);
    ma_io::write_record(record, &annot);
}

pub fn ddda_to_m6a(opts: &mut DddaToM6aOptions) -> Result<(), Error> {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);

    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);
    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        // covert the ddda record to an m6a one
        chunk
            .par_iter_mut()
            .for_each(|record| ddda_to_m6a_record(record, opts));
        // write the output
        chunk
            .into_iter()
            .for_each(|record| out.write(&record).unwrap());
    }
    Ok(())
}
