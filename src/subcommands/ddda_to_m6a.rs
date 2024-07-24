use crate::cli::DddaToM6aOptions;
use crate::utils::basemods;
use crate::utils::basemods::BaseMods;
use crate::*;
use anyhow::Error;
use bio::alphabets::dna::revcomp;
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

    // get the bases that will be modified
    let mut modified_bases_forward = vec![];
    let mut y_count = 0;
    let mut r_count = 0;
    let mut new_forward_seq = vec![];
    for (idx, bp) in forward_seq.iter().enumerate() {
        if bp == &b'Y' {
            modified_bases_forward.push(idx as i64);
            y_count += 1;
            new_forward_seq.push(b'T');
        } else if bp == &b'R' {
            modified_bases_forward.push(idx as i64);
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

    // check if we should set the read to the top or bottom strand
    let is_top_strand = y_count > 0;

    // set things up for the bottom strand
    // let is_bottom_strand = !is_top_strand;
    //if is_bottom_strand {
    //    flag |= 1 << 4;
    //    record.set_flags(flag);
    //    record.set_reverse();
    //}
    if was_reverse {
        // must reverse complement the new sequence if we are on the bottom strand
        new_forward_seq = revcomp(&new_forward_seq);
    }

    // set values for the new modified record
    record.set(&q_name, Some(&cigar), &new_forward_seq, &qual);

    // set up the fake base mods
    let (modified_base, strand) = if is_top_strand {
        (b'T', '-')
    } else {
        (b'A', '+')
    };
    let modification_type = 'a';
    let modified_probabilities_forward = vec![255; modified_bases_forward.len()];

    let fake_base_mods = basemods::BaseMod::new(
        record,
        modified_base,
        strand,
        modification_type,
        modified_bases_forward,
        modified_probabilities_forward,
    );

    // add to the record
    let mut base_mods = basemods::BaseMods::new(record, 0);
    base_mods.base_mods.push(fake_base_mods);
    base_mods.add_mm_and_ml_tags(record);

    // validates that the base mods were added correctly
    let _ = BaseMods::new(record, 0);
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
