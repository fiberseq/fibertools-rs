use super::fiber::FiberseqData;
use super::*;
use rayon::prelude::*;
use rust_htslib::{bam, bam::HeaderView, bam::Read};

pub fn process_bam_chunk(
    records: Vec<bam::Record>,
    out_files: &mut FiberOut,
    head_view: &HeaderView,
) {
    let fiber_data = FiberseqData::from_records(records, head_view, out_files.min_ml_score);

    match &mut out_files.m6a {
        Some(m6a) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_m6a(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, m6a);
            }
        }
        None => {}
    }
    match &mut out_files.cpg {
        Some(cpg) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_cpg(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, cpg);
            }
        }
        None => {}
    }
    match &mut out_files.msp {
        Some(msp) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_msp(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, msp);
            }
        }
        None => {}
    }
    match &mut out_files.nuc {
        Some(nuc) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_nuc(out_files.reference))
                .collect();
            for line in out {
                write_to_file(&line, nuc);
            }
        }
        None => {}
    }
    match &mut out_files.all {
        Some(all) => {
            let out: Vec<String> = fiber_data
                .par_iter()
                .map(|r| r.write_all(out_files.simplify, out_files.quality, out_files.full_float))
                .collect();
            for line in out {
                write_to_file(&line, all);
            }
        }
        None => {}
    }
}

pub fn extract_contained(bam: &mut bam::Reader, mut out_files: FiberOut) {
    let header = bam::Header::from_template(bam.header());
    let head_view = bam::HeaderView::from_header(&header);

    // print the header if in all mode
    match &mut out_files.all {
        Some(all) => {
            write!(
                all,
                "{}",
                FiberseqData::all_header(out_files.simplify, out_files.quality)
            )
            .unwrap();
        }
        None => {}
    }

    // process bam in chunks
    // keeps mem pretty low, about 1GB per thread
    let bam_chunk_iter = BamChunk::new(bam.records(), None);
    for chunk in bam_chunk_iter {
        process_bam_chunk(chunk, &mut out_files, &head_view);
    }
}
