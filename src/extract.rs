use super::fiber::FiberseqData;
use super::*;
use colored::Colorize;
use rayon::{current_num_threads, prelude::*};
use rust_htslib::{bam, bam::HeaderView, bam::Read};
use std::time::Instant;

pub fn process_bam_chunk(
    records: Vec<bam::Record>,
    so_far: usize,
    out_files: &mut FiberOut,
    head_view: &HeaderView,
) {
    let start = Instant::now();
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

    let duration = start.elapsed().as_secs_f64();
    log::info!(
        "Processing {}, {} reads done so far.",
        format!("{:.2?} reads/s", fiber_data.len() as f64 / duration)
            .bright_cyan()
            .bold(),
        format!("{:}", so_far + fiber_data.len())
            .bright_magenta()
            .bold()
    );
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
    let chunk_size = current_num_threads() * 500;
    let bam_chunk_iter = BamChunk {
        bam: bam.records(),
        chunk_size,
    };
    let mut processed_reads = 0;
    for chunk in bam_chunk_iter {
        processed_reads += chunk.len();
        process_bam_chunk(chunk, processed_reads, &mut out_files, &head_view);
    }
}
