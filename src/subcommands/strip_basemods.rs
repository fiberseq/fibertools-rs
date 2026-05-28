use crate::cli::StripBasemodsOptions;
use crate::utils::bamannotations::primary_qual;
use crate::utils::basemods::{CPG_TYPE, M6A_TYPE};
use crate::utils::bio_io::BamChunk;
use crate::utils::ma_io;
use bio::alphabets::dna::revcomp;
use molecular_annotation::MolecularAnnotations;
use rayon::iter::ParallelIterator;
use rayon::prelude::IntoParallelRefMutIterator;
use rust_htslib::bam::Read;
use rust_htslib::bam::Record;

pub fn strip_base_mods(opts: &mut StripBasemodsOptions) {
    let mut bam = opts.input.bam_reader();
    let mut out = opts.input.bam_writer(&opts.out);
    // read in bam data
    let bam_chunk_iter = BamChunk::new(bam.records(), None);

    let filter_mod = match &opts.basemod {
        Some(mod_name) => mod_name,
        None => "",
    };

    // iterate over chunks
    for mut chunk in bam_chunk_iter {
        let records: Vec<&mut Record> = chunk
            .par_iter_mut()
            .map(|record| {
                let mut annot = ma_io::read_record(record).unwrap_or_else(|e| {
                    log::warn!("read_record failed for {:?}: {e}", String::from_utf8_lossy(record.qname()));
                    MolecularAnnotations::from_record(record)
                });

                // Drop whole annotation types.
                if filter_mod == "5mC" || filter_mod == "CpG" {
                    annot.annotation_types.retain(|t| t.name != CPG_TYPE);
                } else if filter_mod == "6mA" || filter_mod == "m6A" {
                    annot.annotation_types.retain(|t| t.name != M6A_TYPE);
                }

                // Drop by canonical-base / strand. The "strand" the legacy
                // BaseMods stored was the MM tag's strand sign (`+` or `-`).
                // After S2a the per-call strand is determined at write
                // time from the forward seq: m6a calls on A → `+`, on T → `-`;
                // 5mC calls on C → `+`, on G → `+` (canonical fiberseq
                // never emits `-` for 5mC). drop_forward removes calls
                // whose write-time strand would be `+`; drop_reverse `-`.
                if opts.drop_forward || opts.drop_reverse {
                    let forward_seq = if record.is_reverse() {
                        revcomp(record.seq().as_bytes())
                    } else {
                        record.seq().as_bytes()
                    };
                    if opts.drop_forward {
                        annot.retain(M6A_TYPE, |a| {
                            forward_seq.get(a.start as usize).copied() != Some(b'A')
                        });
                        annot.retain(CPG_TYPE, |a| {
                            forward_seq.get(a.start as usize).copied() != Some(b'C')
                        });
                    }
                    if opts.drop_reverse {
                        annot.retain(M6A_TYPE, |a| {
                            forward_seq.get(a.start as usize).copied() != Some(b'T')
                        });
                        // canonical fiberseq doesn't emit `-` for 5mC; nothing to drop.
                    }
                }

                if opts.ml_m6a > 0 {
                    annot.retain(M6A_TYPE, |a| {
                        primary_qual(&a.qualities, M6A_TYPE) >= opts.ml_m6a
                    });
                }
                if opts.ml_5mc > 0 {
                    annot.retain(CPG_TYPE, |a| {
                        primary_qual(&a.qualities, CPG_TYPE) >= opts.ml_5mc
                    });
                }

                ma_io::ensure_basemod_encoding(&mut annot);
                ma_io::write_record(record, &annot);
                record
            })
            .collect();

        records
            .into_iter()
            .for_each(|record| out.write(record).unwrap());
    }
}
