use super::cli::DecoratorOptions;
use super::fiber::FiberseqData;
use super::*;
use crate::fiber::FiberseqRecords;
use anyhow;
use rust_htslib::bam::ext::BamRecordExtensions;

const NUC_COLOR: &str = "169,169,169";
const M6A_COLOR: &str = "128,0,128";
const CPG_COLOR: &str = "139,69,19";
const LINKER_COLOR: &str = "147,112,219";
const FIRE_COLORS: [(f32, &str); 8] = [
    (1.0, "139,0,0"),
    (2.0, "175,0,0"),
    (3.0, "200,0,0"),
    (4.0, "225,0,0"),
    (5.0, "255,0,0"),
    (10.0, "255,140,0"),
    (25.0, "225,225,0"),
    (100.0, LINKER_COLOR),
];

pub fn pos_to_string(pos: &[Option<i64>]) -> String {
    pos.iter().flatten().map(|p| p.to_string() + ",").collect()
}

fn decorator_from_positions(
    fiber: &FiberseqData,
    starts: &[Option<i64>],
    lengths: Option<&[Option<i64>]>,
    color: &str,
    element_type: &str,
) -> String {
    // clean the starts and lengths
    let starts: Vec<i64> = starts.iter().flatten().copied().collect();
    let lengths: Vec<i64> = if let Some(l) = lengths {
        l.iter().flatten().copied().collect()
    } else {
        vec![1; starts.len()]
    };
    // skip if not ref positions
    if starts.is_empty() {
        return "".to_string();
    }
    // get fields
    let strand = if fiber.record.is_reverse() { '-' } else { '+' };
    let start = starts[0];
    let end = starts[starts.len() - 1] + lengths[lengths.len() - 1];
    let tag = format!(
        "{}:{}-{}:{}",
        fiber.target_name,
        fiber.record.reference_start(),
        fiber.record.reference_end(),
        String::from_utf8_lossy(fiber.record.qname())
    );
    let bed6 = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t",
        fiber.target_name, start, end, element_type, 0, strand
    );
    let bed12 = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t",
        start,
        end,
        color.to_owned() + ",1",
        lengths.len(),
        join_by_str(lengths, ","),
        join_by_str(starts.iter().map(|p| p - start), ","),
    );
    let decorator = format!(
        "{}\t{}\t{}\t{}\t{}\n",
        tag,
        "block",
        color.to_owned() + ",0",
        "Ignored",
        element_type,
    );
    bed6 + &bed12 + &decorator
}

pub fn bed12_from_fiber(fiber: &FiberseqData) -> String {
    let strand = if fiber.record.is_reverse() { '-' } else { '+' };
    let bed6 = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t",
        fiber.target_name,
        fiber.record.reference_start(),
        fiber.record.reference_end(),
        String::from_utf8_lossy(fiber.record.qname()),
        0,
        strand
    );
    let length = fiber.record.reference_end() - fiber.record.reference_start() - 1;
    let bed12 = format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        fiber.record.reference_start(),
        fiber.record.reference_end(),
        "0,0,0",
        2,
        "1,1",
        "0,".to_owned() + &length.to_string(),
        fiber.get_hp(),
    );
    bed6 + &bed12
}

pub fn fire_decorators(fiber: &FiberseqData) -> String {
    // hashmap with colors and empty vecs
    let mut map = std::collections::HashMap::new();
    for (_, color) in FIRE_COLORS.iter() {
        map.insert(*color, vec![]);
    }

    for ((pos, length), qual) in fiber
        .msp
        .reference_starts
        .iter()
        .zip(fiber.msp.reference_lengths.iter())
        .zip(fiber.msp.qual.iter())
    {
        if let (Some(p), Some(l)) = (pos, length) {
            let l = Some(*l);
            let p = Some(*p);
            let fdr_val = 100.0 - *qual as f32 / 255.0 * 100.0;
            for (fdr, color) in FIRE_COLORS.iter() {
                if fdr_val <= *fdr {
                    map.get_mut(*color).unwrap().push((p, l));
                    break;
                }
            }
        }
    }
    let mut rtn = "".to_string();
    for (color, values) in map.into_iter() {
        let (starts, lengths): (Vec<Option<i64>>, Vec<Option<i64>>) = values.into_iter().unzip();
        rtn += &decorator_from_positions(fiber, &starts, Some(&lengths), color, "FIRE");
    }
    rtn
}

pub fn decorator_from_bam(fiber: &FiberseqData) -> (String, String) {
    let decorators =
        decorator_from_positions(fiber, &fiber.m6a.reference_starts, None, M6A_COLOR, "m6A")
            + &decorator_from_positions(fiber, &fiber.cpg.reference_starts, None, CPG_COLOR, "5mC")
            + &decorator_from_positions(
                fiber,
                &fiber.nuc.reference_starts,
                Some(&fiber.nuc.reference_lengths),
                NUC_COLOR,
                "NUC",
            )
            + &fire_decorators(fiber);
    (bed12_from_fiber(fiber), decorators)
}

pub fn get_decorators_from_bam(dec_opts: &DecoratorOptions) -> Result<(), anyhow::Error> {
    let mut bam = bio_io::bam_reader(&dec_opts.bam, 8);
    let mut bed12 = bio_io::writer(&dec_opts.bed12)?;
    let mut decorator = bio_io::writer(&dec_opts.decorator)?;
    for rec in FiberseqRecords::new(&mut bam, dec_opts.min_ml_score) {
        if rec.record.is_secondary() || rec.record.is_supplementary() || rec.record.is_unmapped() {
            continue;
        }
        let (fiber_bed12, fiber_decorator) = decorator_from_bam(&rec);
        bed12.write_all(fiber_bed12.as_bytes())?;
        decorator.write_all(fiber_decorator.as_bytes())?;
    }
    Ok(())
}
