use super::cli::DecoratorOptions;
use super::fiber::FiberseqData;
use super::*;
use anyhow;
use rust_htslib::bam::ext::BamRecordExtensions;
use std::cmp::Ordering;
use std::fmt::Display;

pub fn get_fire_color(fdr: f32) -> &'static str {
    for (fdr_val, color) in FIRE_COLORS.iter() {
        if fdr <= *fdr_val {
            return color;
        }
    }
    LINKER_COLOR
}

pub fn pos_to_string(pos: &[Option<i64>]) -> String {
    pos.iter().flatten().map(|p| p.to_string() + ",").collect()
}

#[derive(Debug, Clone)]
pub struct Decorator<'a> {
    pub starts: Vec<i64>,
    pub lengths: Vec<i64>,
    pub fiber: &'a FiberseqData,
    pub color: &'static str,
    pub element_type: &'static str,
    pub start: i64,
    pub end: i64,
}

impl<'a> Decorator<'a> {
    pub fn new(
        fiber: &'a FiberseqData,
        starts: &[Option<i64>],
        lengths: Option<&[Option<i64>]>,
        color: &'static str,
        element_type: &'static str,
    ) -> Self {
        // clean the starts and lengths
        let starts: Vec<i64> = starts.iter().flatten().copied().collect();
        let lengths: Vec<i64> = if let Some(l) = lengths {
            l.iter().flatten().copied().collect()
        } else {
            vec![1; starts.len()]
        };

        let start = if starts.is_empty() {
            fiber.record.reference_start()
        } else {
            starts[0]
        };
        let end = if starts.is_empty() {
            fiber.record.reference_end()
        } else {
            starts[starts.len() - 1] + lengths[lengths.len() - 1]
        };

        Self {
            starts,
            lengths,
            fiber,
            color,
            element_type,
            start,
            end,
        }
    }
}

impl Display for Decorator<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // skip if not ref positions
        if self.starts.is_empty() {
            return write!(f, "");
        }
        // get fields
        let strand = if self.fiber.record.is_reverse() {
            '-'
        } else {
            '+'
        };
        let tag = format!(
            "{}:{}-{}:{}",
            self.fiber.target_name,
            self.fiber.record.reference_start(),
            self.fiber.record.reference_end(),
            String::from_utf8_lossy(self.fiber.record.qname())
        );
        let bed6 = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            self.fiber.target_name, self.start, self.end, self.element_type, 0, strand
        );
        let bed12 = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t",
            self.start,
            self.end,
            self.color.to_owned() + ",1",
            self.lengths.len(),
            join_by_str(&self.lengths, ","),
            join_by_str(self.starts.iter().map(|p| p - self.start), ","),
        );
        let decorator = format!(
            "{}\t{}\t{}\t{}\t{}\n",
            tag,
            "block",
            self.color.to_owned() + ",0",
            "Ignored",
            self.element_type,
        );
        write!(f, "{}", bed6 + &bed12 + &decorator)
    }
}

impl<'a> Ord for Decorator<'a> {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.start, self.end).cmp(&(other.start, other.end))
    }
}

impl<'a> PartialOrd for Decorator<'a> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<'a> PartialEq for Decorator<'a> {
    fn eq(&self, other: &Self) -> bool {
        (self.start, self.end) == (other.start, other.end)
    }
}

impl<'a> Eq for Decorator<'a> {}

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

pub fn fire_decorators(fiber: &FiberseqData) -> Vec<Decorator> {
    // hashmap with colors and empty vecs
    let mut map = std::collections::HashMap::new();
    for (_, color) in FIRE_COLORS.iter() {
        map.insert(color, vec![]);
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
            let fire_color = get_fire_color(fdr_val);
            map.get_mut(&fire_color).unwrap().push((p, l));
        }
    }
    let mut rtn = vec![];
    for (color, values) in map.into_iter() {
        let (starts, lengths): (Vec<Option<i64>>, Vec<Option<i64>>) = values.into_iter().unzip();
        rtn.push(Decorator::new(
            fiber,
            &starts,
            Some(&lengths),
            color,
            "FIRE",
        ));
    }
    rtn
}

pub fn decorator_from_bam(fiber: &FiberseqData) -> (String, Vec<Decorator>) {
    let mut decorators = vec![
        Decorator::new(fiber, &fiber.m6a.reference_starts, None, M6A_COLOR, "m6A"),
        Decorator::new(fiber, &fiber.cpg.reference_starts, None, CPG_COLOR, "5mC"),
        Decorator::new(
            fiber,
            &fiber.nuc.reference_starts,
            Some(&fiber.nuc.reference_lengths),
            NUC_COLOR,
            "NUC",
        ),
    ];
    decorators.append(&mut fire_decorators(fiber));
    (bed12_from_fiber(fiber), decorators)
}

pub fn get_decorators_from_bam(dec_opts: &mut DecoratorOptions) -> Result<(), anyhow::Error> {
    let mut bam = dec_opts.input.bam_reader();
    let mut bed12 = bio_io::writer(&dec_opts.bed12)?;
    let mut decorator = bio_io::writer(&dec_opts.decorator)?;

    for rec in dec_opts.input.fibers(&mut bam) {
        let (fiber_bed12, fiber_decorator) = decorator_from_bam(&rec);
        bed12.write_all(fiber_bed12.as_bytes())?;
        for fiber_decorator in fiber_decorator {
            decorator.write_all(fiber_decorator.to_string().as_bytes())?;
        }
    }
    Ok(())
}
