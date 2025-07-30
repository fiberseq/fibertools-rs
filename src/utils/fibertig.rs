use crate::utils::bio_io::get_u32_tag;
use rust_htslib::bam::{record::Aux, Record};
use std::collections::HashMap;

pub const FIBERTIG_DELIMITER: &str = "|\t|";

#[derive(Debug, Clone, PartialEq)]
pub struct BedInterval {
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub extra_columns: Vec<String>,
}

impl BedInterval {
    pub fn from_bed_line(line: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() < 3 {
            return Err("BED line must have at least 3 columns (chrom, start, end)".into());
        }

        let chrom = fields[0].to_string();
        let start: i64 = fields[1].parse()?;
        let end: i64 = fields[2].parse()?;
        
        let extra_columns = if fields.len() > 3 {
            fields[3..].iter().map(|s| s.to_string()).collect()
        } else {
            vec![]
        };

        Ok(BedInterval {
            chrom,
            start,
            end,
            extra_columns,
        })
    }
}

#[derive(Debug, Clone)]
pub struct FiberTigAnnotations {
    pub intervals: Vec<BedInterval>,
}

impl FiberTigAnnotations {
    pub fn new() -> Self {
        Self {
            intervals: Vec::new(),
        }
    }

    pub fn from_bed_file(bed_path: &str) -> Result<HashMap<String, Self>, Box<dyn std::error::Error>> {
        use std::fs::File;
        use std::io::{BufRead, BufReader};

        let file = File::open(bed_path)?;
        let reader = BufReader::new(file);
        let mut annotations_by_chrom: HashMap<String, Self> = HashMap::new();

        for line in reader.lines() {
            let line = line?;
            if line.starts_with('#') || line.trim().is_empty() {
                continue;
            }

            let interval = BedInterval::from_bed_line(&line)?;
            let chrom = interval.chrom.clone();
            
            annotations_by_chrom
                .entry(chrom)
                .or_insert_with(Self::new)
                .intervals
                .push(interval);
        }

        // Sort intervals by start position within each chromosome
        for annotations in annotations_by_chrom.values_mut() {
            annotations.intervals.sort_by_key(|i| i.start);
        }

        Ok(annotations_by_chrom)
    }
}

pub fn set_fibertig_tags(record: &mut Record, intervals: &[BedInterval]) -> Result<(), Box<dyn std::error::Error>> {
    if intervals.is_empty() {
        return Ok(());
    }

    // Set start/end position arrays
    let starts: Vec<u32> = intervals.iter().map(|b| b.start as u32).collect();
    let ends: Vec<u32> = intervals.iter().map(|b| b.end as u32).collect();
    
    let start_aux = Aux::ArrayU32(starts);
    let end_aux = Aux::ArrayU32(ends);
    record.push_aux(b"fs", start_aux)?;
    record.push_aux(b"fe", end_aux)?;
    
    // Set raw BED data string (columns 4+)
    let bed_data: String = intervals
        .iter()
        .map(|b| b.extra_columns.join("\t"))
        .collect::<Vec<_>>()
        .join(FIBERTIG_DELIMITER);
    
    record.push_aux(b"fd", Aux::String(&bed_data))?;
    
    // Mark as FiberTig
    record.push_aux(b"ft", Aux::String("1.0"))?;
    
    Ok(())
}

pub fn get_fibertig_starts(record: &Record) -> Vec<i64> {
    get_u32_tag(record, b"fs")
}

pub fn get_fibertig_ends(record: &Record) -> Vec<i64> {
    get_u32_tag(record, b"fe")
}

pub fn get_fibertig_data(record: &Record) -> Vec<Vec<String>> {
    if let Ok(Aux::String(data)) = record.aux(b"fd") {
        data.split(FIBERTIG_DELIMITER)
            .map(|interval_data| {
                if interval_data.is_empty() {
                    vec![]
                } else {
                    interval_data.split('\t').map(|s| s.to_string()).collect()
                }
            })
            .collect()
    } else {
        vec![]
    }
}

pub fn is_fibertig_record(record: &Record) -> bool {
    matches!(record.aux(b"ft"), Ok(Aux::String(_)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bed_interval_parsing() {
        let line = "chr1\t100\t200\tgene1\t1000\t+";
        let interval = BedInterval::from_bed_line(line).unwrap();
        
        assert_eq!(interval.chrom, "chr1");
        assert_eq!(interval.start, 100);
        assert_eq!(interval.end, 200);
        assert_eq!(interval.extra_columns, vec!["gene1", "1000", "+"]);
    }

    #[test]
    fn test_minimal_bed() {
        let line = "chr1\t100\t200";
        let interval = BedInterval::from_bed_line(line).unwrap();
        
        assert_eq!(interval.chrom, "chr1");
        assert_eq!(interval.start, 100);
        assert_eq!(interval.end, 200);
        assert!(interval.extra_columns.is_empty());
    }
}