use crate::utils::bamranges;

#[derive(Debug)]
enum Threshold {
    Single(f64),
    Range(f64, f64),
}

fn cmp(name: &str, value: f64, operator: &str, threshold: &Threshold) -> bool {
    match threshold {
        Threshold::Single(thr) => match operator {
            ">" => value > *thr,
            "<" => value < *thr,
            ">=" => value >= *thr,
            "<=" => value <= *thr,
            "=" => value == *thr,
            "!=" => value != *thr,
            _ => {
                eprintln!("Invalid operator for {} function: {}", name, operator);
                std::process::exit(1);
            }
        },
        Threshold::Range(min, max) => match operator {
            "=" => (value as f64) >= *min && (value as f64) < *max,
            _ => {
                eprintln!("Invalid operator for {} function with range: {}", name, operator);
                std::process::exit(1);
            }
        },
    }
}

fn len(value: i64, operator: &str, threshold: &Threshold) -> bool {
    cmp("len", value as f64, operator, threshold)
}

fn qual(value: u8, operator: &str, threshold: &Threshold) -> bool {
    cmp("qual", value as f64, operator, threshold)
}

pub struct ParsedExpr {
    fn_name: String,
    pub rng_name: String,
    op: String, // <, <=, >, >=, etc
    threshold: Threshold,
}

pub fn parse_filter(filter_orig: &str) -> ParsedExpr {
    let mut filter = filter_orig.to_string();
    filter.retain(|c| !c.is_whitespace());

    let func_name_end = filter.find('(').unwrap_or(filter.len());
    let func_name = filter[..func_name_end].trim().to_string();

    let gnm_feat_start = filter.find('(').unwrap_or(filter.len()) + 1;
    let gnm_feat_end = filter.find(')').unwrap_or(filter.len());
    let gnm_feat = filter[gnm_feat_start..gnm_feat_end].to_string();
    if !["msp", "nuc", "m6a", "5mC"].contains(&gnm_feat.as_str()) {
        eprintln!("Invalid argument for len function: {}", gnm_feat);
        std::process::exit(1);
    }

    let rest = &filter[gnm_feat_end + 1..].trim();

    let operators = ["!=", ">=", "<=", ">", "<", "="];
    let mut operator = "".to_string();
    let mut threshold = None;
    let mut range = None;

    for &op in operators.iter() {
        if let Some(pos) = rest.find(op) {
            operator = op.to_string();
            let threshold_str = rest[pos + op.len()..].trim();
            if threshold_str.contains(':') {
                let range_parts: Vec<&str> = threshold_str.split(':').collect();
                if range_parts.len() == 2 {
                    range = Some((
                        range_parts[0].trim().parse::<f64>().unwrap(),
                        range_parts[1].trim().parse::<f64>().unwrap(),
                    ));
                }
            } else {
                threshold = Some(threshold_str.parse::<f64>().unwrap());
            }
            break;
        }
    }

    if let Some((_, _)) = range {
        if operator != "=" {
            eprintln!("Range thresholds can only be used with the '=' operator.");
            std::process::exit(1);
        }
    }

    let threshold_value = match range {
        Some((min, max)) => Threshold::Range(min, max),
        None => Threshold::Single(threshold.unwrap()),
    };

    ParsedExpr {
        fn_name: func_name,
        rng_name: gnm_feat,
        op: operator,
        threshold: threshold_value,
    }
}

pub fn apply_filter_to_range(
    parsed: &ParsedExpr,
    mut range: bamranges::Ranges,
) -> Result<bamranges::Ranges, anyhow::Error> {
    let starting_len = range.starts.len();

    let to_keep: Vec<bool> = if parsed.fn_name == "len" {
        range
            .lengths
            .iter()
            .map(|l| len(l.unwrap(), &parsed.op, &parsed.threshold))
            .collect()
    } else if parsed.fn_name == "qual" {
        range
            .qual
            .iter()
            .map(|q| qual(*q, &parsed.op, &parsed.threshold))
            .collect()
    } else {
        anyhow::bail!("Invalid function name: {}", &parsed.fn_name);
    };

    // drop i64 values from the range
    for vec in vec![
        &mut range.starts,
        &mut range.ends,
        &mut range.lengths,
        &mut range.reference_starts,
        &mut range.reference_ends,
        &mut range.reference_lengths,
    ] {
        *vec = vec
            .iter()
            .zip(to_keep.iter())
            .filter_map(|(v, d)| if *d { Some(*v) } else { None })
            .collect();
    }

    // drop u8 values from the range
    range.qual = range
        .qual
        .iter()
        .zip(to_keep.iter())
        .filter_map(|(v, d)| if *d { Some(*v) } else { None })
        .collect();

    // check we dropped the right number of values
    let n_dropped = to_keep.iter().filter(|&&x| !x).count();
    assert_eq!(starting_len, range.starts.len() + n_dropped);

    Ok(range)
}

/// tests
#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::bamranges;

    fn make_fake_range() -> bamranges::Ranges {
        bamranges::Ranges {
            starts: vec![Some(0), Some(10), Some(17)],
            ends: vec![Some(5), Some(15), Some(20)],
            lengths: vec![Some(5), Some(5), Some(3)],
            qual: vec![0, 255, 181],
            reference_starts: vec![Some(0), Some(10), Some(17)],
            reference_ends: vec![Some(5), Some(15), Some(20)],
            reference_lengths: vec![Some(5), Some(5), Some(3)],
            seq_len: 100,
            reverse: false,
        }
    }

    #[test]
    fn test_this_one() {
        let filter = "len(msp)=50:100";
        let mut range = make_fake_range();
        let parser = parse_filter(&filter);
        eprintln!("{:?}", range.starts.len());
        apply_filter_to_range(&parser, &mut range).unwrap();
        eprintln!("{:?}", range.starts.len());
    }
}
