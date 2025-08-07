use crate::fiber::FiberseqData;
use crate::utils::bamannotations;
use crate::utils::input_bam::FiberFilters;

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
                eprintln!("Invalid operator for {name} function: {operator}");
                std::process::exit(1);
            }
        },
        Threshold::Range(min, max) => match operator {
            "=" => value >= *min && value < *max,
            _ => {
                eprintln!("Require = operator for {name} when using ranged args");
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
    feat_name: String,
    op: String, // <, <=, >, >=, etc
    threshold: Threshold,
}

pub fn parse_filter_all(full_expr: &str) -> Vec<ParsedExpr> {
    let mut v = Vec::new();
    let expressions: Vec<String> = full_expr.split(',').map(|s| s.to_string()).collect();
    for expr in expressions.iter() {
        v.push(parse_filter(expr.as_str()));
    }
    v
}

pub fn parse_filter(filter_orig: &str) -> ParsedExpr {
    // Remove all whitespace characters
    let mut filter = filter_orig.to_string();
    filter.retain(|c| !c.is_whitespace());

    // Extract and validate function name
    let func_name_end = filter.find('(').unwrap_or(filter.len());
    let func_name = filter[..func_name_end].trim().to_string();
    if !["len", "qual"].contains(&func_name.as_str()) {
        eprintln!("Invalid function: {func_name}");
        std::process::exit(1);
    }

    // Extract and validate function argument
    let gnm_feat_start = func_name_end + 1;
    let gnm_feat_end = filter.find(')').unwrap_or(filter.len());
    let gnm_feat = filter[gnm_feat_start..gnm_feat_end].to_string();
    if !["msp", "nuc", "m6a", "5mC"].contains(&gnm_feat.as_str()) {
        eprintln!("Invalid argument for {func_name} function: {gnm_feat}");
        std::process::exit(1);
    }

    // Extract and validate the operator and value
    let rest = &filter[gnm_feat_end + 1..].trim();
    let operators = ["!=", ">=", "<=", ">", "<", "="];
    let mut operator = None;
    let mut value_str = "";

    for &op in operators.iter() {
        if let Some(pos) = rest.find(op) {
            operator = Some(op);
            value_str = rest[pos + op.len()..].trim();
            break;
        }
    }

    if operator.is_none() {
        eprintln!("No valid operator found.");
        std::process::exit(1);
    }

    let operator = operator.unwrap().to_string();

    // Validate and parse the value
    let value = if value_str.contains(':') {
        let range_parts: Vec<&str> = value_str.split(':').collect();
        if range_parts.len() != 2 {
            eprintln!("Range format requires exactly 2 values separated by ':', ex: '50:100'");
            std::process::exit(1);
        }
        let start = range_parts[0].trim();
        let end = range_parts[1].trim();
        let start_value = start.parse::<f64>();
        let end_value = end.parse::<f64>();
        if start_value.is_err() || end_value.is_err() {
            eprintln!("Range thresholds must be numeric values, ex: '30:100'");
            std::process::exit(1);
        }
        let start_value = start_value.unwrap();
        let end_value = end_value.unwrap();
        if start_value >= end_value {
            eprintln!("In a range, the start value must be less than the end value.");
            std::process::exit(1);
        }
        Threshold::Range(start_value, end_value)
    } else {
        let value_num = value_str.parse::<f64>();
        if value_num.is_err() {
            eprintln!("Value is not numeric: {value_str}.");
            std::process::exit(1);
        }
        Threshold::Single(value_num.unwrap())
    };

    ParsedExpr {
        fn_name: func_name,
        feat_name: gnm_feat,
        op: operator,
        threshold: value,
    }
}

pub fn apply_filter_to_range(
    parsed: &ParsedExpr,
    range: &mut bamannotations::Ranges,
) -> Result<(), anyhow::Error> {
    let starting_len = range.annotations.len();

    let to_keep: Vec<bool> = if parsed.fn_name == "len" {
        range
            .annotations
            .iter()
            .map(|annotation| len(annotation.length, &parsed.op, &parsed.threshold))
            .collect()
    } else if parsed.fn_name == "qual" {
        range
            .annotations
            .iter()
            .map(|annotation| qual(annotation.qual, &parsed.op, &parsed.threshold))
            .collect()
    } else {
        anyhow::bail!("Invalid function name: {}", &parsed.fn_name);
    };

    // Filter annotations based on the to_keep boolean vector
    range.annotations = range
        .annotations
        .iter()
        .zip(to_keep.iter())
        .filter_map(|(annotation, &keep)| if keep { Some(annotation.clone()) } else { None })
        .collect();

    // check we dropped the right number of values
    let n_dropped = to_keep.iter().filter(|&x| !x).count();
    assert_eq!(starting_len, range.annotations.len() + n_dropped);

    Ok(())
}

pub fn apply_filter_fsd(fsd: &mut FiberseqData, filt: &FiberFilters) -> Result<(), anyhow::Error> {
    if let Some(s) = filt.filter_expression.as_ref() {
        if !s.is_empty() {
            let parsers = parse_filter_all(s.as_str());
            for parser in parsers.iter() {
                match parser.feat_name.as_str() {
                    "msp" => apply_filter_to_range(parser, &mut fsd.msp)?,
                    "nuc" => apply_filter_to_range(parser, &mut fsd.nuc)?,
                    "m6a" => apply_filter_to_range(parser, &mut fsd.m6a)?,
                    "5mC" => apply_filter_to_range(parser, &mut fsd.cpg)?,
                    _ => {
                        return Err(anyhow::anyhow!(
                            "Unknown feature name: {}",
                            parser.feat_name
                        ));
                    }
                }
            }
        }
    }
    Ok(())
}

/// tests
#[cfg(test)]
mod test {
    use super::*;
    use crate::utils::bamannotations;

    fn make_fake_range() -> bamannotations::Ranges {
        use bamannotations::{FiberAnnotation, FiberAnnotations};

        let annotations = vec![
            FiberAnnotation {
                start: 0,
                end: 5,
                length: 5,
                qual: 0,
                reference_start: Some(0),
                reference_end: Some(5),
                reference_length: Some(5),
                extra_columns: None,
            },
            FiberAnnotation {
                start: 10,
                end: 15,
                length: 5,
                qual: 255,
                reference_start: Some(10),
                reference_end: Some(15),
                reference_length: Some(5),
                extra_columns: None,
            },
            FiberAnnotation {
                start: 17,
                end: 20,
                length: 3,
                qual: 181,
                reference_start: Some(17),
                reference_end: Some(20),
                reference_length: Some(3),
                extra_columns: None,
            },
        ];

        FiberAnnotations {
            annotations,
            seq_len: 100,
            reverse: false,
        }
    }

    #[test]
    fn test_this_one() {
        let filter = "len(msp)=50:100";
        let mut range = make_fake_range();
        let parser = parse_filter(filter);
        eprintln!("{:?}", range.annotations.len());
        apply_filter_to_range(&parser, &mut range).unwrap();
        eprintln!("{:?}", range.annotations.len());
    }
}
