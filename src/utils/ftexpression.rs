use crate::fiber::FiberseqData;
use crate::utils::basemods::{CPG_TYPE, M6A_TYPE};
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

pub fn apply_filter_fsd(fsd: &mut FiberseqData, filt: &FiberFilters) -> Result<(), anyhow::Error> {
    if let Some(s) = filt.filter_expression.as_ref() {
        if !s.is_empty() {
            for parser in parse_filter_all(s.as_str()) {
                let type_name = match parser.feat_name.as_str() {
                    "msp" => "msp",
                    "nuc" => "nuc",
                    "m6a" => M6A_TYPE,
                    "5mC" => CPG_TYPE,
                    other => anyhow::bail!("Unknown feature name: {}", other),
                };
                match parser.fn_name.as_str() {
                    "len" => fsd.annotations.retain(type_name, |a| {
                        len(a.length as i64, &parser.op, &parser.threshold)
                    }),
                    "qual" if type_name == "msp" => {
                        // qual(msp) historically meant "filter MSPs by FIRE
                        // precision" (legacy fibertools wrote FIRE precisions
                        // onto the MSP `aq` tag). Post-MA, that quality lives
                        // on the `fire` annotation type instead. Collect
                        // per-MSP precisions in *molecular* order to align
                        // with retain's iteration, then drop msp and fire
                        // in lockstep.
                        let primary = crate::utils::bamannotations::primary_qual;
                        let mol_quals: Vec<u8> = if let Some(f) = fsd
                            .annotations
                            .get_type("fire")
                            .filter(|f| {
                                fsd.annotations
                                    .get_type("msp")
                                    .is_some_and(|m| m.annotations.len() == f.annotations.len())
                            }) {
                            f.annotations
                                .iter()
                                .map(|a| primary(&a.qualities, "fire"))
                                .collect()
                        } else if let Some(m) = fsd.annotations.get_type("msp") {
                            m.annotations
                                .iter()
                                .map(|a| primary(&a.qualities, "msp"))
                                .collect()
                        } else {
                            Vec::new()
                        };
                        let keep: Vec<bool> = mol_quals
                            .iter()
                            .map(|q| qual(*q, &parser.op, &parser.threshold))
                            .collect();
                        let has_fire = fsd.annotations.get_type("fire").is_some();
                        let mut i = 0;
                        fsd.annotations.retain("msp", |_| {
                            let k = keep.get(i).copied().unwrap_or(false);
                            i += 1;
                            k
                        });
                        if has_fire {
                            let mut i = 0;
                            fsd.annotations.retain("fire", |_| {
                                let k = keep.get(i).copied().unwrap_or(false);
                                i += 1;
                                k
                            });
                        }
                    }
                    "qual" => fsd.annotations.retain(type_name, |a| {
                        qual(
                            crate::utils::bamannotations::primary_qual(&a.qualities, type_name),
                            &parser.op,
                            &parser.threshold,
                        )
                    }),
                    other => anyhow::bail!("Invalid function name: {}", other),
                }
            }
        }
    }
    Ok(())
}
