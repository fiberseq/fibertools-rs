use std::env;

#[derive(Debug)]
enum Threshold {
    Single(f64),
    Range(f64, f64),
}

fn len(name: &str, value: i32, operator: &str, threshold: &Threshold) -> bool {
    if !["msp", "fire", "nuc"].contains(&name) {
        eprintln!("Invalid argument for len function: {}", name);
        std::process::exit(1);
    }

    match threshold {
        Threshold::Single(thr) => match operator {
            ">" => (value as f64) > *thr,
            "<" => (value as f64) < *thr,
            ">=" => (value as f64) >= *thr,
            "<=" => (value as f64) <= *thr,
            "=" => (value as f64) == *thr,
            "!=" => (value as f64) != *thr,
            _ => {
                eprintln!("Invalid operator for len function: {}", operator);
                std::process::exit(1);
            }
        },
        Threshold::Range(min, max) => match operator {
            "=" => (value as f64) >= *min && (value as f64) < *max,
            _ => {
                eprintln!("Invalid operator for len function with range: {}", operator);
                std::process::exit(1);
            }
        },
    }
}

fn qual(name: &str, value: f64, operator: &str, threshold: &Threshold) -> bool {
    if !["m6a", "5mC"].contains(&name) {
        eprintln!("Invalid argument for qual function: {}", name);
        std::process::exit(1);
    }

    match threshold {
        Threshold::Single(thr) => match operator {
            ">" => value > *thr,
            "<" => value < *thr,
            ">=" => value >= *thr,
            "<=" => value <= *thr,
            "=" => value == *thr,
            "!=" => value != *thr,
            _ => {
                eprintln!("Invalid operator for qual function: {}", operator);
                std::process::exit(1);
            }
        },
        Threshold::Range(min, max) => match operator {
            "=" => value >= *min && value < *max,
            _ => {
                eprintln!("Invalid operator for qual function with range: {}", operator);
                std::process::exit(1);
            }
        },
    }
}

fn parse_filter(filter: &str) -> (String, String, String, Threshold) {
    let func_name_end = filter.find('(').unwrap_or(filter.len());
    let func_name = filter[..func_name_end].trim().to_string();

    let gnm_feat_start = filter.find('(').unwrap_or(filter.len()) + 1;
    let gnm_feat_end = filter.find(')').unwrap_or(filter.len());
    let gnm_feat = filter[gnm_feat_start..gnm_feat_end].trim().to_string();

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

    (func_name, gnm_feat, operator, threshold_value)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Usage: cargo run -- <filter>");
        std::process::exit(1);
    }

    let filter = &args[1];
    let (func_name, gnm_feat, operator, threshold) = parse_filter(filter);

    let my_gnm_feats = vec![30, 40, 50, 60, 70, 80, 90, 100];
    let my_m6as = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];

    if func_name == "len" {
        let filtered_gnm_feats: Vec<i32> = my_gnm_feats
            .into_iter()
            .filter(|&msp| len(&gnm_feat, msp, &operator, &threshold))
            .collect();
        println!("Filtered GNM_FEATS: {:?}", filtered_gnm_feats);
    } else if func_name == "qual" {
        let filtered_m6as: Vec<f64> = my_m6as
            .into_iter()
            .filter(|&qual_value| qual(&gnm_feat, qual_value, &operator, &threshold))
            .collect();
        println!("Filtered M6AS: {:?}", filtered_m6as);
    } else {
        eprintln!("Unsupported function: {}", func_name);
        std::process::exit(1);
    }
}
