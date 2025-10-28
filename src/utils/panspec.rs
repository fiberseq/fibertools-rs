use crate::utils::bio_io::header_from_hashmap;
use rust_htslib::bam::{self, Header};

pub fn strip_pan_spec_header(header: &bam::Header, pan_spec_delimiter: &char) -> Header {
    let mut hash_map: std::collections::HashMap<
        String,
        Vec<linear_map::LinearMap<String, String>>,
    > = header.to_hashmap();
    for (key, value) in hash_map.iter_mut() {
        if key.eq("SQ") {
            for sn_line in value.iter_mut() {
                let name = sn_line
                    .get_mut("SN")
                    .expect("SN tag not found within an @SQ line");
                let mut del_count = 0;
                let mut new_name = String::new();
                for char in name.chars() {
                    if del_count >= 2 {
                        new_name.push(char);
                    } else if char == *pan_spec_delimiter {
                        del_count += 1;
                    }
                }
                name.clear();
                name.push_str(&new_name);
            }
        }
    }
    let mut new_header = header_from_hashmap(hash_map);

    // Preserve comments from original header
    for comment in header.comments() {
        new_header.push_comment(comment.as_bytes());
    }

    new_header
}

pub fn add_pan_spec_header(header: &bam::Header, pan_spec_prefix: &str) -> Header {
    let mut hash_map = header.to_hashmap();
    for (key, value) in hash_map.iter_mut() {
        if key.eq("SQ") {
            for sn_line in value.iter_mut() {
                let name = sn_line
                    .get_mut("SN")
                    .expect("SN tag not found within an @SQ line");
                let mut new_name = String::new();
                new_name.push_str(pan_spec_prefix);
                new_name.push_str(name);
                name.clear();
                name.push_str(&new_name);
            }
        }
    }
    let mut new_header = header_from_hashmap(hash_map);

    // Preserve comments from original header
    for comment in header.comments() {
        new_header.push_comment(comment.as_bytes());
    }

    new_header
}
