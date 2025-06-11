use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use fibertools_rs::utils::bamranges::Ranges;
use rust_htslib::bam;
use std::time::Duration;

/// Create a mock BAM record for benchmarking using the test BAM data
fn get_test_bam_record() -> bam::Record {
    // Use the existing test data infrastructure
    use fibertools_rs::utils::bio_io;
    use rust_htslib::bam::Read;
    
    let mut bam = bio_io::bam_reader("tests/data/msp_nuc.bam");
    if let Some(Ok(record)) = bam.records().next() {
        record
    } else {
        // Fallback: create a minimal record if test file not available
        bam::Record::new()
    }
}

/// Generate test data of various sizes for benchmarking
fn generate_test_data(num_features: usize) -> (Vec<i64>, Option<Vec<i64>>, Option<Vec<i64>>) {
    let mut starts = Vec::with_capacity(num_features);
    let mut lengths = Vec::with_capacity(num_features);
    
    // Generate realistic fiber-seq feature positions
    let mut pos = 100i64;
    for i in 0..num_features {
        starts.push(pos);
        
        // Vary feature lengths realistically (5-500bp range)
        let length = match i % 4 {
            0 => 5 + (i % 10) as i64,      // Small features (5-15bp) - m6A sites
            1 => 30 + (i % 20) as i64,     // Medium features (30-50bp) - MSPs  
            2 => 100 + (i % 50) as i64,    // Large features (100-150bp) - Nucleosomes
            _ => 200 + (i % 300) as i64,   // Very large features (200-500bp) - Large domains
        };
        
        lengths.push(length);
        pos += length + 10 + (i % 20) as i64; // Add gaps between features
    }
    
    (starts, None, Some(lengths))
}

/// Benchmark Ranges creation with different feature counts
fn bench_ranges_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("ranges_creation");
    group.measurement_time(Duration::from_secs(10));
    
    let feature_counts = vec![1, 5, 10, 50, 100, 500, 1000];
    
    for &count in &feature_counts {
        let record = get_test_bam_record();
        let (starts, ends, lengths) = generate_test_data(count);
        
        group.bench_with_input(
            BenchmarkId::new("new", count),
            &count,
            |b, _| {
                b.iter(|| {
                    let ranges = Ranges::new(
                        black_box(&record),
                        black_box(starts.clone()),
                        black_box(ends.clone()),
                        black_box(lengths.clone()),
                    );
                    black_box(ranges)
                });
            },
        );
    }
    group.finish();
}

/// Benchmark getter method performance
fn bench_getter_methods(c: &mut Criterion) {
    let mut group = c.benchmark_group("getter_methods");
    
    // Create ranges with moderate size
    let record = get_test_bam_record();
    let (starts, ends, lengths) = generate_test_data(100);
    let ranges = Ranges::new(&record, starts, ends, lengths);
    
    group.bench_function("starts", |b| {
        b.iter(|| {
            let result = black_box(&ranges).starts();
            black_box(result)
        });
    });
    
    group.bench_function("ends", |b| {
        b.iter(|| {
            let result = black_box(&ranges).ends();
            black_box(result)
        });
    });
    
    group.bench_function("lengths", |b| {
        b.iter(|| {
            let result = black_box(&ranges).lengths();
            black_box(result)
        });
    });
    
    group.bench_function("reference_starts", |b| {
        b.iter(|| {
            let result = black_box(&ranges).reference_starts();
            black_box(result)
        });
    });
    
    group.bench_function("get_starts", |b| {
        b.iter(|| {
            let result = black_box(&ranges).get_starts();
            black_box(result)
        });
    });
    
    group.bench_function("get_molecular", |b| {
        b.iter(|| {
            let result = black_box(&ranges).get_molecular();
            black_box(result)
        });
    });
    
    group.bench_function("get_reference", |b| {
        b.iter(|| {
            let result = black_box(&ranges).get_reference();
            black_box(result)
        });
    });
    
    group.finish();
}

/// Benchmark iteration performance
fn bench_iteration(c: &mut Criterion) {
    let mut group = c.benchmark_group("iteration");
    
    let feature_counts = vec![10, 100, 1000];
    
    for &count in &feature_counts {
        let record = get_test_bam_record();
        let (starts, ends, lengths) = generate_test_data(count);
        let ranges = Ranges::new(&record, starts, ends, lengths);
        
        group.bench_with_input(
            BenchmarkId::new("for_loop_iteration", count),
            &ranges,
            |b, ranges| {
                b.iter(|| {
                    let mut sum = 0i64;
                    for (start, end, length, qual, reference) in black_box(ranges) {
                        sum += start + end + length + qual as i64;
                        if let Some((ref_start, ref_end, ref_length)) = reference {
                            sum += ref_start + ref_end + ref_length;
                        }
                    }
                    black_box(sum)
                });
            },
        );
    }
    
    group.finish();
}

/// Benchmark filtering operations
fn bench_filtering(c: &mut Criterion) {
    let mut group = c.benchmark_group("filtering");
    
    let record = get_test_bam_record();
    let (starts, ends, lengths) = generate_test_data(500);
    
    group.bench_function("filter_by_qual", |b| {
        b.iter(|| {
            let mut ranges = Ranges::new(black_box(&record), black_box(starts.clone()), black_box(ends.clone()), black_box(lengths.clone()));
            // Set some quality scores
            let qual_scores = (0..500).map(|i| (i % 256) as u8).collect();
            ranges.set_qual(qual_scores);
            
            ranges.filter_by_qual(black_box(100));
            black_box(ranges)
        });
    });
    
    group.bench_function("filter_starts_at_read_ends", |b| {
        b.iter(|| {
            let mut ranges = Ranges::new(black_box(&record), black_box(starts.clone()), black_box(ends.clone()), black_box(lengths.clone()));
            ranges.filter_starts_at_read_ends(black_box(100));
            black_box(ranges)
        });
    });
    
    group.finish();
}

/// Benchmark cloning performance
fn bench_cloning(c: &mut Criterion) {
    let mut group = c.benchmark_group("cloning");
    
    let feature_counts = vec![10, 100, 1000];
    
    for &count in &feature_counts {
        let record = get_test_bam_record();
        let (starts, ends, lengths) = generate_test_data(count);
        let ranges = Ranges::new(&record, starts, ends, lengths);
        
        group.bench_with_input(
            BenchmarkId::new("clone", count),
            &ranges,
            |b, ranges| {
                b.iter(|| {
                    let cloned = black_box(ranges).clone();
                    black_box(cloned)
                });
            },
        );
    }
    
    group.finish();
}

/// Benchmark string serialization
fn bench_to_strings(c: &mut Criterion) {
    let mut group = c.benchmark_group("to_strings");
    
    let record = get_test_bam_record();
    let (starts, ends, lengths) = generate_test_data(100);
    let ranges = Ranges::new(&record, starts, ends, lengths);
    
    group.bench_function("to_strings_reference", |b| {
        b.iter(|| {
            let result = black_box(&ranges).to_strings(black_box(true), black_box(false));
            black_box(result)
        });
    });
    
    group.bench_function("to_strings_molecular", |b| {
        b.iter(|| {
            let result = black_box(&ranges).to_strings(black_box(false), black_box(false));
            black_box(result)
        });
    });
    
    group.finish();
}

/// Benchmark merge_ranges performance
fn bench_merge_ranges(c: &mut Criterion) {
    let mut group = c.benchmark_group("merge_ranges");
    
    let merge_counts = vec![2, 5, 10, 20];
    
    for &count in &merge_counts {
        let mut ranges_vec = Vec::new();
        
        for i in 0..count {
            let record = get_test_bam_record();
            let (starts, ends, lengths) = generate_test_data(50 + i * 10);
            let ranges = Ranges::new(&record, starts, ends, lengths);
            ranges_vec.push(ranges);
        }
        
        group.bench_with_input(
            BenchmarkId::new("merge", count),
            &ranges_vec,
            |b, ranges_vec| {
                b.iter(|| {
                    let refs: Vec<&Ranges> = ranges_vec.iter().collect();
                    let merged = Ranges::merge_ranges(black_box(refs));
                    black_box(merged)
                });
            },
        );
    }
    
    group.finish();
}

/// Benchmark memory usage patterns
fn bench_memory_patterns(c: &mut Criterion) {
    let mut group = c.benchmark_group("memory_patterns");
    
    // Test different allocation patterns
    group.bench_function("many_small_ranges", |b| {
        b.iter(|| {
            let mut ranges_vec = Vec::new();
            let record = get_test_bam_record();
            
            for _ in 0..100 {
                let (starts, ends, lengths) = generate_test_data(5); // Small ranges
                let ranges = Ranges::new(black_box(&record), starts, ends, lengths);
                ranges_vec.push(ranges);
            }
            black_box(ranges_vec)
        });
    });
    
    group.bench_function("few_large_ranges", |b| {
        b.iter(|| {
            let mut ranges_vec = Vec::new();
            let record = get_test_bam_record();
            
            for _ in 0..10 {
                let (starts, ends, lengths) = generate_test_data(500); // Large ranges
                let ranges = Ranges::new(black_box(&record), starts, ends, lengths);
                ranges_vec.push(ranges);
            }
            black_box(ranges_vec)
        });
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_ranges_creation,
    bench_getter_methods,
    bench_iteration,
    bench_filtering,
    bench_cloning,
    bench_to_strings,
    bench_merge_ranges,
    bench_memory_patterns
);

criterion_main!(benches);