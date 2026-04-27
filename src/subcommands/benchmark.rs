use crate::cli::BenchmarkOptions;
use anyhow::Result;
use std::time::Instant;

pub fn run_benchmark(opts: &mut BenchmarkOptions) -> Result<()> {
    log::info!("Starting fiber iterator benchmark");

    let mut bam = opts.input.bam_reader();

    // Benchmark input.fibers()
    let start_time = Instant::now();
    let mut fiber_count = 0;

    for fiber in opts.input.fibers(&mut bam) {
        fiber_count += 1;
        // Force materialization but do nothing with the fiber
        std::hint::black_box(&fiber);
    }

    let elapsed = start_time.elapsed();

    log::info!(
        "Processed {} fibers in {:.3} seconds",
        fiber_count,
        elapsed.as_secs_f64()
    );
    log::info!(
        "Rate: {:.0} fibers/second",
        fiber_count as f64 / elapsed.as_secs_f64()
    );

    Ok(())
}
