use super::common::{fixture, run};
use tempfile::NamedTempFile;

// Every MSP must produce a decoration: called FIREs by their precision and
// non-FIRE (precision-0) MSPs as LINKER, matching pre-MA output. FIRE quals
// live on the `fire` annotation type, so the decorator must overlay them
// onto the MSPs rather than iterate the fire type alone.
#[test]
fn track_decorators_emit_fire_and_linker() {
    let scored = NamedTempFile::with_suffix(".bam").unwrap();
    run(&[
        "fire",
        fixture("all.bam").to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let bed12 = NamedTempFile::with_suffix(".bed").unwrap();
    let out = run(&[
        "track-decorators",
        "--bed12",
        bed12.path().to_str().unwrap(),
        scored.path().to_str().unwrap(),
    ]);
    let count = |el: &str| {
        out.lines()
            .filter(|l| l.split('\t').any(|f| f == el))
            .count()
    };
    let (fire, linker) = (count("FIRE"), count("LINKER"));
    assert!(fire > 0, "no FIRE decorations emitted");
    assert!(
        linker > fire,
        "expected precision-0 MSPs (the majority) to decorate as LINKER; got {linker} LINKER vs {fire} FIRE"
    );
    insta::assert_snapshot!(out);
}
