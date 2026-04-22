#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
"""
Assemble a UCSC Track Hub with one decorator track per model.

Flat layout (avoids double-nesting with hubCheck):
  <hub-dir>/
    hub.txt
    genomes.txt
    trackDb.txt
    chrom.sizes
    bb/<model>.fire-fibers.bb
    bb/<model>.fire-fiber-decorators.bb
"""
import argparse
import shutil
from pathlib import Path


BASELINE_COLOR = "0,0,0"
PALETTE = [
    (166, 54, 3), (217, 95, 14), (54, 144, 192), (34, 94, 168),
    (5, 112, 176), (35, 139, 69), (116, 196, 118), (49, 130, 189),
    (140, 81, 10), (128, 0, 128),
]


def color(i):
    return ",".join(str(c) for c in PALETTE[i % len(PALETTE)])


TRACK_TEMPLATE = """track {model}
shortLabel {short_label}
longLabel {long_label}
type bigBed 12 +
itemRgb on
visibility squish
color {color}
bigDataUrl bb/{model}.fire-fibers.bb
decorator.default.bigDataUrl bb/{model}.fire-fiber-decorators.bb
decorator.default.filterValues.keywords 5mC,m6A,NUC,LINKER,FIRE
decorator.default.filterValuesDefault.keywords LINKER,FIRE
"""


def render_block(model, is_baseline, palette_idx):
    if is_baseline:
        return TRACK_TEMPLATE.format(
            model=model,
            short_label="baseline (input CRAM)",
            long_label="FIRE fibers from the model baked into the input CRAM (no retraining)",
            color=BASELINE_COLOR,
        )
    return TRACK_TEMPLATE.format(
        model=model,
        short_label=model,
        long_label=f"FIRE fibers, trained model={model}",
        color=color(palette_idx),
    )


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--hub-dir", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--short-label", required=True)
    ap.add_argument("--long-label", required=True)
    ap.add_argument("--email", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--models", nargs="+", required=True)
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--chrom-sizes", required=True)
    args = ap.parse_args()

    hub = Path(args.hub_dir)
    (hub / "bb").mkdir(parents=True, exist_ok=True)

    shutil.copyfile(args.chrom_sizes, hub / "chrom.sizes")

    (hub / "hub.txt").write_text(
        f"hub {args.name}\n"
        f"shortLabel {args.short_label}\n"
        f"longLabel {args.long_label}\n"
        f"genomesFile genomes.txt\n"
        f"email {args.email}\n"
    )
    (hub / "genomes.txt").write_text(
        f"genome {args.genome}\n"
        f"trackDb trackDb.txt\n"
    )

    ordered = (["baseline"] if "baseline" in args.models else []) + [
        m for m in args.models if m != "baseline"
    ]

    root = Path(args.results_root)
    blocks = []
    palette_idx = 0
    for model in ordered:
        for name in ("fire-fibers.bb", "fire-fiber-decorators.bb"):
            shutil.copyfile(root / model / name, hub / "bb" / f"{model}.{name}")
        is_baseline = model == "baseline"
        blocks.append(render_block(model, is_baseline, palette_idx))
        if not is_baseline:
            palette_idx += 1

    (hub / "trackDb.txt").write_text("\n".join(blocks))


if __name__ == "__main__":
    main()
