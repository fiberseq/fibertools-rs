#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Mitchell R. Vollger
"""
Assemble a UCSC Track Hub with one decorator track per experiment.

Each experiment contributes a base bigBed 12+ track ("fire-fibers") plus
a decorator overlay ("fire-fiber-decorators"). Default visibility: squish.

Layout:
  <hub-dir>/
    hub.txt
    genomes.txt
    chrom.sizes
    <genome>/
      trackDb.txt
      bb/<exp>.fire-fibers.bb
      bb/<exp>.fire-fiber-decorators.bb
"""
import argparse
import shutil
from pathlib import Path


PALETTE = [
    (166, 54, 3), (217, 95, 14), (54, 144, 192), (34, 94, 168),
    (5, 112, 176), (35, 139, 69), (116, 196, 118), (49, 130, 189),
    (140, 81, 10), (128, 0, 128),
]


def color(i):
    return ",".join(str(c) for c in PALETTE[i % len(PALETTE)])


TRACK_TEMPLATE = """track {exp}
shortLabel {exp}
longLabel FIRE fibers, model={exp}
type bigBed 12 +
itemRgb on
visibility squish
color {color}
bigDataUrl bb/{exp}.fire-fibers.bb
decorator.default.bigDataUrl bb/{exp}.fire-fiber-decorators.bb
decorator.default.filterValues.keywords 5mC,m6A,NUC,LINKER,FIRE
decorator.default.filterValuesDefault.keywords LINKER,FIRE
"""


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
    ap.add_argument("--experiments", nargs="+", required=True)
    ap.add_argument("--results-root", required=True)
    ap.add_argument("--chrom-sizes", required=True, help="2-col chrom.sizes file")
    args = ap.parse_args()

    hub = Path(args.hub_dir)
    gdir = hub / args.genome
    (gdir / "bb").mkdir(parents=True, exist_ok=True)

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
        f"trackDb {args.genome}/trackDb.txt\n"
    )

    root = Path(args.results_root)
    blocks = []
    for i, exp in enumerate(args.experiments):
        for name in ("fire-fibers.bb", "fire-fiber-decorators.bb"):
            src = root / exp / name
            dst = gdir / "bb" / f"{exp}.{name}"
            shutil.copyfile(src, dst)
        blocks.append(TRACK_TEMPLATE.format(exp=exp, color=color(i)))

    (gdir / "trackDb.txt").write_text("\n".join(blocks))
    # convenience copy at top level (snakemake output path expects it here)
    shutil.copyfile(gdir / "trackDb.txt", hub / "trackDb.txt")


if __name__ == "__main__":
    main()
