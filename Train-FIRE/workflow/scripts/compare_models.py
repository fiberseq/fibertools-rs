#!/usr/bin/env python
"""Aggregate per-experiment FIRE.conf.json + metrics.json into one plot and TSV."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="results/experiments dir")
    ap.add_argument("--experiments", nargs="+", required=True)
    ap.add_argument("--out-pdf", required=True)
    ap.add_argument("--out-tsv", required=True)
    args = ap.parse_args()

    root = Path(args.root)
    Path(args.out_pdf).parent.mkdir(parents=True, exist_ok=True)

    rows = []
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))

    for exp in args.experiments:
        conf_path = root / exp / "FIRE.conf.json"
        m_path = root / exp / "metrics.json"
        if not conf_path.exists() or not m_path.exists():
            continue
        metrics = json.loads(m_path.read_text())
        rows.append(metrics)

        conf = json.loads(conf_path.read_text())
        df = pd.DataFrame(conf["data"], columns=conf["columns"])
        df = df.sort_values("mokapot score")
        ax.plot(df["mokapot q-value"], range(len(df), 0, -1), label=exp)

    ax.set_xlabel("q-value")
    ax.set_ylabel("Cumulative PSMs")
    ax.set_xscale("log")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("FIRE model FDR curves")
    plt.tight_layout()
    plt.savefig(args.out_pdf)
    plt.close()

    pd.DataFrame(rows).to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()
