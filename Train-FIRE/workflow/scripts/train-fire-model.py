#!/usr/bin/env python
"""
Train a FIRE XGBoost classifier via mokapot and emit:
  <outdir>/FIRE.xgb.bin
  <outdir>/FIRE.gbdt.json
  <outdir>/FIRE.conf.json
  <outdir>/FIRE.FDR.pdf
  <outdir>/FIRE.feature.importance.pdf
  <outdir>/metrics.json

Adapted from Train-FIRE/train-fire-model.py with CLI-driven hyperparameters
and a configurable output directory so multiple experiments can run in parallel.
"""
from __future__ import print_function

import argparse
import ast
import json
import logging
import os
import struct
import sys
from pathlib import Path
from typing import List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mokapot
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import GridSearchCV
from xgboost import XGBClassifier


RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


def parse_list(s: str, cast=float) -> List:
    """Parse a CLI-supplied list like '[0.5, 1]', '0.5,1', or '0.5 1' into [cast(x), ...]."""
    if s is None:
        return []
    s = s.strip()
    if s.startswith("["):
        return [cast(x) for x in ast.literal_eval(s)]
    if "," in s:
        return [cast(x) for x in s.split(",") if x.strip() != ""]
    return [cast(x) for x in s.split() if x.strip() != ""]


def convert_to_gbdt(input_model: str, output_file: str) -> int:
    model = xgb.Booster()
    model.load_model(input_model)
    tmp_file = output_file + ".gbdt_rs.mid"
    try:
        with open(input_model, "rb") as f:
            model_format = struct.unpack("cccc", f.read(4))
            model_format = b"".join(model_format)
            if model_format == b"bs64":
                raise RuntimeError("bs64 model format not supported")
            elif model_format != b"binf":
                f.seek(0)
            base_score = struct.unpack("f", f.read(4))[0]
    except Exception as e:
        logging.error(f"Reading base_score failed: {e}")
        return 1

    if os.path.exists(tmp_file):
        os.remove(tmp_file)
    model.dump_model(tmp_file, dump_format="json")
    with open(output_file, "w") as f:
        f.write(repr(base_score) + "\n")
        with open(tmp_file) as f2:
            f.write(f2.read())
    os.remove(tmp_file)
    return 0


def save_model(model, test_conf, outdir: Path, max_fdr: float) -> dict:
    conf_df = test_conf.confidence_estimates["psms"]
    conf_df = pd.concat([conf_df, test_conf.decoy_confidence_estimates["psms"]])
    simple = (
        conf_df[["mokapot score", "mokapot q-value"]]
        .drop_duplicates()
        .sort_values(by=["mokapot score", "mokapot q-value"])
    )
    simple.loc[simple["mokapot q-value"] > max_fdr, "mokapot q-value"] = 1.0
    g = simple.sort_values(["mokapot q-value", "mokapot score"]).groupby(
        "mokapot q-value"
    )
    simple = (
        pd.concat([g.head(1), g.tail(1)])
        .drop_duplicates()
        .sort_values("mokapot score")
        .reset_index(drop=True)
    )

    model.estimator.save_model(str(outdir / "FIRE.xgb.bin"))
    convert_to_gbdt(str(outdir / "FIRE.xgb.bin"), str(outdir / "FIRE.gbdt.json"))
    (outdir / "FIRE.conf.json").write_text(simple.to_json(orient="split", index=False))

    model.estimator.get_booster().feature_names = model.features
    ax = xgb.plot_importance(model.estimator)
    ax.figure.set_size_inches(10, 8)
    ax.figure.tight_layout()
    ax.figure.savefig(outdir / "FIRE.feature.importance.pdf")
    plt.close(ax.figure)

    _fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    test_conf.plot_qvalues(c="#24B8A0", ax=ax, label="mokapot", threshold=max_fdr * 2)
    ax.axvline(max_fdr, color="#343131", linestyle="--")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outdir / "FIRE.FDR.pdf")
    plt.close()

    psms = test_conf.confidence_estimates["psms"]
    n_at_fdr = int((psms["mokapot q-value"] <= max_fdr).sum())
    return {
        "n_test_psms": int(len(psms)),
        "n_test_psms_at_fdr": n_at_fdr,
        "max_fdr": max_fdr,
    }


def train_classifier(train_df, test_df, args, scale_pos_weight):
    if args.grid_search:
        mcw = (len(train_df) * np.array(args.min_child_weight_fracs)).astype(int)
        grid = {
            "n_estimators": args.n_estimators_grid,
            "scale_pos_weight": [scale_pos_weight],
            "max_depth": args.max_depth_grid,
            "min_child_weight": mcw.tolist(),
            "colsample_bytree": args.colsample_bytree_grid,
            "gamma": args.gamma_grid,
        }
        xgb_model = GridSearchCV(
            XGBClassifier(
                use_label_encoder=False,
                eval_metric="auc",
                seed=RANDOM_SEED,
                n_jobs=args.inner_jobs,
            ),
            param_grid=grid,
            cv=5,
            scoring="roc_auc",
            verbose=2,
            n_jobs=args.outer_jobs,
        )
    else:
        # no grid search -> give all threads to XGBoost
        xgb_model = XGBClassifier(
            use_label_encoder=False,
            eval_metric="auc",
            n_estimators=args.n_estimators_grid[0],
            max_depth=args.max_depth_grid[0],
            min_child_weight=int(len(train_df) * args.min_child_weight_fracs[0]),
            gamma=args.gamma_grid[0],
            colsample_bytree=args.colsample_bytree_grid[0],
            scale_pos_weight=scale_pos_weight,
            seed=RANDOM_SEED,
            n_jobs=args.outer_jobs * args.inner_jobs,
        )

    train_psms = mokapot.read_pin(train_df)
    model = mokapot.Model(
        xgb_model,
        train_fdr=args.train_fdr,
        subset_max_train=args.subset_max_train,
        direction=args.direction,
        max_iter=15,
    )
    model.fit(train_psms)
    test_psms = mokapot.read_pin(test_df)
    scores = model.predict(test_psms)
    test_conf = test_psms.assign_confidence(scores)
    return model, test_conf


def balance_df(df):
    counts = df.Label.value_counts()
    return df.groupby("Label").apply(lambda g: g.sample(min(counts))).reset_index(drop=True)


def read_features(infile, args):
    df = pd.read_csv(infile, sep="\t")
    df.insert(0, "SpecId", df.index)
    df["Peptide"] = df.SpecId
    df["Proteins"] = df.SpecId
    df["scannr"] = df.SpecId
    assert "Label" in df.columns
    logging.info(f"Label counts raw: {df.Label.value_counts().to_dict()}")

    df = df[(df.msp_len >= args.min_msp_length_for_positive_fire_call) | (df.Label == -1)]
    df = df[(df.msp_len >= args.min_msp_length_for_negative_fire_call) | (df.Label == 1)]
    df = df.groupby(["fiber", "Label"]).sample(n=1).reset_index(drop=True)

    for col in df.columns:
        if "AT" in col or "rle" in col:
            df[col] = 1

    df.drop(columns=["#chrom", "start", "end", "fiber"], inplace=True)
    r = np.random.rand(len(df))
    train = df[r < 0.80]
    test = df[r >= 0.80]
    train = train[
        (train.msp_len >= args.min_msp_length_for_positive_fire_call) | (train.Label == 1)
    ]
    return balance_df(train), balance_df(test)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("training_data")
    ap.add_argument("--outdir", required=True)
    ap.add_argument(
        "--outer-jobs",
        type=int,
        default=1,
        help="joblib workers for GridSearchCV (one process per fold/grid-point fit).",
    )
    ap.add_argument(
        "--inner-jobs",
        type=int,
        default=8,
        help="XGBoost n_jobs (OpenMP threads per single fit).",
    )
    ap.add_argument("--train-fdr", type=float, default=0.05)
    ap.add_argument("--test-fdr", type=float, default=0.05)
    ap.add_argument("--subset-max-train", type=int, default=2_000_000)
    ap.add_argument("--direction", default="msp_len_times_m6a_fc")
    ap.add_argument("--min-msp-length-for-positive-fire-call", type=int, default=85)
    ap.add_argument("--min-msp-length-for-negative-fire-call", type=int, default=85)
    ap.add_argument("--grid-search", action="store_true")
    ap.add_argument("--n-estimators-grid", default="[200, 300]")
    ap.add_argument("--max-depth-grid", default="[9, 15]")
    ap.add_argument("--min-child-weight-fracs", default="[0.001, 0.005]")
    ap.add_argument("--colsample-bytree-grid", default="[0.5, 1.0]")
    ap.add_argument("--gamma-grid", default="[1]")
    args = ap.parse_args()

    args.n_estimators_grid = parse_list(args.n_estimators_grid, int)
    args.max_depth_grid = parse_list(args.max_depth_grid, int)
    args.min_child_weight_fracs = parse_list(args.min_child_weight_fracs, float)
    args.colsample_bytree_grid = parse_list(args.colsample_bytree_grid, float)
    args.gamma_grid = parse_list(args.gamma_grid, float)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        format="[%(levelname)s][%(relativeCreated)d ms]: %(message)s",
        level=logging.INFO,
    )

    train_df, test_df = read_features(args.training_data, args)
    n_pos = int((train_df.Label == 1).sum())
    n_neg = int((train_df.Label == -1).sum())
    scale_pos_weight = (n_neg / max(n_pos, 1)) if n_pos else 1.0

    model, test_conf = train_classifier(train_df, test_df, args, scale_pos_weight)

    metrics = save_model(model, test_conf, outdir, max_fdr=args.test_fdr)
    metrics.update(
        dict(
            experiment=outdir.name,
            n_train=int(len(train_df)),
            n_test=int(len(test_df)),
            n_train_pos=n_pos,
            n_train_neg=n_neg,
            scale_pos_weight=float(scale_pos_weight),
            train_fdr=args.train_fdr,
            test_fdr=args.test_fdr,
            direction=args.direction,
            grid_search=bool(args.grid_search),
            outer_jobs=int(args.outer_jobs),
            inner_jobs=int(args.inner_jobs),
        )
    )
    if args.grid_search and hasattr(model.estimator, "best_params_"):
        metrics["best_params"] = model.estimator.best_params_
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str))
    logging.info(f"Wrote metrics: {metrics}")


if __name__ == "__main__":
    main()
