#!/usr/bin/env python
from __future__ import print_function
import defopt
import sys
import os
import logging
from pathlib import Path
from typing import Optional
import mokapot  # version 0.9
from xgboost import XGBClassifier  # version 0.82 # required for rust package
import pandas as pd
import numpy as np
import xgboost as xgb
import struct
from ctypes import cdll
from ctypes import c_float, c_uint, c_char_p, c_bool
import matplotlib.pyplot as plt


# set random seed
RANDOM_SEED = 42
np.random.seed(RANDOM_SEED)


def convert_to_gbdt(input_model, objective, output_file):
    """
    supported objectives: 'reg:linear', 'binary:logistic', 'reg:logistic',
        'binary:logitraw', 'multi:softmax', 'multi:softprob', 'rank:pairwise'
    """
    model = xgb.Booster()
    model.load_model(input_model)
    tmp_file = output_file + ".gbdt_rs.mid"
    # extract base score
    try:
        with open(input_model, "rb") as f:
            model_format = struct.unpack("cccc", f.read(4))
            model_format = b"".join(model_format)
            if model_format == b"bs64":
                print("This model type is not supported")
            elif model_format != "binf":
                f.seek(0)
            base_score = struct.unpack("f", f.read(4))[0]
    except Exception as e:
        print("error: ", e)
        return 1

    if os.path.exists(tmp_file):
        print(
            "Intermediate file %s exists. Please remove this file or change your output file path"
            % tmp_file
        )
        return 1

    # dump json
    model.dump_model(tmp_file, dump_format="json")

    # add base score to json file
    try:
        with open(output_file, "w") as f:
            f.write(repr(base_score) + "\n")
            with open(tmp_file) as f2:
                for line in f2.readlines():
                    f.write(line)
    except Exception as e:
        print("error: ", e)
        os.remove(tmp_file)
        return 1

    os.remove(tmp_file)
    return 0


def save_mokapot_model_for_fibertools(model, test_conf, max_fdr=0.10):
    conf_df = test_conf.confidence_estimates["psms"]
    conf_df = pd.concat([conf_df, test_conf.decoy_confidence_estimates["psms"]])
    logging.info(f"Test FDR:\n{conf_df}")
    simple_conf_df = (
        conf_df[["mokapot score", "mokapot q-value"]]
        .drop_duplicates()
        .sort_values(by=["mokapot score", "mokapot q-value"])
    )
    high_fdr = simple_conf_df["mokapot q-value"] > max_fdr
    simple_conf_df.loc[high_fdr, "mokapot q-value"] = 1.0

    # keep only the first and last score of each q-value
    g = simple_conf_df.sort_values(["mokapot q-value", "mokapot score"]).groupby(
        "mokapot q-value"
    )
    simple_conf_df = (
        pd.concat([g.head(1), g.tail(1)])
        .drop_duplicates()
        .sort_values("mokapot score")
        .reset_index(drop=True)
    )
    logging.info(f"Test FDR:\n{simple_conf_df}")

    # test_conf.to_txt("tmp_mokapot_results")
    model.estimator.save_model("FIRE.xgb.bin")
    convert_to_gbdt("FIRE.xgb.bin", "binary:logistic", "FIRE.gbdt.json")
    with open("FIRE.conf.json", "w") as f:
        f.write(simple_conf_df.to_json(orient="split", index=False))

    # save the important features
    model.estimator.get_booster().feature_names = model.features
    ax = xgb.plot_importance(model.estimator)
    ax.figure.set_size_inches(10, 8)
    ax.figure.tight_layout()
    ax.figure.savefig("FIRE.feature.importance.pdf")

    # save the fdr curve
    _fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    colors = ("#343131", "#24B8A0")
    test_conf.plot_qvalues(c=colors[1], ax=ax, label="mokapot", threshold=max_fdr * 2)
    # draw a vertial line at the max_fdr
    ax.axvline(max_fdr, color=colors[0], linestyle="--")
    ax.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"FIRE.FDR.pdf")


# grid = {
#    "n_estimators": [50, 100, 150],
#    "scale_pos_weight": [scale_pos_weight],  # [0.5, 1, 2], #np.logspace(0, 2, 3),
#    "max_depth": [3, 6, 9],
#    "min_child_weight": [3, 6, 9, 12],
#    "gamma": [0.1, 1, 10],
# }
def train_classifier(
    train_df,
    test_df,
    min_msp_size=1,
    subset_max_train=200_000,
    test_fdr=0.05,
    train_fdr=0.1,
    n_estimators=100,  # from previous grid search, may need to use 50
    max_depth=6,  # from previous grid search
    min_child_weight=9,  # from previous grid search
    gamma=10,  # from previous grid search
):
    train_df = train_df[train_df.msp_len >= min_msp_size]
    scale_pos_weight = sum(train_df.Label == -1) / sum(train_df.Label == 1)
    logging.info(f"scale_pos_weight: {scale_pos_weight}")

    xgb_model = XGBClassifier(
        use_label_encoder=False,
        eval_metric="auc",
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_child_weight=min_child_weight,
        gamma=gamma,
        scale_pos_weight=scale_pos_weight,
        seed=RANDOM_SEED,
    )

    # train_psms = mokapot.read_pin("train.pin")
    train_psms = mokapot.read_pin(train_df)
    # model = mokapot.PercolatorModel()
    model = mokapot.Model(
        xgb_model,
        train_fdr=train_fdr,
        subset_max_train=subset_max_train,
        direction="log_msp_len",
    )
    model.fit(train_psms)

    # test_psms = mokapot.read_pin("test.pin")
    test_psms = mokapot.read_pin(test_df)
    scores = model.predict(test_psms)
    test_conf = test_psms.assign_confidence(scores)
    return (model, test_conf)


def read_input_features(
    infile, min_msp_length_for_positive_fire_call, test_train_split=0.8
):
    df = pd.read_csv(infile, sep="\t")
    logging.info(f"Columns: {df.columns}")
    # add columns needed for mokapot
    df.insert(0, "SpecId", df.index)
    df["Peptide"] = df.SpecId
    df["Proteins"] = df.SpecId
    df["scannr"] = df.SpecId
    # check the labels
    assert "Label" in df.columns
    logging.info(
        f"Label counts before setting a minimum FIRE size: {df.Label.value_counts()}"
    )
    # df.loc[df.msp_len < min_msp_length_for_positive_fire_call, "Label"] = -1
    df = df[df.msp_len >= min_msp_length_for_positive_fire_call]
    logging.info(
        f"Label counts after setting a minimum FIRE size: {df.Label.value_counts()}"
    )

    df.drop(columns=["#chrom", "start", "end", "fiber"], inplace=True)
    # make random test and train in 90 10 split
    random = np.random.rand(len(df))
    train = df[random < test_train_split]
    test = df[random >= test_train_split]
    logging.info(f"Train size: {train.shape}")
    logging.info(f"Test size: {test.shape}")
    return (train, test)


def main(
    infile: Path,
    *,
    min_msp_size: int = 1,
    subset_max_train: int = 200_000,
    test_fdr: float = 0.05,
    train_fdr: float = 0.1,
    n_estimators: int = 100,
    max_depth: int = 6,
    min_child_weight: int = 6,
    min_msp_length_for_positive_fire_call: int = 85,
    gamma: float = 1,
    verbose: int = 1,
):
    """
    Author Mitchell R. Vollger

    :param infile: Input file, stdin by default
    :param outfile: Output file, stdout by default
    :param count: Number of times to display the greeting
    :param verbose: Set the logging level of the function
    """
    logger = logging.getLogger()
    log_format = "[%(levelname)s][Time elapsed (ms) %(relativeCreated)d]: %(message)s"
    log_level = 10 * (3 - verbose)
    logging.basicConfig(format=log_format)
    logger.setLevel(log_level)

    train_df, test_df = read_input_features(
        infile, min_msp_length_for_positive_fire_call
    )
    model, test_conf = train_classifier(
        train_df,
        test_df,
        min_msp_size=min_msp_size,
        subset_max_train=subset_max_train,
        test_fdr=test_fdr,
        train_fdr=train_fdr,
        n_estimators=n_estimators,
        max_depth=max_depth,
        min_child_weight=min_child_weight,
        gamma=gamma,
    )
    save_mokapot_model_for_fibertools(model, test_conf, max_fdr=test_fdr)
    return (model, test_conf)


if __name__ == "__main__":
    model, test_conf = defopt.run(main, show_types=True, version="0.0.1")
