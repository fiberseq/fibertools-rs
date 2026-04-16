import os
import shlex


def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


def _expand(p):
    return os.path.expanduser(os.path.expandvars(str(p)))


FAI = _expand(config["reference"]["fai"])
FASTA = _expand(config["reference"]["fasta"])
GENOME = config["reference"]["genome"]

S3_ENDPOINT = config.get("s3_endpoint", None)
SAMPLE_RATE = float(config["sample_rate"])


def is_remote(p):
    s = str(p)
    return s.startswith(("s3://", "http://", "https://"))


def resolve_alignment(src):
    """
    Convert an s3://bucket/key URL into https://<endpoint-host>/bucket/key
    using the configured s3_endpoint. htslib can read the resulting HTTPS
    URL directly via libcurl with no auth, which works for public buckets
    where its S3 transport fails ('Resource temporarily unavailable').
    Local paths and existing http(s) URLs pass through unchanged.
    """
    s = _expand(src)
    if s.startswith("s3://"):
        if not S3_ENDPOINT:
            raise ValueError("training_bam is s3:// but config.s3_endpoint is not set")
        return S3_ENDPOINT.rstrip("/") + "/" + s[len("s3://") :]
    return s


TRAINING_BAM = resolve_alignment(config["training_bam"])
TEST_BAM = resolve_alignment(config["test_bam"] or config["training_bam"])


def aln_input(src):
    """If local, track as a snakemake input; if remote (http/https), skip."""
    return [] if is_remote(src) else [src]


TEST_REGION = str(config["test_region"])
N_SITES = int(config["n_sites"])
MERGE_DIST = int(config["merge_distance"])
EXCLUDE_CHROMS_RE = config["exclude_chroms_regex"]
EXCLUDE_BEDS = [_expand(p) for p in config["exclude_beds"]]

TRAIN_DEFAULTS = config["train_defaults"]


def exp_cfg(exp):
    return config["experiments"][exp]


def exp_positive_specs(exp):
    """
    Return list of (path, awk_filter_or_None) for an experiment's positive beds.
    Accepts either bare strings or {path:, awk_filter:} dicts in config.
    """
    out = []
    for item in exp_cfg(exp)["positive_beds"]:
        if isinstance(item, str):
            out.append((_expand(item), None))
        else:
            out.append((_expand(item["path"]), item.get("awk_filter")))
    return out


def exp_positive_paths(exp):
    return [p for p, _ in exp_positive_specs(exp)]


def positive_source_cmds(exp):
    """
    Bash snippet that concatenates every positive bed, applying any
    per-source awk_filter, and emits 3-column output. Uses `zcat -f`
    so both gzipped and plain-text beds work transparently.
    """
    parts = []
    for path, awkf in exp_positive_specs(exp):
        q = shlex.quote(path)
        if awkf:
            parts.append(f"zcat -f {q} | awk {shlex.quote(awkf)} | cut -f1-3")
        else:
            parts.append(f"zcat -f {q} | cut -f1-3")
    return "; ".join(parts)


def exp_negative_exclude_paths(exp):
    """Extra beds whose regions are excluded from negative sampling but not used as positives."""
    out = []
    for item in exp_cfg(exp).get("negative_exclude_beds", []) or []:
        out.append(_expand(item) if isinstance(item, str) else _expand(item["path"]))
    return out


def exp_neg_mask_paths(exp):
    """Union of positives + negative_exclude_beds used to keep negatives off signal."""
    return exp_positive_paths(exp) + exp_negative_exclude_paths(exp)


def exp_train_params(exp):
    """Merge train_defaults with per-experiment overrides."""
    params = dict(TRAIN_DEFAULTS)
    params.update(exp_cfg(exp).get("train", {}) or {})
    return params


def all_exp_region_files():
    """Union of regions across experiments (positives + negatives) needed for feature extraction."""
    paths = []
    for exp in config["experiments"].keys():
        paths.append(f"results/experiments/{exp}/positives.bed.gz")
        paths.append(f"results/experiments/{exp}/negatives.bed.gz")
    return paths
