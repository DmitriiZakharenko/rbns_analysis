#!/usr/bin/env python3
"""
Validate ML dataset: basic checks, label distribution, optional motif check.
Output: results/validation_summary_rbns.tsv, results/dataset_stats_rbns.json
"""

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.io import load_config, get_paths, read_tsv, write_json, write_tsv


def main():
    parser = argparse.ArgumentParser(description="Validate ML dataset")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--input", help="Cleaned dataset TSV")
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    results_dir = paths["results"]
    in_path = Path(args.input or results_dir / "ml_dataset_rbns_clean.tsv")

    if not in_path.exists():
        print(f"Dataset not found: {in_path}")
        return 1

    import pandas as pd
    df = pd.DataFrame(read_tsv(in_path))

    if df.empty:
        print("Empty dataset")
        return 1

    # binding_label is read as string from TSV — cast to int for comparisons
    df["binding_label"] = pd.to_numeric(df["binding_label"], errors="coerce").fillna(-1).astype(int)

    stats = {}
    issues = []

    # Basic checks
    for col in ["target_name", "rna_sequence", "binding_label"]:
        if col not in df.columns:
            issues.append(f"Missing column: {col}")
        else:
            missing = df[col].isna().sum()
            if missing > 0:
                issues.append(f"Missing values in {col}: {missing}")

    n_total = len(df)
    stats["n_total"] = n_total
    stats["n_positive"] = int((df["binding_label"] == 1).sum())
    stats["n_negative"] = int((df["binding_label"] == 0).sum())
    stats["n_rbp"] = df["target_name"].nunique()

    # Pos/neg overlap: no (target, seq) should appear with both labels
    overlap = df.groupby(["target_name", "rna_sequence"])["binding_label"].nunique()
    n_overlap = int((overlap > 1).sum())
    if n_overlap > 0:
        issues.append(f"Pos/neg overlap (same target+seq, both labels): {n_overlap}")
    stats["n_overlap"] = n_overlap

    # Duplicates
    dup = df.duplicated(subset=["target_name", "rna_sequence"]).sum()
    if dup > 0:
        issues.append(f"Duplicate (target_name, rna_sequence): {dup}")
    stats["n_duplicates"] = int(dup)

    # Sequence length
    if "rna_sequence" in df.columns:
        lengths = df["rna_sequence"].astype(str).str.len()
        stats["seq_length_min"] = int(lengths.min())
        stats["seq_length_max"] = int(lengths.max())
        stats["seq_length_mean"] = float(lengths.mean())
        bad_alpha = df["rna_sequence"].astype(str).str.contains(r"[^ACGU]").sum()
        if bad_alpha > 0:
            issues.append(f"Sequences with non-ACGU characters: {bad_alpha}")
        stats["n_bad_alphabet"] = int(bad_alpha)

    # Per-RBP summary
    per_rbp = df.groupby("target_name").agg(
        total=("binding_label", "count"),
        positive=("binding_label", lambda s: (s == 1).sum()),
        negative=("binding_label", lambda s: (s == 0).sum()),
    ).reset_index()
    summary_path = results_dir / "validation_summary_rbns.tsv"
    per_rbp.to_csv(summary_path, sep="\t", index=False)
    stats["per_rbp_path"] = str(summary_path)

    stats_path = results_dir / "dataset_stats_rbns.json"
    write_json(stats_path, stats)
    print("Stats:", json.dumps(stats, indent=2))
    if issues:
        print("Issues:", issues)
    else:
        print("No issues found.")
    print(f"Wrote {stats_path}, {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
