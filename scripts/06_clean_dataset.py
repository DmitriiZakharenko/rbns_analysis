#!/usr/bin/env python3
"""
Clean ML dataset: drop missing, fix types, remove duplicates.
Output: results/ml_dataset_rbns_clean.tsv
"""

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.io import load_config, get_paths, read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser(description="Clean ML dataset")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--input", help="Input TSV (default: results/ml_dataset_rbns.tsv)")
    parser.add_argument("--output", help="Output TSV (default: results/ml_dataset_rbns_clean.tsv)")
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    results_dir = paths["results"]
    in_path = Path(args.input or results_dir / "ml_dataset_rbns.tsv")
    out_path = Path(args.output or results_dir / "ml_dataset_rbns_clean.tsv")

    if not in_path.exists():
        print(f"Input not found: {in_path}")
        return 1

    import pandas as pd
    df = pd.DataFrame(read_tsv(in_path))

    if df.empty:
        print("Empty dataset")
        return 1

    n0 = len(df)
    # Drop missing in key columns
    for col in ["target_name", "rna_sequence", "binding_label"]:
        if col in df.columns:
            df = df[df[col].notna() & (df[col].astype(str).str.strip() != "")]
    df["binding_label"] = pd.to_numeric(df["binding_label"], errors="coerce").fillna(0).astype(int)
    df = df[df["binding_label"].isin([0, 1])]
    # Remove duplicates by (target_name, rna_sequence)
    df = df.drop_duplicates(subset=["target_name", "rna_sequence"], keep="first")
    if "source" in df.columns:
        df["source"] = df["source"].fillna("")
    n1 = len(df)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Cleaned: {n0} -> {n1} rows. Written to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
