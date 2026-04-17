#!/usr/bin/env python3
"""
Trim protein_sequence in RBNS ML dataset to cloned fragment boundaries from mmc2.

Expected mmc2 format (Lambert et al. supplementary table):
  - row 3 (0-based index 2) contains headers
  - columns include:
      RBP
      start_pos (amino acids)
      stop_pos (amino acids)

Usage:
  python scripts/08_trim_protein_fragments.py \
    --mmc2 "1-s2.0-S1097276518303514-mmc2 (1).xlsx" \
    --input results/ml_dataset_rbns_clean.tsv \
    --output results/ml_dataset_rbns_clean.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_fragment_bounds(mmc2_path: Path) -> dict[str, tuple[int, int]]:
    raw = pd.read_excel(mmc2_path, sheet_name=0, header=None)
    if len(raw) < 4:
        raise ValueError("mmc2 table looks too short")

    header = raw.iloc[2].tolist()
    df = raw.iloc[3:].copy()
    df.columns = header

    required = ["RBP", "start_pos (amino acids)", "stop_pos (amino acids)"]
    for col in required:
        if col not in df.columns:
            raise ValueError(f"Missing required column in mmc2: {col}")

    out: dict[str, tuple[int, int]] = {}
    for _, row in df.iterrows():
        rbp = str(row["RBP"]).strip()
        if not rbp or rbp.lower() == "nan":
            continue
        try:
            start = int(float(row["start_pos (amino acids)"]))
            stop = int(float(row["stop_pos (amino acids)"]))
        except (TypeError, ValueError):
            continue
        if start < 1 or stop < start:
            continue
        out[rbp.upper()] = (start, stop)
    return out


def trim_dataset(df: pd.DataFrame, bounds: dict[str, tuple[int, int]]) -> tuple[pd.DataFrame, int, int]:
    if "target_name" not in df.columns or "protein_sequence" not in df.columns:
        raise ValueError("Input dataset must contain target_name and protein_sequence columns")

    mapped = 0
    changed = 0

    def _trim(row: pd.Series) -> str:
        nonlocal mapped, changed
        target = str(row["target_name"]).strip().upper()
        seq = str(row["protein_sequence"]).strip()
        b = bounds.get(target)
        if b is None:
            return seq
        mapped += 1
        start, stop = b
        trimmed = seq[start - 1 : stop]  # mmc2 is 1-based inclusive
        if trimmed != seq:
            changed += 1
        return trimmed

    df = df.copy()
    df["protein_sequence"] = df.apply(_trim, axis=1)
    return df, mapped, changed


def main() -> None:
    parser = argparse.ArgumentParser(description="Trim RBNS protein sequences to mmc2 fragment bounds")
    parser.add_argument("--mmc2", required=True, help="Path to mmc2 xlsx file")
    parser.add_argument("--input", default="results/ml_dataset_rbns_clean.tsv", help="Input dataset TSV")
    parser.add_argument("--output", default="results/ml_dataset_rbns_clean.tsv", help="Output dataset TSV")
    args = parser.parse_args()

    mmc2_path = Path(args.mmc2)
    in_path = Path(args.input)
    out_path = Path(args.output)

    if not mmc2_path.exists():
        raise FileNotFoundError(f"mmc2 file not found: {mmc2_path}")
    if not in_path.exists():
        raise FileNotFoundError(f"input dataset not found: {in_path}")

    bounds = load_fragment_bounds(mmc2_path)
    df = pd.read_csv(in_path, sep="\t")
    trimmed_df, mapped_rows, changed_rows = trim_dataset(df, bounds)

    trimmed_df.to_csv(out_path, sep="\t", index=False)
    n_targets_mapped = trimmed_df.loc[trimmed_df["target_name"].str.upper().isin(bounds.keys()), "target_name"].nunique()

    print(f"Loaded mmc2 bounds for {len(bounds)} RBPs")
    print(f"Rows mapped to mmc2 bounds: {mapped_rows}")
    print(f"Rows with changed protein_sequence: {changed_rows}")
    print(f"Targets mapped in dataset: {n_targets_mapped}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
