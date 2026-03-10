#!/usr/bin/env python3
"""
Build ML dataset: merge positives + sampled negatives per RBP.
Output: results/ml_dataset_rbns.tsv
"""

import argparse
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.io import load_config, ensure_dir, get_paths, read_tsv, write_tsv


def main():
    parser = argparse.ArgumentParser(description="Build ML dataset from positives and input")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--processed-dir", help="Processed TSV root (for input libraries)")
    parser.add_argument("--results-dir", help="Results dir (tables + output)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n-negative", type=int, default=0, help="Per RBP (0=from config)")
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    processed_dir = Path(args.processed_dir) if args.processed_dir else paths["processed"]
    results_dir = Path(args.results_dir) if args.results_dir else paths["results"]
    tables_dir = results_dir / "tables"
    ensure_dir(results_dir)

    n_negative = args.n_negative or config.get("n_negative_per_rbp", 2000)
    random.seed(args.seed)

    import pandas as pd
    exp_path = paths["metadata"] / "rbns_experiments.tsv"
    target_to_control = {}
    if exp_path.exists():
        exp_df = pd.DataFrame(read_tsv(exp_path))
        for _, row in exp_df.iterrows():
            if str(row.get("is_control", "")).lower() == "true":
                continue
            t = str(row.get("target_name", "")).strip()
            c = str(row.get("control_accession", "")).strip()
            if t and c:
                target_to_control[t] = c

    # Find all *_positives.tsv
    positive_files = list(tables_dir.glob("*_positives.tsv"))
    if not positive_files:
        print(f"No *_positives.tsv in {tables_dir}. Run 04_compute_enrichment.py first.")
        return 1

    dataset_rows = []
    for pos_path in positive_files:
        target_name = pos_path.stem.replace("_positives", "")
        pos_df = pd.DataFrame(read_tsv(pos_path))
        if pos_df.empty or "sequence" not in pos_df.columns:
            continue
        pos_seqs = set(pos_df["sequence"].astype(str).unique())
        # Load input for this target (from control experiment dir or same target dir)
        control_acc = target_to_control.get(target_name)
        if control_acc and (processed_dir / control_acc).exists():
            input_dir = processed_dir / control_acc
        else:
            input_dir = processed_dir / target_name
        if not input_dir.exists():
            input_dir = processed_dir / target_name
        input_files = list(input_dir.glob("*input*0*.tsv")) + list(input_dir.glob("*0nM*.tsv"))
        if not input_files:
            input_files = [f for f in input_dir.glob("*.tsv") if "input" in f.stem.lower() or "0nM" in f.stem]
        if not input_files:
            print(f"  {target_name}: no input library, skip negatives")
            neg_seqs = []
        else:
            in_df = pd.DataFrame(read_tsv(input_files[0]))
            if "sequence" in in_df.columns:
                all_input = in_df["sequence"].astype(str).tolist()
                neg_candidates = [s for s in all_input if s not in pos_seqs]
                # Deduplicate while preserving order
                seen = set()
                neg_candidates = [s for s in neg_candidates if not (s in seen or seen.add(s))]
                random.shuffle(neg_candidates)
                neg_seqs = neg_candidates[:n_negative]
            else:
                neg_seqs = []

        # Carry enrichment metadata from positives TSV into the dataset
        pos_meta = {}
        for col in ["R_max", "n_enriched_concs", "n_concs_measured", "high_confidence"]:
            if col in pos_df.columns:
                pos_meta[col] = pos_df.set_index("sequence")[col].to_dict()

        for seq in pos_df["sequence"].astype(str).tolist():
            row = {
                "target_name": target_name,
                "rna_sequence": seq,
                "binding_label": 1,
                "source": "enriched",
                "R_max": pos_meta.get("R_max", {}).get(seq, ""),
                "n_enriched_concs": pos_meta.get("n_enriched_concs", {}).get(seq, ""),
                "n_concs_measured": pos_meta.get("n_concs_measured", {}).get(seq, ""),
                "high_confidence": pos_meta.get("high_confidence", {}).get(seq, ""),
            }
            dataset_rows.append(row)
        for seq in neg_seqs:
            dataset_rows.append({
                "target_name": target_name,
                "rna_sequence": seq,
                "binding_label": 0,
                "source": "background",
                "R_max": "",
                "n_enriched_concs": "",
                "n_concs_measured": "",
                "high_confidence": "",
            })
        print(f"  {target_name}: {len(pos_seqs)} pos, {len(neg_seqs)} neg")

    random.shuffle(dataset_rows)
    out_path = results_dir / "ml_dataset_rbns.tsv"
    header = ["target_name", "rna_sequence", "binding_label", "source",
              "R_max", "n_enriched_concs", "n_concs_measured", "high_confidence"]
    write_tsv(out_path, dataset_rows, header=header)
    print(f"Wrote {len(dataset_rows)} rows to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
