#!/usr/bin/env python3
"""
Compute enrichment R = f_pulldown / f_input per sequence across all available
concentrations, then select top positive sequences.

For each target RBP:
  - R(seq, conc) = f_pulldown(conc) / f_input
  - R_max(seq)   = max R across all concentrations
  - n_enriched   = number of concentrations where R >= min_R (dose-response evidence)

Positives are selected by R_max. Sequences enriched at multiple concentrations
(n_enriched >= 2) are marked as high-confidence.

Output: results/tables/{target_name}_positives.tsv
"""

import argparse
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from utils.io import ensure_dir, get_paths, load_config, read_tsv, write_tsv


def load_counts_tsv(path, min_count=1):
    """Load TSV with sequence/count columns. Returns (dict seq->count, total)."""
    df = pd.DataFrame(read_tsv(path))
    if df.empty or "sequence" not in df.columns or "count" not in df.columns:
        return {}, 0
    df["count"] = pd.to_numeric(df["count"], errors="coerce").fillna(0).astype(int)
    if min_count > 1:
        df = df[df["count"] >= min_count]
        if df.empty:
            return {}, 0
    total = int(df["count"].sum())
    return df.set_index("sequence")["count"].to_dict(), total


def parse_conc_from_stem(stem):
    """Extract integer concentration from filename stem like pulldown_1300nM_ENCFF..."""
    m = re.search(r"_(\d+)nM_", stem)
    return int(m.group(1)) if m else None


def compute_R(pull_counts, pull_total, input_counts, input_total, pseudo):
    """Return dict seq -> R for all sequences in pull_counts."""
    pseudo_f_input = pseudo / (input_total + pseudo)
    result = {}
    for seq, count_pull in pull_counts.items():
        count_in = input_counts.get(seq, 0)
        f_input = (count_in + pseudo) / (input_total + pseudo)
        f_pull = count_pull / pull_total
        R = f_pull / f_input if f_input > 0 else f_pull / pseudo_f_input
        result[seq] = R
    return result


def main():
    parser = argparse.ArgumentParser(description="Compute R(seq) across concentrations and select positives")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--processed-dir", help="Processed TSV root")
    parser.add_argument("--results-dir", help="Results tables dir")
    parser.add_argument("--top-k", type=int, default=0, help="Top K positives per target (0=from config)")
    parser.add_argument("--min-R", type=float, default=0, help="Min R threshold (0=from config)")
    parser.add_argument("--pseudo-count", type=float, default=0, help="Pseudo-count (0=from config)")
    parser.add_argument(
        "--min-pulldown-count", type=int, default=3,
        help="Min count in pulldown to consider a sequence (default: 3, removes noise/singletons)."
    )
    parser.add_argument(
        "--min-enriched-concs", type=int, default=1,
        help="Mark sequence as high_confidence if enriched (R>=min_R) at >= N concentrations (default: 1)."
    )
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    processed_dir = Path(args.processed_dir) if args.processed_dir else paths["processed"]
    results_dir = Path(args.results_dir) if args.results_dir else paths["results"]
    tables_dir = results_dir / "tables"
    ensure_dir(tables_dir)

    top_k = args.top_k or config.get("top_k_positive", 1000)
    min_R = args.min_R or config.get("min_R", 1.5)
    pseudo = args.pseudo_count or config.get("pseudo_count_input", 1.0)
    min_pulldown_count = max(1, args.min_pulldown_count)
    min_enriched_concs = max(1, args.min_enriched_concs)

    # Load target -> control_accession mapping from experiments metadata.
    # ALL pulldown targets are included; control_accession may be empty for some.
    exp_path = paths["metadata"] / "rbns_experiments.tsv"
    target_to_control = {}   # target_name -> control_accession (may be "")
    if exp_path.exists():
        exp_df = pd.DataFrame(read_tsv(exp_path))
        for _, row in exp_df.iterrows():
            if str(row.get("is_control", "")).lower() == "true":
                continue
            t = str(row.get("target_name", "")).strip()
            if not t:
                continue
            c = str(row.get("control_accession", "")).strip()
            # Include all pulldown targets; control_accession may be empty
            if t not in target_to_control:
                target_to_control[t] = c
            elif c and not target_to_control[t]:
                # Prefer non-empty control_accession
                target_to_control[t] = c

    # All subdirs that look like RBP names (not ENCSR control accessions)
    all_dirs = [d.name for d in processed_dir.iterdir() if d.is_dir()]
    if target_to_control:
        # Use metadata-driven list: all pulldown targets that have a processed dir
        targets = [t for t in target_to_control if (processed_dir / t).exists()]
        # Also pick up any non-ENCSR dirs not in metadata (safety net)
        known = set(target_to_control)
        extra = [d for d in all_dirs if not d.startswith("ENCSR") and d not in known]
        if extra:
            print(f"  Note: {len(extra)} dirs not in metadata, processing anyway: {extra}")
            targets += extra
    else:
        targets = [d for d in all_dirs if not d.startswith("ENCSR")]
    if not targets:
        print(f"No target dirs in {processed_dir}")
        return 1

    print(f"Processing {len(targets)} targets, min_pulldown_count={min_pulldown_count}, "
          f"min_R={min_R}, top_k={top_k}, min_enriched_concs={min_enriched_concs}")

    skipped = []
    processed_ok = []

    for target_name in sorted(targets):
        target_dir = processed_dir / target_name

        # Locate input library
        control_acc = target_to_control.get(target_name)
        input_dir = (
            processed_dir / control_acc
            if control_acc and (processed_dir / control_acc).exists()
            else target_dir
        )
        # Use "input" keyword OR exact _0nM_ pattern (underscore-bounded to avoid matching _20nM_, _80nM_)
        input_candidates = sorted(set(
            [f for f in input_dir.glob("*.tsv") if "input" in f.stem.lower()]
            + [f for f in input_dir.glob("*.tsv") if re.search(r"_0nM_", f.stem)]
        ))
        if not input_candidates:
            print(f"  Skip {target_name}: no input library TSV in {input_dir}")
            continue
        input_path = input_candidates[0]
        input_counts, input_total = load_counts_tsv(input_path, min_count=1)
        if input_total == 0:
            print(f"  Skip {target_name}: empty input file {input_path.name}")
            continue

        # Collect ALL pulldown TSVs with their concentrations
        pulldown_by_conc = {}  # conc (int) -> Path
        for f in target_dir.glob("*.tsv"):
            if "input" in f.stem.lower():
                continue
            conc = parse_conc_from_stem(f.stem)
            if conc is not None and conc > 0:
                pulldown_by_conc[conc] = f

        if not pulldown_by_conc:
            reason = f"no pulldown TSVs found in {target_dir}"
            print(f"  Skip {target_name}: {reason}")
            skipped.append((target_name, reason))
            continue

        concentrations = sorted(pulldown_by_conc.keys())
        print(f"  {target_name}: input={input_path.name} | {len(concentrations)} concs: {concentrations}")

        # Compute R for each sequence at each concentration
        # seq -> {conc -> R}
        seq_R_by_conc: dict[str, dict[int, float]] = {}

        for conc in concentrations:
            pull_counts, pull_total = load_counts_tsv(
                pulldown_by_conc[conc], min_count=min_pulldown_count
            )
            if pull_total == 0:
                continue
            R_map = compute_R(pull_counts, pull_total, input_counts, input_total, pseudo)
            for seq, R in R_map.items():
                if seq not in seq_R_by_conc:
                    seq_R_by_conc[seq] = {}
                seq_R_by_conc[seq][conc] = R

        if not seq_R_by_conc:
            reason = f"no sequences after filtering (min_pulldown_count={min_pulldown_count})"
            print(f"  Skip {target_name}: {reason}")
            skipped.append((target_name, reason))
            continue

        # Build result rows
        rows = []
        for seq, R_by_conc in seq_R_by_conc.items():
            R_max = max(R_by_conc.values())
            conc_at_Rmax = max(R_by_conc, key=R_by_conc.get)
            n_enriched = sum(1 for R in R_by_conc.values() if R >= min_R)
            # Monotonic: check R increases with concentration (optional quality signal)
            sorted_R = [R_by_conc[c] for c in sorted(R_by_conc)]
            is_monotonic = all(
                sorted_R[i] <= sorted_R[i + 1] for i in range(len(sorted_R) - 1)
            ) if len(sorted_R) > 1 else True
            rows.append({
                "sequence": seq,
                "R_max": round(R_max, 6),
                "conc_at_Rmax_nM": conc_at_Rmax,
                "n_enriched_concs": n_enriched,
                "n_concs_measured": len(R_by_conc),
                "is_monotonic": int(is_monotonic),
                "high_confidence": int(n_enriched >= min_enriched_concs and R_max >= min_R),
            })

        # Select top_k by R_max; prefer high_confidence
        rows.sort(key=lambda x: (-x["high_confidence"], -x["R_max"]))
        above = [r for r in rows if r["R_max"] >= min_R]
        if len(above) >= top_k:
            rows = above[:top_k]
        else:
            taken = {r["sequence"] for r in above}
            rest = [r for r in rows if r["sequence"] not in taken][: top_k - len(above)]
            rows = (above + rest)[:top_k]
            rows.sort(key=lambda x: (-x["high_confidence"], -x["R_max"]))

        n_hc = sum(r["high_confidence"] for r in rows)
        out_path = tables_dir / f"{target_name}_positives.tsv"
        write_tsv(
            out_path, rows,
            header=["sequence", "R_max", "conc_at_Rmax_nM", "n_enriched_concs",
                    "n_concs_measured", "is_monotonic", "high_confidence"],
        )
        print(
            f"  {target_name}: {len(rows)} positives "
            f"(high_conf={n_hc}, concs={concentrations}) -> {out_path.name}"
        )
        processed_ok.append(target_name)

    print(f"\n=== Step 4 summary ===")
    print(f"  OK:      {len(processed_ok)} targets")
    print(f"  Skipped: {len(skipped)} targets")
    for name, reason in skipped:
        print(f"    - {name}: {reason}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
