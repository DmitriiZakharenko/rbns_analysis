#!/usr/bin/env python3
"""
Compute enrichment R per sequence across all available concentrations,
then select top positive sequences.

Two modes (--mode):
  kmer (default):  score(seq) = max R(kmer) over all k-mers in the sequence.
                   R(kmer) = (kmer_count_pulldown / total_reads_pulldown)
                           / (kmer_count_input    / total_reads_input)
                   Matches the ENCODE RBNS computational pipeline (Lambert et al. 2020).
                   More sensitive: aggregates reads across all sequences sharing a k-mer.

  sequence:        R(seq) = f_pulldown(seq) / f_input(seq)  [direct 20-mer frequency]
                   Higher precision, lower sensitivity. Good for strongly binding RBPs only.

Output: results/tables/{target_name}_positives.tsv
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pandas as pd

from utils.io import ensure_dir, get_paths, load_config, read_tsv, write_tsv

RNA_LEN = 20


# ── I/O helpers ───────────────────────────────────────────────────────────────

def load_counts_tsv(path, min_count=1):
    """Load sequence-count TSV. Returns (dict seq->count, total_count_sum)."""
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
    m = re.search(r"_(\d+)nM_", stem)
    return int(m.group(1)) if m else None


# ── Sequence mode ─────────────────────────────────────────────────────────────

def compute_seq_R(pull_counts, pull_total, input_counts, input_total, pseudo):
    """R(seq) = f_pulldown(seq) / f_input(seq) with pseudo-count on input."""
    result = {}
    for seq, count_pull in pull_counts.items():
        count_in = input_counts.get(seq, 0)
        f_in = (count_in + pseudo) / (input_total + pseudo)
        result[seq] = (count_pull / pull_total) / f_in
    return result


# ── K-mer mode ────────────────────────────────────────────────────────────────

def count_kmers(seq_count_dict: dict, k: int):
    """
    Count k-mer occurrences across all sequences, weighted by read count.

    R(kmer) formula uses:
        f_kmer = kmer_count / (total_reads * (RNA_LEN - k + 1))
    The (RNA_LEN - k + 1) factor cancels in R(kmer) = f_pull / f_in, so we
    return raw weighted counts and total reads (not k-mer occurrences).

    Returns (Counter kmer->weighted_count, total_reads).
    """
    kmer_counts: Counter = Counter()
    total_reads = 0
    for seq, count in seq_count_dict.items():
        for i in range(len(seq) - k + 1):
            kmer_counts[seq[i:i + k]] += count
        total_reads += count
    return kmer_counts, total_reads


def compute_kmer_scores(
    pull_counts: dict,
    pull_total: int,
    input_kmer_counts: Counter,
    input_total: int,
    k: int,
    pseudo: float,
) -> dict:
    """
    Score each pulldown sequence by its max k-mer enrichment ratio.

    R(kmer) = (pull_kmer_count / pull_total) / ((in_kmer_count + pseudo) / input_total)
    score(seq) = max R(kmer) over all k-mers in seq.
    """
    pull_kmer_counts, _ = count_kmers(pull_counts, k)

    kmer_R: dict[str, float] = {}
    for kmer, pk in pull_kmer_counts.items():
        in_count = input_kmer_counts.get(kmer, 0)
        f_pull = pk / pull_total
        f_in = (in_count + pseudo) / input_total
        kmer_R[kmer] = f_pull / f_in

    seq_scores: dict[str, float] = {}
    for seq in pull_counts:
        best = 0.0
        for i in range(len(seq) - k + 1):
            r = kmer_R.get(seq[i:i + k], 0.0)
            if r > best:
                best = r
        seq_scores[seq] = best
    return seq_scores


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compute per-sequence enrichment R and select top positives."
    )
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--processed-dir", help="Processed TSV root")
    parser.add_argument("--results-dir", help="Results tables dir")
    parser.add_argument("--top-k", type=int, default=0,
                        help="Max positives per target (0 = from config, default 1000)")
    parser.add_argument("--min-R", type=float, default=0,
                        help="Min enrichment threshold (0 = from config, default 1.5)")
    parser.add_argument("--pseudo-count", type=float, default=0,
                        help="Pseudo-count for input (0 = from config, default 1.0)")
    parser.add_argument("--min-pulldown-count", type=int, default=2,
                        help="Min read count in pulldown TSV (default: 2)")
    parser.add_argument("--min-enriched-concs", type=int, default=1,
                        help="Min concentrations with R>=min_R for high_confidence (default: 1)")
    parser.add_argument("--min-positives", type=int, default=10,
                        help="Skip target if fewer than N positives pass filter (default: 10)")
    parser.add_argument(
        "--mode", choices=["kmer", "sequence"], default="kmer",
        help="Enrichment mode: 'kmer' (default, ENCODE method) or 'sequence' (full 20-mer)",
    )
    parser.add_argument("--kmer-size", type=int, default=5,
                        help="K-mer size for kmer mode (default: 5)")
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
    min_positives = max(1, args.min_positives)
    mode = args.mode
    k = args.kmer_size

    # Load target -> control_accession from experiments metadata
    exp_path = paths["metadata"] / "rbns_experiments.tsv"
    target_to_control: dict[str, str] = {}
    if exp_path.exists():
        exp_df = pd.DataFrame(read_tsv(exp_path))
        for _, row in exp_df.iterrows():
            if str(row.get("is_control", "")).lower() == "true":
                continue
            t = str(row.get("target_name", "")).strip()
            if not t:
                continue
            c = str(row.get("control_accession", "")).strip()
            if t not in target_to_control:
                target_to_control[t] = c
            elif c and not target_to_control[t]:
                target_to_control[t] = c

    all_dirs = [d.name for d in processed_dir.iterdir() if d.is_dir()]
    if target_to_control:
        targets = [t for t in target_to_control if (processed_dir / t).exists()]
        known = set(target_to_control)
        extra = [d for d in all_dirs if not d.startswith("ENCSR") and d not in known]
        if extra:
            print(f"  Note: {len(extra)} dirs not in metadata: {extra}")
            targets += extra
    else:
        targets = [d for d in all_dirs if not d.startswith("ENCSR")]

    if not targets:
        print(f"No target dirs found in {processed_dir}")
        return 1

    print(f"Mode: {mode}" + (f" (k={k})" if mode == "kmer" else ""))
    print(
        f"Processing {len(targets)} targets | min_pulldown_count={min_pulldown_count} | "
        f"min_R={min_R} | top_k={top_k} | min_positives={min_positives}"
    )

    skipped: list[tuple[str, str]] = []
    processed_ok: list[str] = []

    for target_name in sorted(targets):
        target_dir = processed_dir / target_name
        control_acc = target_to_control.get(target_name, "")
        input_dir = (
            processed_dir / control_acc
            if control_acc and (processed_dir / control_acc).exists()
            else target_dir
        )

        # Find input (0 nM) file
        input_candidates = sorted(set(
            [f for f in input_dir.glob("*.tsv") if "input" in f.stem.lower()]
            + [f for f in input_dir.glob("*.tsv") if re.search(r"_0nM_", f.stem)]
        ))
        if not input_candidates:
            reason = f"no input library TSV in {input_dir.name}"
            print(f"  Skip {target_name}: {reason}")
            skipped.append((target_name, reason))
            continue

        input_path = input_candidates[0]
        input_counts, input_total = load_counts_tsv(input_path, min_count=1)
        if input_total == 0:
            reason = f"empty input file {input_path.name}"
            print(f"  Skip {target_name}: {reason}")
            skipped.append((target_name, reason))
            continue

        # Pre-compute input k-mer counts once per target (kmer mode only)
        input_kmer_counts: Counter | None = None
        if mode == "kmer":
            input_kmer_counts, _ = count_kmers(input_counts, k)

        # Collect pulldown TSVs by concentration
        pulldown_by_conc: dict[int, Path] = {}
        for f in target_dir.glob("*.tsv"):
            if "input" in f.stem.lower():
                continue
            conc = parse_conc_from_stem(f.stem)
            if conc is not None and conc > 0:
                pulldown_by_conc[conc] = f

        if not pulldown_by_conc:
            reason = f"no pulldown TSVs in {target_dir.name}"
            print(f"  Skip {target_name}: {reason}")
            skipped.append((target_name, reason))
            continue

        concentrations = sorted(pulldown_by_conc.keys())
        print(
            f"  {target_name}: input={input_path.name} | "
            f"{len(concentrations)} concs: {concentrations}"
        )

        # Compute per-sequence scores at each concentration
        seq_R_by_conc: dict[str, dict[int, float]] = {}
        for conc in concentrations:
            pull_counts_c, pull_total_c = load_counts_tsv(
                pulldown_by_conc[conc], min_count=min_pulldown_count
            )
            if pull_total_c == 0:
                continue

            if mode == "kmer":
                R_map = compute_kmer_scores(
                    pull_counts_c, pull_total_c,
                    input_kmer_counts, input_total, k, pseudo,
                )
            else:
                R_map = compute_seq_R(
                    pull_counts_c, pull_total_c, input_counts, input_total, pseudo
                )

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
            sorted_R = [R_by_conc[c] for c in sorted(R_by_conc)]
            is_monotonic = (
                all(sorted_R[i] <= sorted_R[i + 1] for i in range(len(sorted_R) - 1))
                if len(sorted_R) > 1 else True
            )
            rows.append({
                "sequence": seq,
                "R_max": round(R_max, 6),
                "conc_at_Rmax_nM": conc_at_Rmax,
                "n_enriched_concs": n_enriched,
                "n_concs_measured": len(R_by_conc),
                "is_monotonic": int(is_monotonic),
                "high_confidence": int(n_enriched >= min_enriched_concs and R_max >= min_R),
            })

        # Keep only sequences with R_max >= min_R; never fill with sub-threshold rows
        above = [r for r in rows if r["R_max"] >= min_R]
        above.sort(key=lambda x: (-x["high_confidence"], -x["R_max"]))
        above = above[:top_k]

        if len(above) < min_positives:
            reason = f"too few positives after filter: {len(above)} < {min_positives}"
            skipped.append((target_name, reason))
            continue

        n_hc = sum(r["high_confidence"] for r in above)
        out_path = tables_dir / f"{target_name}_positives.tsv"
        write_tsv(
            out_path, above,
            header=["sequence", "R_max", "conc_at_Rmax_nM", "n_enriched_concs",
                    "n_concs_measured", "is_monotonic", "high_confidence"],
        )
        print(
            f"  {target_name}: {len(above)} positives "
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
