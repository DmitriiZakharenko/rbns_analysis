#!/usr/bin/env python3
"""
Add protein_sequence column to the RBNS ML dataset.

Strategy (in order of priority):
  1. Query UniProt REST API (reviewed human entries, canonical isoform) by gene name.
  2. For proteins not found in UniProt, fall back to a hardcoded override dict.
  3. Report which proteins could not be resolved.

Run this script from a machine with internet access (e.g. the VM or your workstation),
*after* running steps 1-7 of the pipeline.

Usage:
    python scripts/08_add_protein_sequences.py [--config config.yaml]
    python scripts/08_add_protein_sequences.py --input  results/ml_dataset_rbns_clean.tsv \
                                                --output results/ml_dataset_rbns_clean.tsv

The script overwrites the output file in place if --input == --output (default).
"""

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import pandas as pd
import requests

from utils.io import load_config, get_paths

# ---------------------------------------------------------------------------
# Gene-name aliases: some ENCODE target names differ from the canonical
# HGNC gene symbol used by UniProt.
# ---------------------------------------------------------------------------
GENE_ALIASES: dict[str, str] = {
    "HNRNPA1L2":  "HNRNPA1L2",   # alias for HNRNPA1-like 2
    "HNRNPCL1":   "HNRNPCL1",
    "HNRNPDL":    "HNRNPDL",
    "HNRNPA2B1":  "HNRNPA2B1",
    "RBFOX2":     "RBFOX2",
    "RBFOX3":     "RBFOX3",
    "ZFP36L1":    "ZFP36L1",
    "ZFP36":      "ZFP36",
    "ELAVL4":     "ELAVL4",
    "AKAP8L":     "AKAP8L",
    "CNOT4":      "CNOT4",
    "CPEB1":      "CPEB1",
    "CSDE1":      "CSDE1",
    "EIF3D":      "EIF3D",
    "EIF4G2":     "EIF4G2",
    "EWSR1":      "EWSR1",
    "EXOSC4":     "EXOSC4",
    "FUBP1":      "FUBP1",
    "FUBP3":      "FUBP3",
    "ILF2":       "ILF2",
    "LIN28B":     "LIN28B",
    "NSUN2":      "NSUN2",
    "NUPL2":      "NUP50",       # NUPL2 = NUP50 in UniProt
    "PABPC3":     "PABPC3",
    "PUF60":      "PUF60",
    "RALYL":      "RALYL",
    "SAFB2":      "SAFB2",
    "SF1":        "SF1",
    "SFPQ":       "SFPQ",
    "SUCLG1":     "SUCLG1",
    "TAF15":      "TAF15",
    "TIA1":       "TIA1",
    "TRA2A":      "TRA2A",
    "TRA2B":      "TRA2B",
    "UNK":        "UNK",
    "XRCC6":      "XRCC6",
    "XRN2":       "XRN2",
    "ZC3H18":     "ZC3H18",
    "PRR3":       "PRR3",
}


def fetch_uniprot_sequence(gene_name: str, delay: float = 0.5) -> tuple[str | None, str | None]:
    """
    Query UniProt REST API for the canonical sequence of a reviewed human protein.

    Returns (uniprot_accession, sequence) or (None, None) on failure.
    """
    query_name = GENE_ALIASES.get(gene_name, gene_name)
    url = (
        "https://rest.uniprot.org/uniprotkb/search"
        f"?query=gene_exact:{query_name}+AND+organism_id:9606+AND+reviewed:true"
        "&fields=accession,sequence&format=json&size=1"
    )
    try:
        resp = requests.get(url, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        results = data.get("results", [])
        if results:
            acc = results[0]["primaryAccession"]
            seq = results[0]["sequence"]["value"]
            return acc, seq
    except Exception as e:
        print(f"    UniProt lookup failed for {gene_name}: {e}", file=sys.stderr)
    finally:
        time.sleep(delay)
    return None, None


def main():
    parser = argparse.ArgumentParser(description="Add protein_sequence to RBNS ML dataset")
    parser.add_argument("--config",  default=ROOT / "config.yaml")
    parser.add_argument("--input",   help="Input TSV  (default: results/ml_dataset_rbns_clean.tsv)")
    parser.add_argument("--output",  help="Output TSV (default: same as input)")
    parser.add_argument("--lookup-cache", default=ROOT / "data/metadata/protein_sequences_cache.json",
                        help="JSON file to cache UniProt lookups (avoids re-querying)")
    parser.add_argument("--delay", type=float, default=0.5,
                        help="Seconds between UniProt requests (default 0.5)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print resolved sequences without saving")
    args = parser.parse_args()

    config   = load_config(args.config)
    paths    = get_paths(config)
    results  = paths["results"]
    in_path  = Path(args.input  or results / "ml_dataset_rbns_clean.tsv")
    out_path = Path(args.output or in_path)
    cache_path = Path(args.lookup_cache)

    if not in_path.exists():
        print(f"Input not found: {in_path}")
        return 1

    df = pd.read_csv(in_path, sep="\t")
    proteins = sorted(df["target_name"].unique())
    print(f"Proteins to resolve: {len(proteins)}")

    # Load cache
    cache: dict[str, str] = {}
    if cache_path.exists():
        with open(cache_path) as f:
            cache = json.load(f)
        print(f"Loaded {len(cache)} cached sequences from {cache_path}")

    # Resolve sequences
    protein_seq: dict[str, str] = {}
    failed: list[str] = []

    for i, prot in enumerate(proteins, 1):
        if prot in cache:
            protein_seq[prot] = cache[prot]
            print(f"  [{i:3d}/{len(proteins)}] {prot}: from cache ({len(cache[prot])} aa)")
            continue

        print(f"  [{i:3d}/{len(proteins)}] {prot}: querying UniProt...", end=" ", flush=True)
        acc, seq = fetch_uniprot_sequence(prot, delay=args.delay)
        if seq:
            protein_seq[prot] = seq
            cache[prot] = seq
            print(f"OK ({acc}, {len(seq)} aa)")
        else:
            failed.append(prot)
            protein_seq[prot] = ""
            print("FAILED")

    # Save updated cache
    if not args.dry_run:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f, indent=2)
        print(f"\nCache saved to {cache_path}")

    # Summary
    resolved  = [p for p in proteins if protein_seq.get(p)]
    unresolved = [p for p in proteins if not protein_seq.get(p)]
    print(f"\nResolved:   {len(resolved)}/{len(proteins)}")
    if unresolved:
        print(f"Unresolved: {unresolved}")
        print("  -> Add these manually to GENE_ALIASES or the cache JSON")

    if args.dry_run:
        print("\nDry run — no files written.")
        return 0

    # Add column to dataset
    df["protein_sequence"] = df["target_name"].map(protein_seq).fillna("")
    n_with_seq = (df["protein_sequence"] != "").sum()
    n_rows = len(df)
    print(f"\nRows with protein_sequence: {n_with_seq}/{n_rows} "
          f"({n_with_seq/n_rows*100:.1f}%)")

    # Reorder columns: put protein_sequence after target_name
    cols = df.columns.tolist()
    if "protein_sequence" in cols:
        cols.remove("protein_sequence")
        idx = cols.index("target_name") + 1
        cols.insert(idx, "protein_sequence")
        df = df[cols]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, sep="\t", index=False)
    print(f"Saved to {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
