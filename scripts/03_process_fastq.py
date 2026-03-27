#!/usr/bin/env python3
"""
Process FASTQ files: extract sequences, T->U, count unique sequences per file.
Output: data/processed/{target_name}/{library_type}_{concentration_nM}nM.tsv
"""

import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.io import load_config, get_paths, ensure_dir, write_tsv
from utils.sequences import to_rna, is_valid_rna


def read_fastq_seqs(path, expected_len=None):
    """Yield sequence from each FASTQ record. expected_len: int or (min,max) or None."""
    open_fn = gzip.open if str(path).endswith(".gz") else open
    with open_fn(path, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip().upper()
            f.readline()  # +
            f.readline()  # qual
            if not seq:
                continue
            if expected_len is not None:
                if isinstance(expected_len, (list, tuple)):
                    if len(seq) < expected_len[0] or len(seq) > expected_len[1]:
                        continue
                elif len(seq) != expected_len:
                    continue
            yield seq


RNA_LEN = 20  # RBNS random region is always 20 nt; longer reads contain adapter at 3' end


def process_one_fastq(fastq_path, out_path, to_rna_alphabet=True, expected_len=None, min_count=1):
    """Count unique sequences and write TSV. Returns (n_reads, n_unique_written).

    Reads longer than RNA_LEN are trimmed to RNA_LEN (removes 3' adapter).
    expected_len: kept for API compatibility but ignored (always uses RNA_LEN=20).
    min_count: only write sequences with count >= N (default 2 removes singletons).
    """
    counts = Counter()
    n_reads = 0
    try:
        open_fn = gzip.open if str(fastq_path).endswith(".gz") else open
        with open_fn(fastq_path, "rt") as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip().upper()
                f.readline()  # +
                f.readline()  # qual
                if not seq:
                    continue
                # Trim to RNA_LEN if read is longer (adapter at 3' end)
                if len(seq) > RNA_LEN:
                    seq = seq[:RNA_LEN]
                if len(seq) != RNA_LEN:
                    continue
                if to_rna_alphabet:
                    seq = to_rna(seq)
                if is_valid_rna(seq, allowed_lengths=(RNA_LEN,)):
                    n_reads += 1
                    counts[seq] += 1
    except (EOFError, gzip.BadGzipFile) as e:
        print(f"  Warning: {fastq_path.name} appears incomplete or corrupted: {e}")
        if n_reads == 0:
            return 0, 0
    if min_count > 1:
        rows = [{"sequence": s, "count": c} for s, c in counts.most_common() if c >= min_count]
    else:
        rows = [{"sequence": s, "count": c} for s, c in counts.most_common()]
    ensure_dir(Path(out_path).parent)
    write_tsv(out_path, rows, header=["sequence", "count"])
    return n_reads, len(rows)


def main():
    parser = argparse.ArgumentParser(description="Process RBNS FASTQ to sequence-count TSV")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--raw-dir", help="Raw FASTQ root (default: from config)")
    parser.add_argument("--processed-dir", help="Processed output root (default: from config)")
    parser.add_argument("--expected-len", type=int, default=None,
                        help="Expected read length (default: auto-detected per file)")
    parser.add_argument("--limit", type=int, default=0, help="Max FASTQ files to process (0=all)")
    parser.add_argument(
        "--min-count", type=int, default=2,
        help="Only write sequences with count >= N (default: 2, removes singletons which are mostly noise)."
    )
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    raw_dir = Path(args.raw_dir) if args.raw_dir else paths["raw"]
    processed_dir = Path(args.processed_dir) if args.processed_dir else paths["processed"]
    ensure_dir(processed_dir)

    fastq_files = list(raw_dir.rglob("*.fastq.gz")) + list(raw_dir.rglob("*.fastq"))
    if not fastq_files:
        print(f"No FASTQ under {raw_dir}")
        return 1
    if args.limit:
        fastq_files = fastq_files[: args.limit]

    total_reads = 0
    total_unique = 0
    for fp in fastq_files:
        # Target name = parent dir name
        target_name = fp.parent.name
        stem = fp.stem.replace(".fastq", "")
        out_sub = processed_dir / target_name
        ensure_dir(out_sub)
        out_path = out_sub / f"{stem}.tsv"
        n_reads, n_uniq = process_one_fastq(fp, out_path, expected_len=args.expected_len, min_count=args.min_count)
        total_reads += n_reads
        total_unique += n_uniq
        print(f"  {fp.name} -> {out_path.name}  reads={n_reads}  unique={n_uniq}")

    print(f"Total: {len(fastq_files)} files, {total_reads} reads, {total_unique} unique sequences")
    return 0


if __name__ == "__main__":
    sys.exit(main())
