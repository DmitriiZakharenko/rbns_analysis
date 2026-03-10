#!/usr/bin/env python3
"""
Download FASTQ files listed in data/metadata/rbns_files.tsv.
Uses requests with streaming for reliable downloads with retry.
Optional parallel downloads (--workers) to improve throughput.
"""

import argparse
import gzip
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.io import load_config, ensure_dir, get_paths, read_tsv

CHUNK_SIZE = 1024 * 1024  # 1 MB


def download_file(url, path, max_retries=3, timeout=300):
    """
    Download file via requests with streaming.
    Returns (ok, bytes_downloaded, error_message).
    """
    for attempt in range(1, max_retries + 1):
        try:
            with requests.get(url, stream=True, timeout=(30, timeout), allow_redirects=True) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                downloaded = 0
                with open(path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                        f.write(chunk)
                        downloaded += len(chunk)

            if not path.exists() or path.stat().st_size == 0:
                if path.exists():
                    path.unlink()
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return False, 0, "Empty file after download"

            # Validate gzip
            try:
                with gzip.open(path, "rb") as f:
                    f.read(64)
            except (gzip.BadGzipFile, EOFError, OSError) as e:
                path.unlink()
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                    continue
                return False, 0, f"Invalid gzip: {e}"

            return True, path.stat().st_size, None

        except requests.exceptions.RequestException as e:
            if path.exists():
                path.unlink()
            if attempt < max_retries:
                time.sleep(2 ** attempt)
                continue
            return False, 0, str(e)[:200]

    return False, 0, "Max retries exceeded"


def build_file_path(row, raw_dir):
    """Build canonical (subdir, filename, path) from metadata row."""
    acc = row.get("file_accession", "unknown")
    target_name = str(row.get("target_name", "unknown")).strip() or "unknown"
    target_name = target_name.replace("(control)", "control").replace("/", "_")
    if target_name in ("control", "unknown"):
        exp_acc = row.get("experiment_accession", "")
        if exp_acc and str(exp_acc).startswith("ENCSR"):
            target_name = str(exp_acc)

    conc_val = row.get("concentration_nM", "")
    try:
        conc_str = str(int(float(conc_val))) if conc_val and str(conc_val).strip() not in ("", "nan", "None") else ""
    except (ValueError, TypeError):
        conc_str = ""
    if not conc_str and str(row.get("library_type", "")).lower() == "input":
        conc_str = "0"
    if not conc_str:
        conc_str = "unknown"

    lib = row.get("library_type", "pulldown")
    subdir = raw_dir / target_name
    fname = f"{lib}_{conc_str}nM_{acc}.fastq.gz"
    return subdir, fname, subdir / fname


def main():
    parser = argparse.ArgumentParser(description="Download RBNS FASTQ files")
    parser.add_argument("--config", default=ROOT / "config.yaml")
    parser.add_argument("--metadata", help="rbns_files.tsv path")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--limit", type=int, default=0, help="Max files (0 = all)")
    parser.add_argument("--workers", type=int, default=1, help="Parallel downloads (1=sequential, 4–8 often faster over distance)")
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    raw_dir = paths["raw"]
    logs_dir = paths["logs"]
    ensure_dir(raw_dir)
    ensure_dir(logs_dir)

    meta_path = Path(args.metadata) if args.metadata else paths["metadata"] / "rbns_files.tsv"
    if not meta_path.exists():
        print(f"Metadata not found: {meta_path}. Run 01_fetch_encode_metadata.py first.")
        return 1

    rows = read_tsv(meta_path)
    if not rows:
        print("No files in metadata.")
        return 1

    # Filter to FASTQ only
    rows = [r for r in rows if "fastq" in str(r.get("file_type", "")).lower()]
    if not rows:
        print("No FASTQ files in metadata.")
        return 1

    if args.limit:
        rows = rows[:args.limit]

    total = len(rows)
    workers = max(1, min(args.workers, 16))
    print(f"Processing {total} files -> {raw_dir} (workers={workers})")

    log_path = logs_dir / "download.log"
    summary = {"ok": 0, "skip": 0, "fail": 0, "failed_files": []}
    log_lock = threading.Lock()

    # Collect items to download (skip existing)
    to_download = []
    for idx, row in enumerate(rows, 1):
        subdir, fname, path = build_file_path(row, raw_dir)
        if path.exists() and path.stat().st_size > 100:
            summary["skip"] += 1
            if idx % 50 == 0:
                print(f"  [{idx}/{total}] skip {row.get('file_accession', '?')} (exists)")
            continue
        ensure_dir(subdir)
        if args.dry_run:
            print(f"  [{idx}/{total}] would download {row.get('file_accession')} -> {path.relative_to(raw_dir)}")
            continue
        to_download.append((idx, row))

    n_todo = len(to_download)
    if not to_download:
        print("Nothing to download (all files present or dry-run).")
    else:
        print(f"Downloading {n_todo} files...")

    def do_one(item):
        idx, row = item
        acc = row.get("file_accession", "?")
        url = row.get("download_url", "")
        subdir, fname, path = build_file_path(row, raw_dir)
        ok, size, err = download_file(url, path)
        return (idx, acc, fname, url, ok, size, err)

    with open(log_path, "w") as log:
        log.write(f"# Download started at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"# Total files: {total}\n\n")
        log.flush()

        if workers <= 1 or not to_download:
            for idx, row in to_download:
                acc = row.get("file_accession", "?")
                url = row.get("download_url", "")
                subdir, fname, path = build_file_path(row, raw_dir)
                if idx % 10 == 0 or idx == 1:
                    print(f"  [{idx}/{total}] downloading {acc}...")
                ok, size, err = download_file(url, path)
                if ok:
                    summary["ok"] += 1
                    log.write(f"OK\t{fname}\t{size}\n")
                else:
                    summary["fail"] += 1
                    summary["failed_files"].append({"accession": acc, "url": url, "error": err})
                    log.write(f"FAIL\t{fname}\t{err}\n")
                log.flush()
        else:
            done = 0
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {executor.submit(do_one, item): item for item in to_download}
                for fut in as_completed(futures):
                    idx, acc, fname, url, ok, size, err = fut.result()
                    with log_lock:
                        if ok:
                            summary["ok"] += 1
                            log.write(f"OK\t{fname}\t{size}\n")
                        else:
                            summary["fail"] += 1
                            summary["failed_files"].append({"accession": acc, "url": url, "error": err})
                            log.write(f"FAIL\t{fname}\t{err}\n")
                        log.flush()
                    done += 1
                    if done % 10 == 0 or done == 1:
                        print(f"  completed {done}/{n_todo} (last: {acc})")

        log.write(f"\n# Finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        log.write(f"# OK: {summary['ok']}, SKIP: {summary['skip']}, FAIL: {summary['fail']}\n")
        log.flush()

    # Write summary JSON
    summary_path = logs_dir / "download_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\nDone: {summary['ok']} downloaded, {summary['skip']} skipped, {summary['fail']} failed")
    if summary["fail"]:
        print(f"Failed files listed in: {summary_path}")
    print(f"Log: {log_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
