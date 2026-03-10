#!/usr/bin/env python3
"""
Fetch RBNS experiment and file metadata from ENCODE.
Output: data/metadata/rbns_experiments.tsv, data/metadata/rbns_files.tsv
"""

import argparse
import sys
from pathlib import Path

# Project root
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from utils.encode_api import (
    ENCODE_BASE,
    concentration_from_replicate_obj,
    get_experiment,
    get_file,
    get_file_download_url,
    get_replicate,
    search_experiments,
)
from utils.io import ensure_dir, get_paths, load_config, write_tsv


def main():
    parser = argparse.ArgumentParser(description="Fetch RBNS metadata from ENCODE")
    parser.add_argument("--config", default=ROOT / "config.yaml", help="Config YAML path")
    parser.add_argument("--out-dir", help="Output dir for TSV (default: from config)")
    parser.add_argument("--limit", type=int, default=0, help="Max experiments to process (0 = all)")
    args = parser.parse_args()

    config = load_config(args.config)
    paths = get_paths(config)
    out_dir = Path(args.out_dir) if args.out_dir else paths["metadata"]
    ensure_dir(out_dir)

    print("Searching ENCODE for RBNS experiments...")
    experiments = search_experiments()
    if args.limit:
        experiments = experiments[: args.limit]
        print(f"Processing first {len(experiments)} experiments (--limit={args.limit})")
    else:
        print(f"Found {len(experiments)} experiments")

    # Separate pulldown vs control by control_type or absence of target
    experiments_list = []
    files_list = []
    total = len(experiments)

    for idx, exp in enumerate(experiments, 1):
        if idx % 10 == 0 or idx == 1:
            print(f"Processing experiment {idx}/{total}...")
        accession = exp.get("accession") or exp.get("@id", "").split("/")[-2]
        if not accession:
            continue
        # Control experiments have control_type; pulldown have target
        control_type = exp.get("control_type")
        target = exp.get("target")
        is_control = control_type is not None or target is None
        target_name = ""
        if target:
            if isinstance(target, dict):
                target_name = target.get("label") or target.get("name") or ""
            else:
                target_name = str(target).split("/")[-2] if "/" in str(target) else str(target)
        # Fetch full experiment (files and replicates are list of @id strings)
        # Also need full object to check controlled_by for pulldown->control links
        try:
            full = get_experiment(accession, frame="object")
            if not isinstance(full, dict):
                print(f"Warning: {accession} returned non-dict, skipping")
                continue
        except Exception as e:
            print(f"Warning: could not fetch {accession}: {e}")
            continue

        # Control experiment description may list RBPs (e.g. "input control of NUPL2 and PRR3")
        controls = full.get("controls", [])
        control_accession = ""
        if controls:
            c = controls[0]
            control_accession = c.get("accession") or (c.get("@id", "").split("/")[-2] if isinstance(c, dict) else str(c).split("/")[-2])
        # Link pulldown to its control (check in full object)
        if not is_control and not control_accession and full.get("controlled_by"):
            controlled_by = full.get("controlled_by", [])
            if controlled_by:
                c = controlled_by[0]
                control_accession = c.get("accession") or (c.get("@id", "").split("/")[-2] if isinstance(c, dict) else str(c).split("/")[-2])

        experiments_list.append({
            "experiment_accession": accession,
            "target_name": target_name,
            "control_accession": control_accession,
            "is_control": "true" if is_control else "false",
            "description": (full.get("description") or exp.get("description") or "")[:200],
        })
        file_refs = full.get("files", [])
        replicate_refs = full.get("replicates", [])
        # Build replicate @id -> concentration by fetching each replicate (has rbns_protein_concentration)
        rep_concentration = {}
        for rep_ref in replicate_refs:
            if not isinstance(rep_ref, str) or "/replicates/" not in rep_ref:
                continue
            rep_obj = get_replicate(rep_ref)
            conc = concentration_from_replicate_obj(rep_obj, is_control=is_control)
            if conc is not None:
                # Normalize key for lookup: strip and ensure leading slash, no trailing slash
                key = rep_ref.strip().rstrip("/")
                if not key.startswith("/"):
                    key = "/" + key
                rep_concentration[key] = conc

        # For each file, fetch file metadata to get replicate link, then look up concentration
        for file_ref in file_refs:
            if not isinstance(file_ref, str) or "/files/" not in file_ref:
                continue
            file_acc = file_ref.strip("/").split("/")[-1]
            if not file_acc or "ENCFF" not in file_acc:
                continue
            try:
                file_obj = get_file(file_acc)
            except Exception:
                continue
            if not isinstance(file_obj, dict):
                continue
            file_type = (file_obj.get("file_type") or "").lower()
            if "fastq" not in file_type and file_type != "fastq":
                continue
            rep_ref = file_obj.get("replicate")
            conc = None
            # Replicate may be embedded (dict) when file was fetched with frame=object
            if isinstance(rep_ref, dict):
                conc = concentration_from_replicate_obj(rep_ref, is_control=is_control)
            elif isinstance(rep_ref, str):
                key = rep_ref.strip().rstrip("/")
                if not key.startswith("/"):
                    key = "/" + key
                conc = rep_concentration.get(key)
                if conc is None and "/replicates/" in rep_ref:
                    rep_obj = get_replicate(rep_ref)
                    conc = concentration_from_replicate_obj(rep_obj, is_control=is_control)
            if conc is None and is_control:
                conc = 0
            download_url = get_file_download_url(file_acc)
            library_type = "input" if (conc == 0 or (is_control and conc is not None)) else "pulldown"
            files_list.append({
                "file_accession": file_acc,
                "experiment_accession": accession,
                "target_name": target_name or "(control)",
                "concentration_nM": conc if conc is not None else "",
                "library_type": library_type,
                "download_url": download_url,
                "file_type": "fastq",
            })

    # Deduplicate files by file_accession (same file can appear in multiple experiments)
    seen = set()
    unique_files = []
    for row in files_list:
        if row["file_accession"] in seen:
            continue
        seen.add(row["file_accession"])
        unique_files.append(row)

    # Fill control_accession from control experiment descriptions (API often leaves controlled_by empty)
    import re
    # lowercase target -> (accession, is_input_control)
    target_lc_to_control = {}
    for row in experiments_list:
        if row.get("is_control") != "true":
            continue
        acc = row.get("experiment_accession", "")
        desc = (row.get("description") or "").strip()
        if not acc or "control" not in desc.lower():
            continue
        m = re.search(r"\b(?:control|controlled)\s+(?:experiment\s+)?(?:of|for)\s+([^.]+)", desc, re.I)
        if not m:
            continue
        is_input = "input" in desc.lower()
        names_str = m.group(1).strip()
        # Split on " and ", " or ", commas
        parts = re.split(r"\s+(?:and|or)\s+|,\s*", names_str, flags=re.I)
        for part in parts:
            # Handle composite names like "IGF2BP1/IMP1" -> ["IGF2BP1", "IMP1"]
            subnames = [s.strip() for s in part.split("/") if s.strip()]
            for t in subnames:
                if not t or len(t) > 50:
                    continue
                t_lc = t.lower()
                prev = target_lc_to_control.get(t_lc)
                if prev is None or (is_input and not prev[1]):
                    target_lc_to_control[t_lc] = (acc, is_input)

    for row in experiments_list:
        if row.get("is_control") == "false" and not (row.get("control_accession") or "").strip():
            t = (row.get("target_name") or "").strip()
            match = target_lc_to_control.get(t.lower())
            if match:
                row["control_accession"] = match[0]

    write_tsv(out_dir / "rbns_experiments.tsv", experiments_list)
    write_tsv(out_dir / "rbns_files.tsv", unique_files, header=list(unique_files[0].keys()) if unique_files else None)
    print(f"Wrote {len(experiments_list)} experiments, {len(unique_files)} unique FASTQ files to {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
