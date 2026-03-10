"""ENCODE REST API helpers: search experiments, get file metadata and download URLs."""

import re
import time
from pathlib import Path

import requests

# Base URL for ENCODE (no trailing slash)
ENCODE_BASE = "https://www.encodeproject.org"


def _get(url, params=None, delay=0.2):
    """GET request with optional rate limiting."""
    time.sleep(delay)
    r = requests.get(url, params=params, timeout=60)
    r.raise_for_status()
    # Try JSON first, fallback to text
    try:
        return r.json()
    except ValueError:
        return r.text


def search_experiments(assay_title="RNA Bind-n-Seq", status="released", limit="all", base=ENCODE_BASE):
    """Return list of experiment objects from ENCODE search."""
    url = f"{base}/search/"
    params = {
        "type": "Experiment",
        "assay_title": assay_title,
        "status": status,
        "format": "json",
        "limit": limit,
    }
    data = _get(url, params=params)
    if isinstance(data, dict) and "@graph" in data:
        return data["@graph"]
    return []


def get_experiment(accession, frame="object", base=ENCODE_BASE):
    """Fetch one experiment with optional frame=object for embedded objects."""
    # Accession may be with or without ENCSR prefix
    acc = accession if accession.startswith("ENCSR") else f"ENCSR{accession}"
    url = f"{base}/experiments/{acc}/"
    return _get(url, params={"frame": frame, "format": "json"})


def get_file(file_accession, frame="object", base=ENCODE_BASE):
    """Fetch file metadata. Use frame=object to get embedded replicate with rbns_protein_concentration."""
    acc = file_accession if file_accession.startswith("ENCFF") else f"ENCFF{file_accession}"
    url = f"{base}/files/{acc}/"
    return _get(url, params={"format": "json", "frame": frame})


def get_replicate(replicate_ref, base=ENCODE_BASE):
    """Fetch replicate by @id (e.g. /replicates/uuid/). Returns dict with rbns_protein_concentration, technical_replicate_number."""
    if not replicate_ref or not isinstance(replicate_ref, str):
        return None
    ref = replicate_ref.strip().strip("/")
    # ref is like "replicates/uuid" or "replicates/uuid/" -> take uuid (last part)
    path = ref.split("/")[-1] if "/" in ref else ref
    url = f"{base}/replicates/{path}/"
    try:
        data = _get(url, params={"frame": "object", "format": "json"})
        return data if isinstance(data, dict) else None
    except Exception:
        return None


def concentration_from_replicate_obj(rep_obj, is_control=False):
    """
    Extract concentration from replicate object (from get_replicate).
    Uses rbns_protein_concentration; fallback: RBNS convention by technical_replicate_number.
    """
    if not rep_obj or not isinstance(rep_obj, dict):
        return None
    # ENCODE RBNS replicates have rbns_protein_concentration
    conc = rep_obj.get("rbns_protein_concentration")
    if conc is not None:
        try:
            return int(conc)
        except (TypeError, ValueError):
            pass
    # Fallback: technical_replicate_number -> concentration (RBNS standard)
    if is_control:
        return 0
    tech = rep_obj.get("technical_replicate_number")
    if tech is not None:
        mapping = {1: 5, 2: 20, 3: 80, 4: 320, 5: 1300}
        return mapping.get(int(tech))
    return None


def get_file_download_url(file_accession, base=ENCODE_BASE):
    """Return direct download URL for a file (FASTQ)."""
    acc = file_accession if file_accession.startswith("ENCFF") else file_accession
    return f"{base}/files/{acc}/@@download/{acc}.fastq.gz"


def extract_concentration_from_text(text):
    """Parse 'Protein concentration 80nM' or '0nM' -> 0, 5, 20, 80, 320, 1300 or None."""
    if not text:
        return None
    text = str(text)
    # Match "80nM" or "0 nM" etc
    m = re.search(r"(\d+)\s*nM", text, re.I)
    if m:
        return int(m.group(1))
    return None


def extract_concentration_from_replicate(replicate_obj, base=ENCODE_BASE):
    """Try to get concentration from replicate (library biosample_summary or description). Legacy fallback."""
    if not replicate_obj or isinstance(replicate_obj, str):
        return None
    if not isinstance(replicate_obj, dict):
        return None
    # Prefer rbns_protein_concentration (set by concentration_from_replicate_obj when using get_replicate)
    conc = replicate_obj.get("rbns_protein_concentration")
    if conc is not None:
        try:
            return int(conc)
        except (TypeError, ValueError):
            pass
    # Replicate may have library as @id or embedded
    lib = replicate_obj.get("library")
    if isinstance(lib, dict):
        desc = lib.get("biosample_summary") or lib.get("description") or ""
        return extract_concentration_from_text(desc)
    if isinstance(lib, str):
        lib_id = lib.split("/")[-2] if "/" in lib else lib
        url = f"{base}/libraries/{lib_id}/"
        try:
            lib_obj = _get(url, params={"format": "json"})
            if isinstance(lib_obj, dict):
                desc = lib_obj.get("biosample_summary") or lib_obj.get("description") or ""
                return extract_concentration_from_text(desc)
        except Exception:
            pass
    return None
