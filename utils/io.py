"""I/O helpers: config, paths, TSV read/write."""

import json
from pathlib import Path

import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def load_config(config_path=None):
    """Load YAML config; default: config.yaml in project root."""
    if config_path is None:
        config_path = PROJECT_ROOT / "config.yaml"
    path = Path(config_path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with open(path) as f:
        data = yaml.safe_load(f)
    return data or {}


def get_paths(config=None):
    """Return paths as Path objects relative to project root. Single source of truth for all scripts."""
    if config is None:
        config = load_config()
    root = PROJECT_ROOT
    raw = config.get("paths", {})
    return {
        key: (root / value) if isinstance(value, str) else Path(value)
        for key, value in raw.items()
    }


def ensure_dir(path):
    """Create directory and parents if needed."""
    Path(path).mkdir(parents=True, exist_ok=True)


def read_tsv(path):
    """Read TSV into list of dicts. Use pandas in callers if needed: pd.DataFrame(read_tsv(path))."""
    path = Path(path)
    if not path.exists():
        return []
    with open(path) as f:
        lines = f.readlines()
    if not lines:
        return []
    header = lines[0].strip().split("\t")
    return [dict(zip(header, line.strip().split("\t"))) for line in lines[1:] if line.strip()]


def write_tsv(path, rows, header=None):
    """Write list of dicts to TSV. If header is None, use keys of first row."""
    ensure_dir(Path(path).parent)
    if not rows and header:
        with open(path, "w") as f:
            f.write("\t".join(header) + "\n")
        return
    if not rows:
        return
    if header is None:
        header = list(rows[0].keys())
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(k, "")) for k in header) + "\n")


def write_json(path, obj):
    """Write JSON to file."""
    ensure_dir(Path(path).parent)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
