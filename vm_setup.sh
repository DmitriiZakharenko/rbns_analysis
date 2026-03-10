#!/usr/bin/env bash
set -euo pipefail

# VM setup: run ONCE after uploading the project to the VM.
# Run from project root, e.g.: cd /vol/space/rbns_analysis && bash vm_setup.sh

echo "=== RBNS VM Setup ==="

sudo apt-get update -qq
sudo apt-get install -y -qq python3 python3-pip python3-venv git tmux htop

# Use directory containing this script as project root (so /vol/space/rbns_analysis works)
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"
echo "Project root: $ROOT"

python3 -m venv .venv
source .venv/bin/activate

pip install --upgrade pip
pip install -r requirements.txt

echo ""
echo "=== Verifying ==="
python3 -c "import requests, yaml, pandas; print('OK: requests', requests.__version__, 'pandas', pandas.__version__)"

mkdir -p data/raw data/processed data/logs results

echo ""
echo "=== Metadata check ==="
if [[ -f data/metadata/rbns_files.tsv ]]; then
    NFILES=$(tail -n +2 data/metadata/rbns_files.tsv | wc -l)
    echo "rbns_files.tsv: $NFILES files ready for download"
else
    echo "WARNING: metadata not found — step 1 will fetch it from ENCODE API."
fi

echo ""
echo "=== Setup complete ==="
echo "Next steps:"
echo "  source .venv/bin/activate"
echo "  tmux new -s rbns"
echo "  bash run_pipeline.sh --from 2   # metadata already done, start from download"
