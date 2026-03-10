#!/usr/bin/env bash
set -euo pipefail

#
# RBNS Pipeline Runner
# Runs all steps sequentially: metadata -> download -> process -> enrichment -> ML dataset -> clean -> validate
# Usage:
#   bash run_pipeline.sh              # run all steps
#   bash run_pipeline.sh --from 2     # start from step 2 (download)
#   bash run_pipeline.sh --only 2     # run only step 2
#

cd "$(dirname "$0")"

LOGDIR="data/logs"
mkdir -p "$LOGDIR"
MASTER_LOG="$LOGDIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

log() { echo "[$(date '+%H:%M:%S')] $*" | tee -a "$MASTER_LOG"; }
die() { log "FATAL: $*"; exit 1; }

FROM_STEP=1
ONLY_STEP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        --from)  FROM_STEP="$2"; shift 2 ;;
        --only)  ONLY_STEP="$2"; shift 2 ;;
        *)       die "Unknown arg: $1" ;;
    esac
done

should_run() {
    local step=$1
    if [[ $ONLY_STEP -gt 0 ]]; then
        [[ $step -eq $ONLY_STEP ]]
    else
        [[ $step -ge $FROM_STEP ]]
    fi
}

run_step() {
    local num=$1; shift
    local desc=$1; shift
    if ! should_run "$num"; then
        log "SKIP step $num: $desc"
        return 0
    fi
    log "========================================"
    log "START step $num: $desc"
    log "========================================"
    local t0=$SECONDS
    "$@" 2>&1 | tee -a "$MASTER_LOG"
    local rc=${PIPESTATUS[0]}
    local elapsed=$(( SECONDS - t0 ))
    if [[ $rc -ne 0 ]]; then
        die "step $num FAILED (exit=$rc, ${elapsed}s)"
    fi
    log "DONE  step $num: ${elapsed}s"
    echo ""
}

log "Pipeline started. Log: $MASTER_LOG"
log "Python: $(python3 --version 2>&1)"
log "Disk free: $(df -h . | tail -1 | awk '{print $4}')"

run_step 1 "Fetch ENCODE metadata" \
    python3 scripts/01_fetch_encode_metadata.py

run_step 2 "Download FASTQ files" \
    python3 scripts/02_download_fastq.py

run_step 3 "Process FASTQ -> sequence counts" \
    python3 scripts/03_process_fastq.py --min-count 2

run_step 4 "Compute enrichment R(seq) across all concentrations" \
    python3 scripts/04_compute_enrichment.py --min-pulldown-count 2 --min-enriched-concs 1

run_step 5 "Build ML dataset" \
    python3 scripts/05_build_ml_dataset.py

run_step 6 "Clean ML dataset" \
    python3 scripts/06_clean_dataset.py

run_step 7 "Validate ML dataset" \
    python3 scripts/07_validate_dataset.py

log "========================================"
log "ALL STEPS COMPLETE"
log "========================================"

if [[ -f "data/logs/download_summary.json" ]]; then
    log "Download summary:"
    python3 -m json.tool data/logs/download_summary.json 2>/dev/null | tee -a "$MASTER_LOG" || true
fi
if [[ -f "results/dataset_stats_rbns.json" ]]; then
    log "Dataset stats:"
    python3 -m json.tool results/dataset_stats_rbns.json 2>/dev/null | tee -a "$MASTER_LOG" || true
fi

log "Pipeline log: $MASTER_LOG"
