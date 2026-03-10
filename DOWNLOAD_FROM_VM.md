# How to Download Results from the VM to Your Computer

After the pipeline has run on a VM, you only need to copy **results** and optionally **metadata** and **processed** data. You do **not** need to copy raw FASTQ (tens of GB) unless you want to reprocess from scratch.

## What to Copy (priority order)

| Priority | Content | Path on VM | Approx. size | Needed for |
|----------|---------|------------|--------------|------------|
| **1** | Final ML dataset | `results/ml_dataset_rbns_clean.tsv` | &lt; 100 MB | Training models |
| **2** | Dataset stats and validation | `results/dataset_stats_rbns.json`, `results/validation_summary_rbns.tsv` | &lt; 1 MB | Checks |
| **3** | Positive tables per RBP | `results/tables/*_positives.tsv` | &lt; 50 MB | Enrichment details |
| **4** | Processed sequence counts | `data/processed/` | ~1–3 GB | Re-running enrichment without raw FASTQ |
| **5** | Metadata | `data/metadata/*.tsv` | &lt; 5 MB | Reproducibility |
| — | Raw FASTQ | `data/raw/` | ~25–55 GB | Only if you want to re-run from FASTQ |

**Recommended**: Copy at least **1** and **2** (and optionally **3**). Omit `data/raw/` unless you need it.

---

## Option A: `scp` (from your laptop)

Replace `USER`, `VM_IP`, and `REMOTE_PATH` with your VM user, IP/hostname, and project path.

```bash
# Create a folder on your computer
mkdir -p ~/rbns_results

# Copy only results (dataset + stats + tables)
scp -r USER@VM_IP:REMOTE_PATH/rbns_analysis/results ~/rbns_results/

# Optional: copy metadata
scp -r USER@VM_IP:REMOTE_PATH/rbns_analysis/data/metadata ~/rbns_results/
```

Example:

```bash
scp -r ubuntu@12.34.56.78:/home/ubuntu/rbns_analysis/results ~/rbns_results/
```

---

## Option B: `rsync` (efficient, resumable)

Good for large `results/` or if you add `data/processed/` later.

```bash
mkdir -p ~/rbns_results
rsync -avz --progress USER@VM_IP:REMOTE_PATH/rbns_analysis/results/ ~/rbns_results/results/
# Optional: metadata
rsync -avz USER@VM_IP:REMOTE_PATH/rbns_analysis/data/metadata/ ~/rbns_results/metadata/
```

If the connection drops, run the same `rsync` again; it will resume.

---

## Option C: Create a tarball on the VM, then download

On the **VM**:

```bash
cd /path/to/rbns_analysis
tar -czvf rbns_results.tar.gz results data/metadata
# Optional: include processed data (larger)
# tar -czvf rbns_results_with_processed.tar.gz results data/metadata data/processed
```

Then from your **laptop**:

```bash
scp USER@VM_IP:REMOTE_PATH/rbns_analysis/rbns_results.tar.gz ~/Downloads/
cd ~/Downloads && tar -xzvf rbns_results.tar.gz
```

---

## Option D: Cloud storage (S3, GCS, etc.)

If the VM is in the same cloud as your bucket:

**On the VM** (example: AWS S3):

```bash
aws s3 cp --recursive results/ s3://YOUR_BUCKET/rbns_analysis/results/
aws s3 cp data/metadata/ s3://YOUR_BUCKET/rbns_analysis/metadata/ --recursive
```

Then download from the bucket to your computer (browser, CLI, or sync tool).

---

## Summary

- **Minimum**: copy `results/` (dataset + stats + tables) with `scp` or `rsync`.
- **Recommended**: `results/` + `data/metadata/`.
- **Full backup** (no raw FASTQ): add `data/processed/`; use `rsync` or a tarball.
- You do **not** need to download `data/raw/` (FASTQ) unless you plan to re-run the pipeline from scratch.
