# RBNS (RNA Bind-N-Seq): Analysis and ML Dataset Preparation

Step-by-step methodology from raw ENCODE data to a dataset for RBP-RNA binding prediction.

---

## Table of Contents

1. [Project Overview](#1-project-overview)
2. [Data Sources and Scale](#2-data-sources-and-scale)
3. [Fetching Metadata](#3-fetching-metadata)
4. [Downloading FASTQ](#4-downloading-fastq)
5. [Processing FASTQ](#5-processing-fastq)
6. [Enrichment Strategy: Positives and Negatives](#6-enrichment-strategy-positives-and-negatives)
7. [Building the ML Dataset](#7-building-the-ml-dataset)
8. [Quality Validation](#8-quality-validation)
9. [Code and File Structure](#9-code-and-file-structure)
10. [References](#10-references)

---

## 1. Project Overview

### Goal

Build a labelled dataset for ML prediction of RBP-RNA binding:

- **Positives**: sequences enriched in pulldown relative to input at one or more protein concentrations.
- **Negatives**: sequences from the 0 nM input (unselected background library).
- **Features**: RBP name, 20 nt RNA sequence, binding label (0/1), enrichment metrics.

### RBNS vs HTR-SELEX

| Aspect | HTR-SELEX | RBNS |
|--------|-----------|------|
| Source | ENA | ENCODE |
| Structure | Selection cycles, barcodes | 0 nM input + several protein concentrations |
| Concentrations | Not explicit | e.g. 5, 20, 80, 320, 1300 nM (varies per RBP) |
| Enrichment | By cycle frequency | R = f_pulldown / f_input per sequence |
| Read length | 26-40 nt | 20 nt (majority) |

---

## 2. Data Sources and Scale

- **Publication**: Lambert et al., *Nature* 2020 -- [doi:10.1038/s41586-020-2077-3](https://www.nature.com/articles/s41586-020-2077-3)
- **ENCODE collection**: [ENCSR876DCD](https://www.encodeproject.org/publication-data/ENCSR876DCD/)
- **RBNS Computational Pipeline (PDF)**: [ENCODE document](https://www.encodeproject.org/documents/c8b3442a-7e63-4847-af11-c72597bf65b3/@@download/attachment/RBNS_Computational_Pipeline_Aug_2016_update_Dec2018.pdf)

### Actual scale (this pipeline run)

| Resource | Value |
|----------|-------|
| FASTQ files | 601 |
| Raw data size | ~150 GB |
| Processed TSV files | 601 (one per FASTQ) |
| RBPs with enrichment data | 88 |
| RBPs skipped (no signal / no control) | ~23 |

### Disk requirements

| Stage | Approximate size |
|-------|-----------------|
| Metadata only | < 10 MB |
| After step 2 (FASTQ downloaded) | ~150 GB |
| After step 3 (processed TSV, min_count=2) | ~5-10 GB |
| Results only (tables + ML dataset) | < 200 MB |

A VM with at least 200 GB data disk is recommended. See `VM_SETUP_GUIDE.md`.

---

## 3. Fetching Metadata

**Script**: `01_fetch_encode_metadata.py`

**Goal**: retrieve all released RBNS experiments from the ENCODE API and build two tables:

1. `data/metadata/rbns_experiments.tsv` -- experiment list with columns:
   `experiment_accession`, `target_name`, `is_control`, `control_accession`, `description`

2. `data/metadata/rbns_files.tsv` -- file list with columns:
   `experiment_accession`, `target_name`, `file_accession`, `download_url`,
   `concentration_nM`, `library_type` (`input` or `pulldown`), `canonical_name`

**Key design decisions**:

- Pulldown experiments are paired with a control (0 nM input) experiment via ENCODE's `controls` field.
- Many RBPs share a single control experiment. The script resolves this by parsing experiment descriptions when the API's `controlled_by` field is empty.
- The `control_accession` column is filled using a fallback: matching target names in control experiment descriptions (case-insensitive, handles composite names like "IGF2BP1/IMP1").

---

## 4. Downloading FASTQ

**Script**: `02_download_fastq.py`

- Reads `data/metadata/rbns_files.tsv`.
- Downloads each file to `data/raw/{experiment_accession}/{canonical_name}.fastq.gz`.
- Uses `requests` with streaming, automatic retry (exponential backoff), and gzip validation.
- Skips already-downloaded files (resumable).
- Supports `--workers N` for parallel downloads (default 4). Recommended: 4-8 on a VM.
- Logs all outcomes to `data/logs/download.log` and `data/logs/download_summary.json`.

---

## 5. Processing FASTQ

**Script**: `03_process_fastq.py`

For each FASTQ file:

1. Read all reads (supports `.fastq.gz`).
2. Convert DNA to RNA alphabet: `T -> U`.
3. Keep only sequences of the expected length (default 20 nt) with valid ACGU alphabet.
4. Count unique sequences and their read counts.
5. Filter out sequences with count < `min_count` (default **2** -- removes likely sequencing errors).
6. Write to `data/processed/{experiment_accession}/{canonical_name}.tsv`:

```
sequence    count
AAACCCGGGUUUAAACCCGG    1500
...
```

**Note on min_count for input libraries**: Input (0 nM) libraries are processed with `min_count=2` in step 3 (consistent with pulldown files). However, step 4 loads input files with `min_count=1` to preserve the full background representation for frequency normalization.

---

## 6. Enrichment Strategy: Positives and Negatives

### 6.1 Negatives

- **Source**: the 0 nM input library for the matched control experiment.
- **Rationale**: at 0 nM there is no protein -- the library is the unselected background pool.
- **Selection**: up to `n_negative_per_rbp` (default 2000) unique sequences, randomly sampled (seed=42), excluding any sequence already in the positive set (no pos/neg overlap by construction).

### 6.2 Positives -- implemented strategy

We use **all available pulldown concentrations**, not a single fixed concentration. This captures dose-response evidence and is more robust to missing concentrations.

For each sequence at each pulldown concentration:

```
f_pulldown(conc) = count(seq, conc) / total_counts(conc)
f_input          = (count(seq, input) + pseudo) / (total_counts(input) + pseudo)
R(seq, conc)     = f_pulldown(conc) / f_input
```

`pseudo = 1.0` (pseudo-count handles sequences absent from the input).

Derived metrics per sequence (stored in `{target}_positives.tsv`):

| Metric | Description |
|--------|-------------|
| `R_max` | Maximum R across all concentrations |
| `conc_at_Rmax_nM` | Concentration at which R is maximal |
| `n_enriched_concs` | Number of concentrations with R >= min_R |
| `n_concs_measured` | Total pulldown concentrations available |
| `is_monotonic` | 1 if R increases monotonically with concentration |
| `high_confidence` | 1 if R_max >= min_R and n_enriched_concs >= 1 |

**Positive selection**: sort sequences by R_max (high-confidence first), take top `top_k_positive` (default 1000). Remaining slots are filled with highest-R sequences that did not reach `min_R` if needed.

**Parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `min_R` | 1.5 | Minimum enrichment ratio |
| `top_k_positive` | 1000 | Max positives per RBP |
| `--min-pulldown-count` | 2 | Min count in pulldown TSV |
| `--min-enriched-concs` | 1 | Min concentrations with R >= min_R for high_confidence |
| `pseudo_count_input` | 1.0 | Pseudo-count for frequency normalization |

### 6.3 RBPs excluded from the dataset

| Reason | RBPs | Count |
|--------|------|-------|
| No matched input control in ENCODE | APOBEC3C, PTBP3, ZC3H10, TRNAU1AP, TROVE2, ZFP36L2, ... | ~13 |
| No enrichment signal after filtering | CELF1, HNRNPD, HNRNPF, IGF2BP1, MBNL1, MSI1, RBFOX2, RBM47, SRSF2 | 9 |

These are limitations of the ENCODE RBNS collection, not pipeline errors. All skip reasons are logged in `data/logs/step4_final.log`.

---

## 7. Building the ML Dataset

**Script**: `05_build_ml_dataset.py`

For each RBP with a `{target}_positives.tsv` file:

1. Load positive sequences (up to `top_k_positive`).
2. Load the input (0 nM) library; exclude positive sequences from the pool.
3. Random-sample up to `n_negative_per_rbp` negatives (seed=42).
4. Append rows to the dataset with columns:
   `target_name`, `rna_sequence`, `binding_label`, `source`, `R_max`, `n_enriched_concs`, `n_concs_measured`, `high_confidence`.

5. Shuffle all rows (seed=42) and write to `results/ml_dataset_rbns.tsv`.

**Script**: `06_clean_dataset.py`

- Drops rows with missing `target_name`, `rna_sequence`, or `binding_label`.
- Casts `binding_label` to integer (0 or 1).
- Removes duplicates by `(target_name, rna_sequence)`.
- Output: `results/ml_dataset_rbns_clean.tsv`.

---

## 8. Quality Validation

**Script**: `07_validate_dataset.py`

Checks performed on `ml_dataset_rbns_clean.tsv`:

- No missing values in required columns.
- `binding_label` is 0 or 1 (integer).
- Sequence length = 20 nt for all rows.
- Alphabet is ACGU only (no T or N).
- No duplicate `(target_name, rna_sequence)` pairs.
- No `(target_name, rna_sequence)` with both labels 0 and 1.
- Per-RBP positive and negative counts.

Output: `results/dataset_stats_rbns.json`, `results/validation_summary_rbns.tsv`.

### Results from this pipeline run

```json
{
  "n_total": 246777,
  "n_positive": 74135,
  "n_negative": 172642,
  "n_rbp": 88,
  "n_overlap": 0,
  "n_duplicates": 0,
  "seq_length_min": 20,
  "seq_length_max": 20,
  "n_bad_alphabet": 0
}
```

High-confidence positives: 64,688 (87.3% of all positives). Median R_max: 15.8.

---

## 9. Code and File Structure

```
rbns_analysis/
|-- README.md
|-- METHODOLOGY.md          # this document
|-- config.yaml             # parameters: min_R, top_k, paths
|-- config.example.yaml     # documented template
|-- requirements.txt
|-- run_pipeline.sh         # one-command runner (--from N, --only N)
|-- vm_setup.sh             # VM setup (venv, dependencies)
|-- check_vm_disk.sh        # verify /vol/space layout before run
|-- VM_SETUP_GUIDE.md
|-- GIT_AND_RELEASE.md
|-- data/
|   |-- metadata/
|   |   |-- rbns_experiments.tsv
|   |   `-- rbns_files.tsv
|   |-- raw/                # downloaded FASTQ (not committed to git, ~150 GB)
|   |-- processed/          # sequence-count TSV per file (not committed, ~5-10 GB)
|   `-- logs/               # pipeline and download logs (not committed)
|-- results/
|   |-- ml_dataset_rbns.tsv
|   |-- ml_dataset_rbns_clean.tsv
|   |-- dataset_stats_rbns.json
|   |-- validation_summary_rbns.tsv
|   `-- tables/
|       `-- {RBP}_positives.tsv
|-- scripts/
|   |-- 01_fetch_encode_metadata.py
|   |-- 02_download_fastq.py
|   |-- 03_process_fastq.py
|   |-- 04_compute_enrichment.py
|   |-- 05_build_ml_dataset.py
|   |-- 06_clean_dataset.py
|   `-- 07_validate_dataset.py
`-- utils/
    |-- io.py               # config, paths, TSV read/write
    |-- encode_api.py       # ENCODE REST API helpers
    `-- sequences.py        # sequence utilities
```

---

## 10. References

| Resource | URL |
|----------|-----|
| Lambert et al., Nature 2020 | https://www.nature.com/articles/s41586-020-2077-3 |
| ENCODE RBNS collection | https://www.encodeproject.org/publication-data/ENCSR876DCD/ |
| RBNS Computational Pipeline (PDF) | https://www.encodeproject.org/documents/c8b3442a-7e63-4847-af11-c72597bf65b3/@@download/attachment/RBNS_Computational_Pipeline_Aug_2016_update_Dec2018.pdf |
| RBNS Experimental Protocol (PDF) | https://www.encodeproject.org/documents/aa71cabf-aaee-4358-a834-c6ee002938b8/@@download/attachment/RBNSExperimentalProtocol_Feb2016_96well.pdf |
| ENCODE REST API | https://www.encodeproject.org/help/rest-api/ |
