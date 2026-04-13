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
| RBPs included (with control + k-mer R_max ≥ 1.5) | 96 |
| RBPs excluded — no input control in ENCODE | 11 |
| RBPs excluded — no detectable k-mer enrichment | 3 |

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

### 6.2 Positives — k-mer enrichment strategy

We use **k-mer enrichment** (k=5, matching the ENCODE RBNS computational pipeline) rather than full per-sequence frequency ratios. This aggregates signal across all reads sharing a k-mer, giving ~1000× more statistical power than 20-mer analysis for a 150M-read library.

For each k-mer at each pulldown concentration:

```
R(kmer, conc) = (kmer_count_pulldown(conc) / total_reads_pulldown(conc))
              / (kmer_count_input          / total_reads_input)
```

For each 20-mer sequence:

```
score(seq, conc) = max R(kmer, conc)  over all k-mers within seq
R_max(seq)       = max score(seq, conc) over all concentrations
```

`pseudo = 1.0` is added to `kmer_count_input` to handle k-mers absent from input.

This approach detects binding signal for **96 RBPs** (vs 53 with full-20mer method), capturing RBPs whose binding preference manifests at the k-mer level but is distributed across too many distinct 20-mers to show strong per-sequence enrichment.

Derived metrics per sequence (stored in `{target}_positives.tsv`):

| Metric | Description |
|--------|-------------|
| `R_max` | Maximum R across all concentrations |
| `conc_at_Rmax_nM` | Concentration at which R is maximal |
| `n_enriched_concs` | Number of concentrations with R >= min_R |
| `n_concs_measured` | Total pulldown concentrations available |
| `is_monotonic` | 1 if R increases monotonically with concentration |
| `high_confidence` | 1 if R_max >= min_R and n_enriched_concs >= 1 |

**Positive selection**: keep only sequences with `R_max >= min_R`, sort by high-confidence then R_max, take up to `top_k_positive` (default 1000). Targets with fewer than `min_positives` (default 10) after filtering are skipped.

**Known limitation — R_max ceiling effect**: for 83 of 96 RBPs, all 1000 positive sequences share the *same* R_max value. This happens because R_max = max R over all k-mers in the sequence, and when there is one dominant binding k-mer (e.g. UUUUU for HNRNPC) virtually all top sequences contain it and receive that k-mer's enrichment score. As a result, R_max is effectively a binary indicator of "contains the binding k-mer" and does **not** discriminate between positive sequences within a single RBP. It remains a valid confidence threshold for defining the positive set, but should not be used as a continuous feature in ML models trained per-RBP.

**Parameters**:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `min_R` | 1.5 | Minimum enrichment ratio |
| `top_k_positive` | 1000 | Max positives per RBP |
| `--min-pulldown-count` | 2 | Min count in pulldown TSV |
| `--min-enriched-concs` | 1 | Min concentrations with R >= min_R for high_confidence |
| `pseudo_count_input` | 1.0 | Pseudo-count for frequency normalization |

### 6.3 RBPs excluded from the dataset

| Reason | Examples | Count |
|--------|----------|-------|
| No matched 0 nM input control in ENCODE | APOBEC3C, EIF4H, KHSRP, PABPN1L, PPP1R10, RBM20, SF3B6, TRNAU1AP, TROVE2, ZC3H10, ZFP36L2 | 11 |
| Matched control, but no k-mer with R ≥ 1.5 at any concentration | SRSF10, SYNCRIP, PTBP3 | 3 |

**On the first group:** these RBPs have only pulldown FASTQ files in ENCODE, with no 0 nM input library. Without a background reference the R value is unreliable.

**On the second group:** even the k-mer enrichment signal is too weak to produce 10+ positives with R_max ≥ 1.5. These RBPs may bind RNA in a structure-dependent or context-dependent manner that is not captured by 20-mer random library assays.

### 6.4 Note on read lengths

Some ENCODE RBNS experiments use reads of 40 nt (e.g. CELF1, RBFOX2, MBNL1, SRSF2 and others from an earlier experimental batch). The RBNS random RNA region is always **20 nt**; the additional length corresponds to adapter sequence at the 3' end.

The pipeline automatically trims reads longer than 20 nt to the first 20 nucleotides (the random region). This is handled in `03_process_fastq.py` via the constant `RNA_LEN = 20`.

---

## 7. Building the ML Dataset

**Script**: `05_build_ml_dataset.py`

For each RBP with a `{target}_positives.tsv` file:

1. Load positive sequences (up to `top_k_positive`).
2. Load the input (0 nM) library; exclude positive sequences from the pool.
3. Random-sample up to `n_negative_per_rbp` negatives (seed=42, `random.shuffle` — reproducible).
4. Append rows to the dataset with columns:
   `target_name`, `rna_sequence`, `binding_label`, `source`, `R_max`, `n_enriched_concs`, `n_concs_measured`, `high_confidence`.

5. Shuffle all rows (seed=42) and write to `results/ml_dataset_rbns.tsv`.

**Note on negative counts**: six RBPs end up with fewer than 2000 negatives (NSUN2: 600, SUCLG1: 720, PCBP4: 1660, PRR3/NUPL2: 1843, IGF2BP3: 1976). This is because their shared 0 nM control experiments have smaller processed input pools (shallow sequencing depth or low unique-sequence diversity in those libraries). The pipeline takes all available unique sequences up to the target count.

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
  "n_total": 284642,
  "n_positive": 96000,
  "n_negative": 188642,
  "n_rbp": 96,
  "n_overlap": 0,
  "n_duplicates": 0,
  "seq_length_min": 20,
  "seq_length_max": 20,
  "seq_length_mean": 20.0,
  "n_bad_alphabet": 0
}
```

All 96,000 positive sequences (1000 per RBP) have k-mer R_max ≥ 1.5. All 96 RBPs have a verified matched 0 nM input control. No quality issues detected.

---

## 9. Protein Sequences (optional step 8)

The final dataset (`ml_dataset_rbns_clean.tsv`) currently contains **no `protein_sequence` column** because the ENCODE API only provides target gene names, not amino-acid sequences. Script `08_add_protein_sequences.py` resolves this:

- Queries the **UniProt REST API** (`reviewed:true`, organism 9606) by gene name for each of the 96 RBPs.
- Caches results in `data/metadata/protein_sequences_cache.json` (re-run safe).
- Inserts a `protein_sequence` column after `target_name` in the output TSV.
- Run from a machine with internet access (the VM or your workstation):

```bash
python scripts/08_add_protein_sequences.py
```

**Important caveat**: UniProt returns the *canonical full-length* human protein sequence. The ENCODE RBNS experiments were performed with recombinant proteins of unspecified construct length (some may be RNA-binding-domain-only constructs). Using the full-length sequence is the standard and reproducible choice; keep this in mind when comparing sequence features across datasets.

---

## 10. Code and File Structure

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

## 11. References

| Resource | URL |
|----------|-----|
| Lambert et al., Nature 2020 | https://www.nature.com/articles/s41586-020-2077-3 |
| ENCODE RBNS collection | https://www.encodeproject.org/publication-data/ENCSR876DCD/ |
| RBNS Computational Pipeline (PDF) | https://www.encodeproject.org/documents/c8b3442a-7e63-4847-af11-c72597bf65b3/@@download/attachment/RBNS_Computational_Pipeline_Aug_2016_update_Dec2018.pdf |
| RBNS Experimental Protocol (PDF) | https://www.encodeproject.org/documents/aa71cabf-aaee-4358-a834-c6ee002938b8/@@download/attachment/RBNSExperimentalProtocol_Feb2016_96well.pdf |
| ENCODE REST API | https://www.encodeproject.org/help/rest-api/ |
