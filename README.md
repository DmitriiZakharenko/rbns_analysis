# RBNS (RNA Bind-N-Seq) — ML Dataset Pipeline

End-to-end pipeline for building an RBP-RNA binding dataset from
**RNA Bind-N-Seq (RBNS)** experiments published in ENCODE (Lambert et al., *Nature* 2020).

- **Positives**: RNA sequences enriched in protein pulldown across multiple concentrations.
- **Negatives**: sequences from the 0 nM input (background) library.
- **Output**: `results/ml_dataset_rbns_clean.tsv`

---

## Dataset Summary

| Metric | Value |
|--------|-------|
| FASTQ files processed | 601 |
| RBPs with validated binding data | **96** |
| Total examples | **284,642** |
| Positives (binding_label=1) | **96,000** (33.7%) |
| Negatives (binding_label=0) | **188,642** (66.3%) |
| High-confidence positives | **96,000** (100%; all R_max ≥ 1.5) |
| Median R_max (k-mer enrichment score) | **3.69** |
| Mean R_max | **6.43** |
| Min / max R_max (positives) | **1.61** / **51.6** |
| Sequence length | 20 nt (uniform) |
| Alphabet | RNA (A, C, G, U) |
| Duplicates / label overlaps | 0 / 0 |

Of 110 RBPs with FASTQ data in ENCODE:
- **96 RBPs** have a matched 0 nM input control and detectable k-mer enrichment → included.
- **11 RBPs** have no matched input control in ENCODE → excluded (R cannot be computed reliably).
- **3 RBPs** have a matched control but no k-mer with R ≥ 1.5 at any concentration → excluded.

**Note on read lengths**: some ENCODE RBNS experiments use 40 nt reads (e.g. CELF1, RBFOX2, MBNL1).
The random RNA region is always 20 nt; the pipeline trims longer reads to 20 nt automatically.

See `results/dataset_stats_rbns.json` and `results/validation_summary_rbns.tsv` for full statistics.

---

## Quick Start

**Requirements**: Python 3.9+

```bash
pip install -r requirements.txt
```

**Run the full pipeline:**

```bash
bash run_pipeline.sh
```

**Start from a specific step** (e.g. if FASTQ already downloaded):

```bash
bash run_pipeline.sh --from 3
```

**Run a single step:**

```bash
bash run_pipeline.sh --only 4
```

Configuration: `config.yaml`.

---

## Pipeline Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_fetch_encode_metadata.py` | Fetch experiment and file metadata from ENCODE API |
| 2 | `02_download_fastq.py` | Download FASTQ files from ENCODE |
| 3 | `03_process_fastq.py` | Convert FASTQ to sequence count TSV (T to U, min_count=2) |
| 4 | `04_compute_enrichment.py` | K-mer enrichment (default k=5) and top positives per RBP; `--mode sequence` for per-20-mer R |
| 5 | `05_build_ml_dataset.py` | Merge positives and sample negatives per RBP |
| 6 | `06_clean_dataset.py` | Drop missing values, fix types, remove duplicates |
| 7 | `07_validate_dataset.py` | Quality checks, label distribution, stats JSON |
| 8 *(optional)* | `08_add_protein_sequences.py` | Add `protein_sequence` column from UniProt API (run with internet access) |

---

## Enrichment Strategy

For each RBP, each pulldown concentration, and each k-mer (k=5 by default):

```
R(kmer, conc) = (kmer_count_pulldown / total_reads_pulldown)
              / (kmer_count_input    / total_reads_input)

score(seq, conc) = max R(kmer, conc)  over all k-mers in seq
R_max(seq)       = max score(seq, conc) over all concentrations
```

This matches the ENCODE RBNS computational pipeline (Lambert et al. 2020):
aggregating reads across all sequences sharing a k-mer gives far more
statistical power than per-sequence frequency ratios.

For each 20-mer sequence, `R_max = max over concentrations of [max k-mer R within the sequence]`.
Positives: sequences with `R_max >= min_R` (default 1.5), up to `top_k_positive` (default 1000) per RBP.
RBPs with fewer than `min_positives` (default 10) such sequences are excluded.

Quality metrics per sequence:

| Metric | Description |
|--------|-------------|
| `R_max` | Maximum enrichment ratio across concentrations |
| `n_enriched_concs` | Number of concentrations where R >= min_R |
| `n_concs_measured` | Total concentrations available for this RBP |
| `is_monotonic` | 1 if R increases monotonically with concentration |
| `high_confidence` | 1 if R_max >= min_R and n_enriched_concs >= 1 |

---

## Output Files

```
results/
├── ml_dataset_rbns_clean.tsv    # final cleaned dataset
├── ml_dataset_rbns.tsv          # raw merged dataset
├── dataset_stats_rbns.json      # summary statistics
├── validation_summary_rbns.tsv  # per-RBP counts
└── tables/
    └── {RBP}_positives.tsv      # per-RBP enriched sequences with metrics
```

### Dataset columns

| Column | Description |
|--------|-------------|
| `target_name` | RBP name (e.g. `IGF2BP2`) |
| `protein_sequence` | Canonical UniProt amino-acid sequence — **added by step 8** (absent until `08_add_protein_sequences.py` is run) |
| `rna_sequence` | 20 nt RNA sequence (A/C/G/U) |
| `binding_label` | 1 = positive, 0 = negative |
| `source` | `enriched` or `background` |
| `R_max` | Max k-mer enrichment ratio across concentrations (positives only; NaN for negatives) |
| `n_enriched_concs` | Concentrations with R >= min_R |
| `n_concs_measured` | Total concentrations where sequence was observed (min_count ≥ 2) |
| `high_confidence` | 1 if R_max ≥ min_R and n_enriched_concs ≥ 1 |

**Note on R_max**: for 83/96 RBPs all 1000 positives share the same R_max value (the enrichment score of the dominant binding k-mer). R_max is a reliable confidence indicator for the positive threshold but is not a discriminating feature between positives within a single RBP.

---

## Notes on Reproducibility

Raw FASTQ files (~150 GB total) are not included in this repository.
To reproduce from scratch, run the full pipeline starting from step 1.
Step 2 downloads all required FASTQ files directly from ENCODE.

---

## Data Sources

- **ENCODE RBNS collection**: [ENCSR876DCD](https://www.encodeproject.org/publication-data/ENCSR876DCD/)
- **Publication**: Lambert et al., *Nature* 2020 — [doi:10.1038/s41586-020-2077-3](https://www.nature.com/articles/s41586-020-2077-3)
- **RBNS Computational Pipeline**: [ENCODE document (PDF)](https://www.encodeproject.org/documents/c8b3442a-7e63-4847-af11-c72597bf65b3/@@download/attachment/RBNS_Computational_Pipeline_Aug_2016_update_Dec2018.pdf)
- **ENCODE REST API**: [https://www.encodeproject.org/help/rest-api/](https://www.encodeproject.org/help/rest-api/)
