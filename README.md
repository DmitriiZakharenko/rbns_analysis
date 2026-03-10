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
| RBPs with binding data | **88** |
| Total examples | **246,777** |
| Positives (binding_label=1) | **74,135** (30.1%) |
| Negatives (binding_label=0) | **172,642** (69.9%) |
| High-confidence positives | **64,688** (87.3% of positives) |
| Median R_max (enrichment) | **15.8** |
| Sequence length | 20 nt (uniform) |
| Alphabet | RNA (A, C, G, U) |
| Duplicates / label overlaps | 0 / 0 |

RBPs excluded:
- **9 RBPs** — no enrichment signal after filtering (CELF1, HNRNPD, HNRNPF, IGF2BP1, MBNL1, MSI1, RBFOX2, RBM47, SRSF2)
- **Several RBPs** — no matched input control library in ENCODE (e.g. APOBEC3C, PTBP3)

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
| 4 | `04_compute_enrichment.py` | Compute R = f_pulldown / f_input across all concentrations |
| 5 | `05_build_ml_dataset.py` | Merge positives and sample negatives per RBP |
| 6 | `06_clean_dataset.py` | Drop missing values, fix types, remove duplicates |
| 7 | `07_validate_dataset.py` | Quality checks, label distribution, stats JSON |

---

## Enrichment Strategy

For each RBP and each pulldown concentration:

```
R(seq, conc) = f_pulldown(conc) / f_input
R_max(seq)   = max R across all concentrations
```

Positives: sequences with `R_max >= min_R` (default 1.5), top `top_k_positive` (default 1000) per RBP.

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
| `rna_sequence` | 20 nt RNA sequence (A/C/G/U) |
| `binding_label` | 1 = positive, 0 = negative |
| `source` | `enriched` or `background` |
| `R_max` | Maximum enrichment ratio (positives only) |
| `n_enriched_concs` | Concentrations with R >= min_R |
| `n_concs_measured` | Total concentrations measured |
| `high_confidence` | 1 if high-confidence positive |

---

## Notes on Reproducibility

Raw FASTQ files (~260 GB total) are not included in this repository.
To reproduce from scratch, run the full pipeline starting from step 1.
Step 2 downloads all required FASTQ files directly from ENCODE.

---

## Data Sources

- **ENCODE RBNS collection**: [ENCSR876DCD](https://www.encodeproject.org/publication-data/ENCSR876DCD/)
- **Publication**: Lambert et al., *Nature* 2020 — [doi:10.1038/s41586-020-2077-3](https://www.nature.com/articles/s41586-020-2077-3)
- **RBNS Computational Pipeline**: [ENCODE document (PDF)](https://www.encodeproject.org/documents/c8b3442a-7e63-4847-af11-c72597bf65b3/@@download/attachment/RBNS_Computational_Pipeline_Aug_2016_update_Dec2018.pdf)
- **ENCODE REST API**: [https://www.encodeproject.org/help/rest-api/](https://www.encodeproject.org/help/rest-api/)
