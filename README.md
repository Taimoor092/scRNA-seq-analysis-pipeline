# scRNA-seq Analysis Pipeline 🧬

[![R Version](https://img.shields.io/badge/R-%3E%3D4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5-276DC3?style=flat-square)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](LICENSE)
[![Status](https://img.shields.io/badge/status-active-success?style=flat-square)]()

A reproducible, modular end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline built in R using **Seurat v5**. Covers everything from raw count matrix loading through quality control, normalization, dimensionality reduction, clustering, cell type annotation, and publication-ready visualization.

---

## 📁 Project Structure

```
scRNA-seq-pipeline/
├── README.md
├── config/
│   └── params.yaml              # All tunable parameters in one place
├── R/
│   ├── 01_load_data.R           # Load 10x or matrix data
│   ├── 02_quality_control.R     # QC filtering (MT%, nFeature, nCount)
│   ├── 03_normalization.R       # Normalization & feature selection
│   ├── 04_dim_reduction.R       # PCA, UMAP, tSNE
│   ├── 05_clustering.R          # Graph-based clustering
│   ├── 06_annotation.R          # Marker genes & cell type annotation
│   ├── 07_visualization.R       # Publication-ready plots
│   └── utils.R                  # Shared helper functions
├── scripts/
│   ├── run_pipeline.R           # Master script — runs all steps
│   └── install_packages.R       # One-shot package installer
├── data/
│   ├── raw/                     # Input: 10x matrix files or .rds
│   └── processed/               # Output: Seurat objects per step
├── results/
│   ├── figures/                 # All plots (PDF + PNG)
│   ├── tables/                  # DEG tables, marker genes, metadata
│   └── reports/                 # HTML reports
└── docs/
    └── methods.md               # Methods text for manuscript/thesis
```

---

## 🚀 Quick Start

### 1. Clone the repository

```bash
git clone https://github.com/YOUR_USERNAME/scRNA-seq-pipeline.git
cd scRNA-seq-pipeline
```

### 2. Install dependencies

```r
source("scripts/install_packages.R")
```

### 3. Add your data

Place your 10x Genomics output files in `data/raw/`:

```
data/raw/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

Or use the built-in demo dataset (PBMC 3k from 10x Genomics):

```r
# In config/params.yaml, set:
# use_demo_data: true
```

### 4. Configure parameters

Edit `config/params.yaml` to match your experiment:

```yaml
# Key parameters to review before running
qc:
  min_features: 200     # Minimum genes per cell
  max_features: 5000    # Maximum genes (filter doublets)
  max_mt_percent: 20    # Maximum mitochondrial %
clustering:
  resolution: 0.5       # Higher = more clusters
```

### 5. Run the full pipeline

```r
source("scripts/run_pipeline.R")
```

Or run steps individually:

```r
source("R/01_load_data.R")
source("R/02_quality_control.R")
# ... etc
```

---

## 📊 Pipeline Overview

```
Raw Count Matrix
      │
      ▼
 01. Load Data ──────────── 10x Genomics / MEX format
      │
      ▼
 02. Quality Control ─────── MT%, nFeature, nCount filtering
      │
      ▼
 03. Normalization ──────── NormalizeData + FindVariableFeatures
      │
      ▼
 04. Dimensionality Reduction ── PCA → UMAP / tSNE
      │
      ▼
 05. Clustering ─────────── FindNeighbors + FindClusters
      │
      ▼
 06. Annotation ─────────── Marker genes → Cell type labels
      │
      ▼
 07. Visualization ──────── UMAP, violin, dot, feature plots
      │
      ▼
  Results & Report
```

---

## 📈 Example Output

| Output | Description | Location |
|--------|-------------|----------|
| UMAP plot | Cells colored by cluster | `results/figures/umap_clusters.pdf` |
| UMAP by cell type | Annotated cell populations | `results/figures/umap_celltypes.pdf` |
| Violin plots | QC metrics before/after filtering | `results/figures/qc_violin.pdf` |
| Dot plot | Top marker genes per cluster | `results/figures/dotplot_markers.pdf` |
| Marker table | DEGs per cluster (Wilcoxon) | `results/tables/cluster_markers.csv` |
| Seurat object | Final annotated object | `data/processed/seurat_annotated.rds` |

---

## 🧰 Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| Seurat | ≥ 5.0 | Core scRNA-seq analysis |
| SeuratObject | ≥ 5.0 | Data structures |
| ggplot2 | ≥ 3.4 | Visualization |
| dplyr | ≥ 1.1 | Data manipulation |
| patchwork | ≥ 1.1 | Plot composition |
| yaml | ≥ 2.3 | Config parsing |
| Matrix | ≥ 1.6 | Sparse matrix handling |
| BiocManager | ≥ 1.30 | Bioconductor packages |
| SingleR | ≥ 2.0 | Automated cell annotation |
| celldex | ≥ 1.10 | Reference datasets for SingleR |
| scran | ≥ 1.28 | Additional normalization |

---

## 🔬 Demo Dataset

This pipeline uses the **PBMC 3k dataset** from 10x Genomics as the default demo:

- **3,000 peripheral blood mononuclear cells** from a healthy donor
- Sequenced on Illumina NextSeq 500
- Available free from [10x Genomics](https://www.10xgenomics.com/resources/datasets)

Expected cell types: T cells (CD4+, CD8+), B cells, NK cells, Monocytes, Dendritic cells

---

## ⚙️ Configuration Reference

All parameters live in `config/params.yaml`. No need to edit R scripts for routine analyses.

```yaml
project:
  name: "PBMC_analysis"
  use_demo_data: true     # Set false to use your own data

data:
  input_dir: "data/raw"
  output_dir: "data/processed"
  sample_name: "sample1"

qc:
  min_features: 200
  max_features: 5000
  min_counts: 500
  max_mt_percent: 20

normalization:
  method: "LogNormalize"  # or "SCTransform"
  scale_factor: 10000
  n_variable_features: 2000

dim_reduction:
  n_pcs: 30
  run_tsne: false
  umap_dims: 1:20

clustering:
  resolution: 0.5
  algorithm: 1            # 1=Louvain, 2=Louvain refined, 3=SLM, 4=Leiden

annotation:
  method: "singler"       # or "manual"
  reference: "HumanPrimaryCellAtlasData"

output:
  figures_dir: "results/figures"
  tables_dir: "results/tables"
  save_formats: ["pdf", "png"]
  figure_width: 10
  figure_height: 8
```

---

## 📖 Methods

A ready-to-use methods section for your thesis or paper is in [`docs/methods.md`](docs/methods.md).

---

## 🤝 Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

---

## 📄 License

MIT License — see [LICENSE](LICENSE) for details.

---

## 📚 References

- Hao et al. (2024). *Dictionary learning for integrative, multimodal and scalable single-cell analysis.* Nature Methods.
- Stuart et al. (2019). *Comprehensive Integration of Single-Cell Data.* Cell.
- 10x Genomics PBMC 3k dataset: https://www.10xgenomics.com/resources/datasets

---

*Built with ❤️ for reproducible bioinformatics | [Your Name] | MSc Bioinformatics*
