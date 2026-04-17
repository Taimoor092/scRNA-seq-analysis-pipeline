# Methods

## Single-Cell RNA Sequencing Data Analysis

### Data Loading and Pre-processing

Raw gene expression matrices were loaded using the `Read10X()` function from
the Seurat package (v5.0) in R (v≥4.3). A Seurat object was created with
genes expressed in fewer than 3 cells excluded from downstream analysis.

### Quality Control

Cell-level quality control (QC) was performed based on three metrics:
(1) number of detected genes per cell (`nFeature_RNA`),
(2) total UMI counts per cell (`nCount_RNA`), and
(3) percentage of reads mapping to mitochondrial genes (`percent.mt`).

Cells were retained if they had between **200 and 5,000** detected genes,
more than **500** UMI counts, and fewer than **20%** mitochondrial reads.
These thresholds were selected based on the distribution of QC metrics
visualised via violin plots and scatter plots prior to filtering.

### Normalisation and Feature Selection

Gene expression counts were normalised using the `NormalizeData()` function
with the "LogNormalize" method and a scale factor of 10,000, yielding
log-normalised counts: log(counts / total counts per cell × 10,000 + 1).
The top **2,000 highly variable genes** (HVGs) were identified using the
variance-stabilising transformation (VST) method via `FindVariableFeatures()`.
Data were subsequently scaled to zero mean and unit variance across cells
using `ScaleData()`.

### Dimensionality Reduction

Principal component analysis (PCA) was performed on the HVG-scaled expression
matrix using `RunPCA()`. The number of principal components (PCs) for
downstream analyses was selected by inspecting an elbow plot of the standard
deviation explained per PC. UMAP (Uniform Manifold Approximation and
Projection) was computed using `RunUMAP()` with the selected PCs, a minimum
distance of 0.3, and 30 nearest neighbours.

### Clustering

Cell-to-cell KNN graphs were constructed using `FindNeighbors()` based on
the selected PCA dimensions. Cell clusters were identified using the Louvain
community detection algorithm (`FindClusters()`, resolution = 0.5).

### Cell Type Annotation

Cluster identities were assigned using automated annotation with
**SingleR** (v2.0), comparing cluster expression profiles against the
Human Primary Cell Atlas reference dataset (`celldex::HumanPrimaryCellAtlasData()`).
The majority cell type label per cluster was assigned as the cluster identity.

### Differential Expression and Marker Genes

Marker genes for each cluster were identified using the Wilcoxon rank-sum
test via `FindAllMarkers()`, retaining genes expressed in at least 10% of
cells in a cluster and with an absolute log₂ fold-change ≥ 0.25.
Significance threshold: adjusted p-value < 0.05 (Bonferroni correction).

### Visualisation

All plots were generated using **ggplot2** (v3.4) and Seurat's built-in
plotting functions. Figure panels were assembled using **patchwork**.
UMAP embeddings, violin plots, dot plots, feature plots, and heatmaps
were produced for publication.

### Reproducibility

All analyses were performed in R (v≥4.3.0). The random seed was set to 42
throughout. The complete analysis pipeline, including configuration files
and all R scripts, is publicly available at:
`https://github.com/YOUR_USERNAME/scRNA-seq-pipeline`

### Software Versions

| Package | Version |
|---------|---------|
| R | ≥ 4.3.0 |
| Seurat | ≥ 5.0.0 |
| SingleR | ≥ 2.0.0 |
| celldex | ≥ 1.10.0 |
| ggplot2 | ≥ 3.4.0 |
| patchwork | ≥ 1.1.0 |
| dplyr | ≥ 1.1.0 |
