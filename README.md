# scRNA-seq HVG Parameter Evaluation Toolkit

A systematic framework for evaluating and selecting optimal highly variable gene (HVG) parameters in single-cell RNA-seq analysis. Rather than guessing a single layer/flavor/n_top_genes combination, this toolkit sweeps across a configurable parameter grid and scores every configuration against two biologically meaningful criteria:

- **Marker gene coverage** — what fraction of known cell-type markers are retained in the HVG set?
- **Nuisance gene contamination** — how many HVG slots are wasted on mitochondrial, ribosomal, hemoglobin, or other uninformative genes?

A composite quality score then recommends the best configuration, which can be applied to the AnnData object in a single call.

Built on top of [Scanpy](https://scanpy.readthedocs.io/) and designed for use in Jupyter notebooks.

## Table of Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Module Reference](#module-reference)
  - [marker_genes.py](#marker_genespy)
  - [hvg_evaluation.py](#hvg_evaluationpy)
- [Customisation Guide](#customisation-guide)
  - [Writing Your Own Marker Dictionary](#writing-your-own-marker-dictionary)
  - [Changing the Parameter Grid](#changing-the-parameter-grid)
  - [Adding Custom Nuisance Gene Categories](#adding-custom-nuisance-gene-categories)
  - [Tuning the Recommendation Weights](#tuning-the-recommendation-weights)
- [Output Description](#output-description)
- [Example Notebook](#example-notebook)
- [Dependencies](#dependencies)

## Overview

### The Problem

Scanpy provides several HVG selection flavors (`seurat`, `seurat_v3`, `cell_ranger`, `pearson_residuals`), each operating on different data representations (raw counts vs. normalised values). The choice of:

- **Layer** (raw counts, log-normalised, scran-normalised, SoupX-corrected, etc.)
- **Flavor** (statistical method for ranking gene variability)
- **n_top_genes** (how many HVGs to keep)
- **batch_key** (batch-aware vs. batch-naive selection)

...can meaningfully affect downstream clustering and cell-type resolution. This toolkit makes that choice data-driven.

### The Approach

```
┌─────────────────────────────────────────────────────────┐
│                    Parameter Grid                        │
│  layers × flavors × n_top_genes × batch_keys            │
└──────────────────────┬──────────────────────────────────┘
                       │
                       ▼
              ┌────────────────┐
              │  select_hvgs() │  (one call per combination)
              └───────┬────────┘
                      │
           ┌──────────┴──────────┐
           ▼                     ▼
  ┌─────────────────┐   ┌────────────────────┐
  │ Marker coverage │   │ Nuisance gene count │
  │ per cell type   │   │ mito/ribo/hb/custom │
  └────────┬────────┘   └─────────┬──────────┘
           │                      │
           └──────────┬───────────┘
                      ▼
            ┌───────────────────┐
            │ Composite scoring │
            │ & recommendation  │
            └─────────┬─────────┘
                      ▼
            ┌───────────────────┐
            │ apply_optimal_hvgs│
            │ → adata.var[      │
            │ 'highly_variable']│
            └───────────────────┘
```

For each parameter combination, the toolkit:

1. Runs HVG selection via the appropriate Scanpy API
2. Measures what fraction of your reference marker genes are captured
3. Counts how many nuisance genes (mitochondrial, ribosomal, hemoglobin, plus any custom categories) infiltrate the HVG set
4. Computes a composite quality score and recommends the best configuration

## Installation

Clone the repository and ensure the dependencies are available in your Python environment:

```bash
git clone https://github.com/<your-username>/hvg-evaluation.git
cd hvg-evaluation
pip install scanpy anndata matplotlib seaborn numpy pandas
```

Both `marker_genes.py` and `hvg_evaluation.py` should be placed in the same directory as your analysis notebook (or on your `PYTHONPATH`).

## Quick Start

```python
import scanpy as sc
from marker_genes import COARSE_MARKERS
from hvg_evaluation import (
    run_hvg_sweep,
    recommend_optimal,
    apply_optimal_hvgs,
    plot_marker_coverage_heatmap,
    plot_nuisance_heatmap,
    plot_quality_summary,
    plot_missing_markers_table,
)

# Load your AnnData object
adata = sc.read_h5ad("my_data.h5ad")

# 1. Run the parameter sweep
hvg_sweep_df = run_hvg_sweep(
    adata,
    marker_dict=COARSE_MARKERS,
)

# 2. Visualise results
plot_marker_coverage_heatmap(hvg_sweep_df)
plot_nuisance_heatmap(hvg_sweep_df)
plot_quality_summary(hvg_sweep_df)

# 3. Get the optimal recommendation
best = recommend_optimal(
    hvg_sweep_df,
    marker_weight=0.6,
    nuisance_weight=0.4,
    min_marker_coverage=0.5,
    max_nuisance_pct=10.0,
)

# 4. Apply to your AnnData
optimal_hvgs = apply_optimal_hvgs(adata, best)

# 5. Inspect which markers are still missing
missing_df = plot_missing_markers_table(adata, optimal_hvgs)
```

After step 4, `adata.var['highly_variable']` is set and ready for `sc.pp.pca()` and downstream analysis.

## Module Reference

### marker_genes.py

The marker gene dictionary serves as the biological ground truth for evaluating HVG quality. The provided file now contains a **small, generic example dictionary** (broad epithelial / stromal / immune / vascular compartments) intended as a **template**. In almost all real analyses you should **copy and customise** these markers for your own tissue, species and question.

#### Exported Objects

| Name | Type | Description |
|---|---|---|
| `COARSE_MARKERS` | `Dict[str, List[str]]` | Example broad compartments (epithelial, fibroblasts, endothelial, immune, etc.), 8–20 genes each |
| `filter_markers_to_adata()` | function | Filters a marker dict to genes present in the dataset |

#### Dictionary Format

Every marker dictionary follows the same structure:

```python
MARKERS: Dict[str, List[str]] = {
    "Cell_Type_Name": ["GENE1", "GENE2", "GENE3", ...],
    "Another_Type":   ["GENEA", "GENEB", ...],
}
```

- **Keys** are cell type names (strings). Use underscores for multi-word names.
- **Values** are lists of gene symbols (strings). Use the same gene symbol convention as your AnnData's `var_names` (typically HGNC symbols for human data).
- Aim for **8-25 genes** per cell type at the coarse level. Too few gives unstable coverage estimates; too many dilutes specificity.
- Genes present in multiple cell types are allowed (biological co-expression is real), but minimising overlap improves discriminative power.

#### Tips for choosing marker genes

- **Start from biology you trust**: atlases, CellMarker, PanglaoDB and well-cited papers for your tissue.
- **Favour strong enrichment**: genes that are clearly higher in one cell type than all others (e.g. `CD3D` for T cells, `PECAM1` for endothelial cells).
- **Avoid housekeeping / ubiquitous genes**: e.g. `ACTB`, `GAPDH`, `RPL*`, `RPS*` – they add noise to the coverage fraction without helping identity calls.
- **Avoid cell-cycle, stress and immediate-early genes**: e.g. `MKI67`, `TOP2A`, `FOS`, `JUN`, `HSPA1A`. These reflect state rather than identity and are better treated as nuisance.
- **Consider removing extreme structural genes from the marker dict** if they already dominate your dataset and are not helpful for distinguishing *between* cell types in your question (e.g. sarcomeric genes in muscle-only data, collagen super-high expressers in stromal-only data). Those can instead be tracked as custom nuisance categories.

You should treat `COARSE_MARKERS` as a **living file**: for a new project, copy `marker_genes.py`, rename the cell-type keys and adjust the gene lists for your biology.

### hvg_evaluation.py

The core evaluation engine. Imports `COARSE_MARKERS` and `filter_markers_to_adata` from `marker_genes.py`.

#### Constants

| Name | Default Value | Description |
|---|---|---|
| `DEFAULT_LAYER_FLAVORS` | See below | Layer-flavor pairings to test |
| `DEFAULT_N_TOP_GENES` | `[1000, 1500, 2000, 2500, 3000, 4000, 5000]` | HVG counts to sweep |
| `DEFAULT_BATCH_KEYS` | `["donor", None]` | Batch keys to test (`None` = batch-naive) |

Default layer-flavor pairings:

```python
DEFAULT_LAYER_FLAVORS = [
    ("soupX_counts", "seurat_v3"),        # raw counts → seurat_v3
    ("soupX_counts", "pearson_residuals"), # raw counts → pearson_residuals
    ("log1p_norm",   "seurat"),            # log-normalised → seurat
    ("log1p_norm",   "cell_ranger"),       # log-normalised → cell_ranger
    ("scran_log1p",  "seurat"),            # scran-normalised → seurat
    ("scran_log1p",  "cell_ranger"),       # scran-normalised → cell_ranger
]
```

Layer names that are not found in your AnnData are silently skipped — only layers present in `adata.layers` are evaluated.

#### Core Functions

**`select_hvgs(adata, layer, flavor, n_top_genes, batch_key=None) → Set[str]`**

Runs HVG selection for a single parameter set and returns the gene names. Handles Scanpy API differences between flavors internally:
- `seurat_v3` and `pearson_residuals` accept a `layer` parameter directly
- `seurat` and `cell_ranger` require the layer to be placed in `adata.X`

The input AnnData is not mutated (a copy is made internally).

**`compute_marker_coverage(hvg_set, marker_dict, var_names) → dict`**

Calculates per-cell-type and overall marker coverage. Returns a dictionary with:
- `per_type`: dict of `{cell_type: {n_in_data, n_in_hvg, fraction, missing}}`
- `overall_n_in_data`, `overall_n_in_hvg`, `overall_fraction`

**`classify_nuisance_genes(adata, custom_nuisance=None) → Dict[str, Set[str]]`**

Identifies nuisance gene sets. Three built-in categories are always included:
- **Mitochondrial** (`MT-` prefix or `adata.var["mt"]` boolean column)
- **Ribosomal** (`RPS`/`RPL` prefix or `adata.var["ribo"]` boolean column)
- **Hemoglobin** (`HB[ABDEQGMZ]` prefix or `adata.var["hb"]` boolean column)

Custom categories can be added via the `custom_nuisance` parameter (see [Adding Custom Nuisance Gene Categories](#adding-custom-nuisance-gene-categories)).

**`compute_nuisance_counts(hvg_set, nuisance_sets) → dict`**

Counts how many genes from each nuisance category appear in a given HVG set. Returns per-category counts (`n_mito`, `pct_mito`, etc.) and aggregate totals (`n_nuisance`, `pct_nuisance`).

**`run_hvg_sweep(adata, marker_dict=None, layer_flavors=None, n_top_genes_range=None, batch_keys=None, custom_nuisance=None, verbose=True) → pd.DataFrame`**

The main entry point. Sweeps over the full parameter grid and returns a tidy DataFrame with one row per configuration. Each row contains:
- Parameter columns: `layer`, `flavor`, `batch_key`, `n_top_genes`
- Coverage columns: `overall_fraction`, `overall_n_in_hvg`, `frac_<CellType>` per marker dict entry
- Nuisance columns: `n_mito`, `pct_mito`, `n_ribo`, `pct_ribo`, ... `n_nuisance`, `pct_nuisance`
- A human-readable `condition` label

**`recommend_optimal(results_df, marker_weight=0.6, nuisance_weight=0.4, min_marker_coverage=0.0, max_nuisance_pct=100.0, verbose=True) → pd.Series`**

Scores every configuration with a composite quality metric:

```
score = marker_weight × overall_fraction − nuisance_weight × (pct_nuisance / 100)
```

Optional hard thresholds (`min_marker_coverage`, `max_nuisance_pct`) filter out unacceptable configurations before scoring. Returns the best row as a `pd.Series` with an added `hvg_quality_score` field.

**`apply_optimal_hvgs(adata, best) → Set[str]`**

Re-runs HVG selection with the recommended parameters and writes `adata.var['highly_variable']` in place. Returns the set of selected gene names.

#### Plotting Functions

All plotting functions accept an optional `save_path` parameter to save the figure to disk.

| Function | Description |
|---|---|
| `plot_marker_coverage_heatmap()` | Heatmap: overall marker fraction, conditions × n_top_genes |
| `plot_marker_coverage_lines()` | Line plot: overall marker fraction vs. n_top_genes per condition |
| `plot_per_celltype_coverage()` | Grouped bar chart: per-cell-type coverage at a given n_top_genes |
| `plot_nuisance_heatmap()` | Heatmap: nuisance gene percentage, conditions × n_top_genes |
| `plot_nuisance_lines()` | Line plot: nuisance % vs. n_top_genes per condition |
| `plot_nuisance_breakdown()` | Stacked bar chart: per-category nuisance counts for every condition |
| `plot_quality_summary()` | Two-panel summary: marker coverage (higher=better) vs. nuisance (lower=better) |
| `plot_missing_markers_table()` | Returns a DataFrame listing which markers are missing from a given HVG set |

## Customisation Guide

### Writing Your Own Marker Dictionary

To use this toolkit for a different tissue or organism, either **edit the provided `marker_genes.py` directly** or create your own with at minimum a `COARSE_MARKERS` dictionary and the `filter_markers_to_adata` helper:

```python
"""Marker genes for [your tissue/organism]."""
from typing import Dict, List

COARSE_MARKERS: Dict[str, List[str]] = {
    "Epithelial": [
        "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1",
        "CLDN4", "CLDN7", "MUC1",
    ],
    "Fibroblasts": [
        "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA",
        "FAP", "THY1", "VIM",
    ],
    "Endothelial": [
        "PECAM1", "CDH5", "VWF", "ERG", "FLT1",
        "KDR", "CLDN5", "EMCN",
    ],
    "T_Cells": [
        "CD3D", "CD3E", "CD3G", "TRAC",
        "IL7R", "CD4", "CD8A", "CD8B",
    ],
    # ... add more cell types as needed
}


def filter_markers_to_adata(
    marker_dict: Dict[str, List[str]],
    var_names,
) -> Dict[str, List[str]]:
    """Filter marker dictionary to genes present in the dataset."""
    var_set = set(var_names)
    filtered = {}
    for cell_type, genes in marker_dict.items():
        valid = [g for g in genes if g in var_set]
        if len(valid) >= 2:
            filtered[cell_type] = valid
    return filtered
```

**Guidelines:**

- Use the same gene symbol format as your AnnData `var_names` (e.g., HGNC for human, MGI for mouse)
- Aim for 8-25 genes per cell type — enough for stable scoring, not so many that specificity is diluted
- Prefer well-validated, highly specific markers over broadly expressed genes
- Minimise overlap between cell types to improve discriminative power
- Periodically re-check marker sets after annotation: if a gene turns out to be broadly expressed or driven by batch / QC artefacts rather than biology, remove it from the marker dict.
- The `filter_markers_to_adata` function automatically drops genes absent from the dataset and removes cell types with fewer than 2 detectable markers

### Changing the Parameter Grid

Override any axis of the sweep by passing arguments to `run_hvg_sweep()`:

```python
# Custom layers and flavors
my_layer_flavors = [
    ("raw_counts", "seurat_v3"),
    ("raw_counts", "pearson_residuals"),
    ("log_norm",   "seurat"),
]

# Custom HVG counts
my_n_top = [500, 1000, 2000, 3000, 5000]

# Custom batch keys (e.g., sample-level)
my_batch_keys = ["sample_id", None]

hvg_sweep_df = run_hvg_sweep(
    adata,
    marker_dict=MY_MARKERS,
    layer_flavors=my_layer_flavors,
    n_top_genes_range=my_n_top,
    batch_keys=my_batch_keys,
)
```

**Layer-flavor compatibility rules:**
- `seurat_v3` and `pearson_residuals` expect **raw counts** (integers)
- `seurat` and `cell_ranger` expect **normalised values** (typically log-transformed)

Mismatched combinations will fail gracefully and be skipped with a warning.

### Adding Custom Nuisance Gene Categories

Three categories are always checked (mitochondrial, ribosomal, hemoglobin). To add **tissue- and experiment-specific** nuisance genes, pass a dictionary mapping category names to gene lists:

```python
custom_nuisance = {
    "myofibre": [
        "TTN", "NEB", "DES", "DMD", "ACTA1", "RYR1",
        "MYH7", "MYH1", "MYH2",
    ],
    "epithelial_contamination": [
        "KRT1", "KRT5", "KRT10", "KRT14",
    ],
}

hvg_sweep_df = run_hvg_sweep(
    adata,
    custom_nuisance=custom_nuisance,
)
```

Custom categories appear as additional `n_<name>` and `pct_<name>` columns in the results DataFrame, and are included in the aggregate `n_nuisance` / `pct_nuisance` totals used by the recommendation scoring.

#### Tips for selecting nuisance genes

- **Think about “wasted HVG slots”**: genes that are always highly expressed but tell you little about **which** cell type a cluster represents.
- Typical examples:
  - Structural genes that dominate one compartment but are not informative within it (e.g. contractile apparatus genes in muscle, keratins in skin, collagen super-expressers in fibroblast-only datasets).
  - Ambient contamination signatures (e.g. epithelial keratins in mostly stromal samples, plasma proteins like `ALB`, `APOA1` in solid tissues).
  - Technical artefact signatures specific to your protocol, if any are known.
- Avoid adding **true identity markers** (e.g. `CD3D` for T cells) to nuisance lists, or you will unfairly penalise configurations that correctly keep these markers.
- Start with a small, conservative nuisance list, inspect which nuisance genes are retained at the recommended configuration, and iteratively refine.

### Tuning the Recommendation Weights

The `recommend_optimal()` function uses a composite score:

```
score = marker_weight × overall_fraction − nuisance_weight × (pct_nuisance / 100)
```

Both `overall_fraction` and `pct_nuisance / 100` are on a 0-1 scale, so the weights directly control the trade-off:

```python
# Default: balanced towards marker coverage
best = recommend_optimal(hvg_sweep_df, marker_weight=0.6, nuisance_weight=0.4)

# Strict: require high marker coverage and low nuisance
best = recommend_optimal(
    hvg_sweep_df,
    marker_weight=0.6,
    nuisance_weight=0.4,
    min_marker_coverage=0.7,  # reject configs below 70% coverage
    max_nuisance_pct=5.0,     # reject configs above 5% nuisance
)

# Prioritise marker coverage above all else
best = recommend_optimal(hvg_sweep_df, marker_weight=1.0, nuisance_weight=0.0)

# Prioritise clean HVG sets (minimal nuisance)
best = recommend_optimal(hvg_sweep_df, marker_weight=0.3, nuisance_weight=0.7)
```

## Output Description

### Sweep DataFrame Columns

| Column | Description |
|---|---|
| `layer` | AnnData layer used |
| `flavor` | Scanpy HVG flavor |
| `batch_key` | Batch column or `"none"` |
| `n_top_genes` | Requested number of HVGs |
| `n_hvgs_actual` | Actual number of HVGs returned |
| `overall_fraction` | Fraction of all unique marker genes found in HVGs (0-1) |
| `overall_n_in_data` | Number of unique marker genes present in the dataset |
| `overall_n_in_hvg` | Number of those marker genes captured by HVGs |
| `frac_<CellType>` | Per-cell-type marker coverage fraction |
| `n_mito`, `pct_mito` | Mitochondrial gene count and percentage in HVGs |
| `n_ribo`, `pct_ribo` | Ribosomal gene count and percentage |
| `n_hb`, `pct_hb` | Hemoglobin gene count and percentage |
| `n_<custom>`, `pct_<custom>` | Custom nuisance category counts (if provided) |
| `n_nuisance`, `pct_nuisance` | Aggregate nuisance gene count and percentage |
| `condition` | Human-readable label combining layer + flavor + batch_key |

### Recommendation Output

`recommend_optimal()` returns a `pd.Series` containing all sweep columns plus `hvg_quality_score`. It also prints a formatted summary:

```
============================================================
  RECOMMENDED HVG PARAMETERS
============================================================
  Layer:         soupX_counts
  Flavor:        pearson_residuals
  Batch key:     donor
  n_top_genes:   3000
  ──────────────────────────────────
  Marker coverage:  82.5% (132/160)
  Nuisance genes:   2.3% (69 of 3000 HVGs)
  Quality score:    0.4858  (marker_w=0.6, nuisance_w=0.4)
============================================================
```

## Example Notebook

The included `02_DimReduction.ipynb` demonstrates the full workflow in the context of a human skeletal muscle / tendon single-cell analysis:

1. Loads a concatenated, QC-filtered AnnData object with multiple normalisation layers
2. Runs the HVG parameter sweep across 6 layer-flavor combinations × 7 n_top_genes values × 2 batch settings (84 total configurations)
3. Generates diagnostic plots (heatmaps, line plots, bar charts)
4. Recommends and applies the optimal HVG selection
5. Proceeds to PCA, Harmony integration, and cell-type annotation

## Dependencies

| Package | Minimum Version | Purpose |
|---|---|---|
| `scanpy` | 1.9+ | HVG selection, single-cell analysis |
| `anndata` | 0.8+ | Data structure |
| `numpy` | 1.21+ | Numerical operations |
| `pandas` | 1.4+ | DataFrames for sweep results |
| `matplotlib` | 3.5+ | Plotting |
| `seaborn` | 0.12+ | Heatmaps and statistical plots |

## License

MIT
