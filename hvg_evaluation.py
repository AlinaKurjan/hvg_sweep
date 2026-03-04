"""Systematic evaluation of highly variable gene (HVG) selection parameters.

Sweeps over batch_key, n_top_genes, layer, and HVG flavor, measuring marker
gene coverage and nuisance gene contamination against a reference marker
dictionary (e.g. COARSE_MARKERS).
"""

from __future__ import annotations

import re
import warnings
from itertools import product
from typing import Any, Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from anndata import AnnData

from marker_genes import COARSE_MARKERS, filter_markers_to_adata

# ── Layer-flavor compatibility ───────────────────────────────────────────────
# Default pairings when the caller does not specify explicit combinations.
DEFAULT_LAYER_FLAVORS: List[Tuple[str, str]] = [
    ("soupX_counts", "seurat_v3"),
    ("soupX_counts", "pearson_residuals"),
    ("log1p_norm", "seurat"),
    ("log1p_norm", "cell_ranger"),
    ("scran_log1p", "seurat"),
    ("scran_log1p", "cell_ranger"),
]

DEFAULT_N_TOP_GENES: List[int] = [1000, 1500, 2000, 2500, 3000, 4000, 5000]
DEFAULT_BATCH_KEYS: List[Optional[str]] = ["donor", None]

_COUNTS_FLAVORS = {"seurat_v3", "pearson_residuals"}
_NORM_FLAVORS = {"seurat", "cell_ranger"}

# ── Nuisance gene categories ────────────────────────────────────────────────
# Prefix-based patterns (fallback when adata.var lacks boolean flags).
_MITO_RE = re.compile(r"^MT-", re.IGNORECASE)
_RIBO_RE = re.compile(r"^RP[SL]\d", re.IGNORECASE)
_HB_RE = re.compile(r"^HB[ABDEQGMZ]\d?$", re.IGNORECASE)

_BUILTIN_NUISANCE = ("mito", "ribo", "hb")


# ── Nuisance gene helpers ────────────────────────────────────────────────────

def classify_nuisance_genes(
    adata: AnnData,
    custom_nuisance: Dict[str, List[str]] | None = None,
) -> Dict[str, Set[str]]:
    """Build sets of nuisance gene names.

    Always includes the three built-in categories (mitochondrial, ribosomal,
    hemoglobin), using boolean columns in ``adata.var`` when available and
    falling back to prefix-based regex matching otherwise.

    Additional user-defined categories can be supplied via
    ``custom_nuisance``; only genes actually present in ``adata.var_names``
    are retained.

    Parameters
    ----------
    adata
        Annotated data matrix.
    custom_nuisance
        Optional dict of ``{category_name: [gene1, gene2, ...]}`` for
        user-defined nuisance gene sets (e.g. tissue-specific contaminants).

    Returns
    -------
    dict
        ``{category: set_of_gene_names, ...}``
    """
    var = adata.var
    names = np.asarray(adata.var_names)
    var_set = set(names)

    result: Dict[str, Set[str]] = {}

    if "mt" in var.columns and var["mt"].dtype == bool:
        result["mito"] = set(names[var["mt"].values])
    else:
        result["mito"] = {g for g in names if _MITO_RE.match(g)}

    if "ribo" in var.columns and var["ribo"].dtype == bool:
        result["ribo"] = set(names[var["ribo"].values])
    else:
        result["ribo"] = {g for g in names if _RIBO_RE.match(g)}

    if "hb" in var.columns and var["hb"].dtype == bool:
        result["hb"] = set(names[var["hb"].values])
    else:
        result["hb"] = {g for g in names if _HB_RE.match(g)}

    if custom_nuisance:
        for cat_name, gene_list in custom_nuisance.items():
            result[cat_name] = {g for g in gene_list if g in var_set}

    return result


def compute_nuisance_counts(
    hvg_set: Set[str],
    nuisance_sets: Dict[str, Set[str]],
) -> Dict[str, Any]:
    """Count nuisance genes present in an HVG set.

    Iterates over all categories in *nuisance_sets* (both built-in and
    custom), so the output adapts to whatever categories are provided.

    Parameters
    ----------
    hvg_set
        Gene names selected as HVGs.
    nuisance_sets
        Output of :func:`classify_nuisance_genes`.

    Returns
    -------
    dict
        Per-category counts (``n_<cat>``, ``pct_<cat>``, ``genes_<cat>``),
        plus aggregate ``n_nuisance`` / ``pct_nuisance`` /
        ``nuisance_genes``.  The ``categories`` key lists category names
        in order.
    """
    n_hvg = len(hvg_set)
    out: Dict[str, Any] = {}
    total_nuisance: Set[str] = set()
    categories = list(nuisance_sets.keys())

    for cat in categories:
        overlap = hvg_set & nuisance_sets.get(cat, set())
        out[f"n_{cat}"] = len(overlap)
        out[f"pct_{cat}"] = len(overlap) / n_hvg * 100 if n_hvg else 0.0
        out[f"genes_{cat}"] = sorted(overlap)
        total_nuisance |= overlap

    out["n_nuisance"] = len(total_nuisance)
    out["pct_nuisance"] = len(total_nuisance) / n_hvg * 100 if n_hvg else 0.0
    out["nuisance_genes"] = sorted(total_nuisance)
    out["categories"] = categories
    return out


# ── Core functions ───────────────────────────────────────────────────────────

def select_hvgs(
    adata: AnnData,
    layer: str,
    flavor: str,
    n_top_genes: int,
    batch_key: Optional[str] = None,
) -> Set[str]:
    """Run HVG selection and return the set of selected gene names.

    Handles API differences between scanpy flavors:
      - ``seurat_v3``: ``sc.pp.highly_variable_genes`` with ``layer`` param.
      - ``pearson_residuals``: ``sc.experimental.pp.highly_variable_genes``.
      - ``seurat`` / ``cell_ranger``: operate on ``adata.X``; the function
        temporarily swaps X with the requested layer on a shallow copy.

    Parameters
    ----------
    adata
        Annotated data matrix (not mutated).
    layer
        Name of the layer to use for HVG computation.
    flavor
        One of ``'seurat_v3'``, ``'pearson_residuals'``, ``'seurat'``,
        ``'cell_ranger'``.
    n_top_genes
        Number of top HVGs to select.
    batch_key
        Obs column for batch-aware selection, or ``None`` for batch-naive.

    Returns
    -------
    set of str
        Gene names flagged as highly variable.
    """
    if layer not in adata.layers:
        raise KeyError(f"Layer '{layer}' not found in adata.layers")

    ad = adata.copy()

    if flavor == "pearson_residuals":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.experimental.pp.highly_variable_genes(
                ad,
                flavor="pearson_residuals",
                n_top_genes=n_top_genes,
                layer=layer,
                batch_key=batch_key,
            )
    elif flavor == "seurat_v3":
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(
                ad,
                flavor="seurat_v3",
                n_top_genes=n_top_genes,
                layer=layer,
                batch_key=batch_key,
            )
    elif flavor in ("seurat", "cell_ranger"):
        ad.X = ad.layers[layer].copy()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.pp.highly_variable_genes(
                ad,
                flavor=flavor,
                n_top_genes=n_top_genes,
                batch_key=batch_key,
            )
    else:
        raise ValueError(f"Unknown flavor '{flavor}'")

    return set(ad.var_names[ad.var["highly_variable"]])


def compute_marker_coverage(
    hvg_set: Set[str],
    marker_dict: Dict[str, List[str]],
    var_names,
) -> Dict[str, Any]:
    """Compute per-cell-type and overall marker coverage within an HVG set.

    Parameters
    ----------
    hvg_set
        Set of gene names selected as HVGs.
    marker_dict
        Cell type -> gene list mapping (e.g. ``COARSE_MARKERS``).
    var_names
        All gene names present in the AnnData (``adata.var_names``).

    Returns
    -------
    dict
        Keys: ``'per_type'`` (dict of cell-type -> {n_in_data, n_in_hvg,
        fraction}), ``'overall_n_in_data'``, ``'overall_n_in_hvg'``,
        ``'overall_fraction'``.
    """
    var_set = set(var_names)
    per_type: Dict[str, Dict[str, Any]] = {}
    all_in_data: Set[str] = set()
    all_in_hvg: Set[str] = set()

    for cell_type, genes in marker_dict.items():
        in_data = [g for g in genes if g in var_set]
        in_hvg = [g for g in in_data if g in hvg_set]
        frac = len(in_hvg) / len(in_data) if in_data else np.nan
        per_type[cell_type] = {
            "n_in_data": len(in_data),
            "n_in_hvg": len(in_hvg),
            "fraction": frac,
            "missing": sorted(set(in_data) - hvg_set),
        }
        all_in_data.update(in_data)
        all_in_hvg.update(g for g in in_data if g in hvg_set)

    overall_frac = len(all_in_hvg) / len(all_in_data) if all_in_data else np.nan
    return {
        "per_type": per_type,
        "overall_n_in_data": len(all_in_data),
        "overall_n_in_hvg": len(all_in_hvg),
        "overall_fraction": overall_frac,
    }


def run_hvg_sweep(
    adata: AnnData,
    marker_dict: Dict[str, List[str]] | None = None,
    layer_flavors: List[Tuple[str, str]] | None = None,
    n_top_genes_range: List[int] | None = None,
    batch_keys: List[Optional[str]] | None = None,
    custom_nuisance: Dict[str, List[str]] | None = None,
    *,
    verbose: bool = True,
) -> pd.DataFrame:
    """Run full parameter-grid HVG evaluation and return a tidy DataFrame.

    Parameters
    ----------
    adata
        Annotated data matrix.
    marker_dict
        Cell type -> gene list mapping.  Defaults to ``COARSE_MARKERS``.
    layer_flavors
        List of ``(layer, flavor)`` tuples.  Defaults to
        ``DEFAULT_LAYER_FLAVORS``.
    n_top_genes_range
        HVG counts to test.  Defaults to ``DEFAULT_N_TOP_GENES``.
    batch_keys
        Batch keys to test (``None`` = batch-naive).  Defaults to
        ``DEFAULT_BATCH_KEYS``.
    custom_nuisance
        Optional dict of ``{category_name: [gene1, gene2, ...]}`` for
        user-defined nuisance gene sets (e.g. tissue-specific
        contaminants like myofibre genes in tendon samples).  These are
        tracked alongside the built-in mito/ribo/hb categories.
    verbose
        Print progress.

    Returns
    -------
    pd.DataFrame
        One row per parameter combination.  Columns include ``layer``,
        ``flavor``, ``batch_key``, ``n_top_genes``, ``n_hvgs_actual``,
        ``overall_coverage``, ``overall_fraction``, per-category nuisance
        counts/percentages (``n_<cat>``, ``pct_<cat>``), aggregate
        ``n_nuisance`` / ``pct_nuisance``, and one ``frac_<CellType>``
        column per marker-dict entry.
    """
    if marker_dict is None:
        marker_dict = COARSE_MARKERS
    if layer_flavors is None:
        layer_flavors = DEFAULT_LAYER_FLAVORS
    if n_top_genes_range is None:
        n_top_genes_range = DEFAULT_N_TOP_GENES
    if batch_keys is None:
        batch_keys = DEFAULT_BATCH_KEYS

    filtered_markers = filter_markers_to_adata(marker_dict, adata.var_names)
    cell_types = list(filtered_markers.keys())

    nuisance_sets = classify_nuisance_genes(adata, custom_nuisance)
    nui_categories = list(nuisance_sets.keys())
    if verbose:
        for cat in nui_categories:
            print(f"  Nuisance pool — {cat}: {len(nuisance_sets[cat])} genes "
                  f"in dataset")

    available_layers = set(adata.layers.keys())
    layer_flavors = [
        (l, f) for l, f in layer_flavors if l in available_layers
    ]
    if not layer_flavors:
        raise ValueError(
            f"None of the requested layers are present. "
            f"Available: {sorted(available_layers)}"
        )

    grid = list(product(layer_flavors, batch_keys, n_top_genes_range))
    n_total = len(grid)

    rows: List[Dict[str, Any]] = []
    for i, ((layer, flavor), bk, n_top) in enumerate(grid, 1):
        label = f"{layer}+{flavor}, batch={bk}, n={n_top}"
        if verbose:
            print(f"  [{i:3d}/{n_total}] {label}", end="")

        try:
            hvgs = select_hvgs(adata, layer, flavor, n_top, batch_key=bk)
        except Exception as exc:
            if verbose:
                print(f"  FAILED: {exc}")
            continue

        cov = compute_marker_coverage(hvgs, filtered_markers, adata.var_names)
        nui = compute_nuisance_counts(hvgs, nuisance_sets)

        row: Dict[str, Any] = {
            "layer": layer,
            "flavor": flavor,
            "batch_key": bk if bk is not None else "none",
            "n_top_genes": n_top,
            "n_hvgs_actual": len(hvgs),
            "overall_n_in_data": cov["overall_n_in_data"],
            "overall_n_in_hvg": cov["overall_n_in_hvg"],
            "overall_fraction": cov["overall_fraction"],
        }
        for cat in nui_categories:
            row[f"n_{cat}"] = nui[f"n_{cat}"]
            row[f"pct_{cat}"] = nui[f"pct_{cat}"]
        row["n_nuisance"] = nui["n_nuisance"]
        row["pct_nuisance"] = nui["pct_nuisance"]

        for ct in cell_types:
            row[f"frac_{ct}"] = cov["per_type"][ct]["fraction"]
        rows.append(row)

        if verbose:
            per_cat = ", ".join(
                f"{cat}={nui[f'n_{cat}']}" for cat in nui_categories
            )
            print(f"  -> markers {cov['overall_n_in_hvg']}/"
                  f"{cov['overall_n_in_data']} "
                  f"({cov['overall_fraction']:.1%}), "
                  f"nuisance {nui['n_nuisance']} ({per_cat})")

    df = pd.DataFrame(rows)
    df["condition"] = (
        df["layer"] + " + " + df["flavor"]
        + " | batch=" + df["batch_key"].astype(str)
    )
    return df


# ── Plotting ─────────────────────────────────────────────────────────────────

def plot_marker_coverage_heatmap(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (14, 8),
    cmap: str = "YlGnBu",
    annot: bool = True,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of overall marker coverage: conditions (rows) x n_top_genes (cols).

    Parameters
    ----------
    results_df
        Output of :func:`run_hvg_sweep`.
    """
    pivot = results_df.pivot_table(
        index="condition",
        columns="n_top_genes",
        values="overall_fraction",
        aggfunc="first",
    )
    pivot = pivot.sort_index()

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot,
        annot=annot,
        fmt=".2f",
        cmap=cmap,
        vmin=0,
        vmax=1,
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title("Overall marker coverage fraction across HVG parameter grid")
    ax.set_xlabel("n_top_genes")
    ax.set_ylabel("")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_marker_coverage_lines(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Line plot of overall marker coverage vs n_top_genes, one line per condition."""
    fig, ax = plt.subplots(figsize=figsize)

    for cond, grp in results_df.groupby("condition"):
        grp_sorted = grp.sort_values("n_top_genes")
        ax.plot(
            grp_sorted["n_top_genes"],
            grp_sorted["overall_fraction"],
            marker="o",
            label=cond,
        )

    ax.set_xlabel("n_top_genes")
    ax.set_ylabel("Overall marker coverage fraction")
    ax.set_title("Marker coverage vs HVG count")
    ax.set_ylim(0, 1.05)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_per_celltype_coverage(
    results_df: pd.DataFrame,
    n_top_genes: int = 3000,
    figsize: Tuple[float, float] | None = None,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Grouped bar chart of per-cell-type marker coverage for a given n_top_genes."""
    subset = results_df[results_df["n_top_genes"] == n_top_genes].copy()
    if subset.empty:
        raise ValueError(f"No results for n_top_genes={n_top_genes}")

    frac_cols = [c for c in subset.columns if c.startswith("frac_")]
    ct_names = [c.replace("frac_", "") for c in frac_cols]

    melted = subset.melt(
        id_vars=["condition"],
        value_vars=frac_cols,
        var_name="cell_type_col",
        value_name="fraction",
    )
    melted["cell_type"] = melted["cell_type_col"].str.replace("frac_", "", regex=False)

    n_ct = len(ct_names)
    n_cond = melted["condition"].nunique()
    if figsize is None:
        figsize = (max(14, n_ct * 0.8), max(6, n_cond * 0.6 + 2))

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(
        data=melted,
        x="cell_type",
        y="fraction",
        hue="condition",
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Marker coverage fraction")
    ax.set_title(f"Per-cell-type marker coverage (n_top_genes={n_top_genes})")
    ax.set_ylim(0, 1.05)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_nuisance_heatmap(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (14, 8),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Heatmap of nuisance gene percentage: conditions (rows) x n_top_genes (cols).

    Shows the fraction of each HVG set occupied by mitochondrial, ribosomal,
    and hemoglobin genes — slots that are typically uninformative for cell-type
    discrimination.

    Parameters
    ----------
    results_df
        Output of :func:`run_hvg_sweep`.
    """
    pivot = results_df.pivot_table(
        index="condition",
        columns="n_top_genes",
        values="pct_nuisance",
        aggfunc="first",
    )
    pivot = pivot.sort_index()

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        pivot,
        annot=True,
        fmt=".1f",
        cmap="YlOrRd",
        vmin=0,
        linewidths=0.5,
        ax=ax,
    )
    ax.set_title("Nuisance genes in HVGs (% of total HVGs): mito + ribo + hb")
    ax.set_xlabel("n_top_genes")
    ax.set_ylabel("")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_nuisance_breakdown(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (12, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Stacked bar chart of nuisance gene counts per condition at each n_top_genes.

    Each bar is split into per-category segments.  Categories are discovered
    dynamically from ``n_*`` columns in the DataFrame (excluding
    ``n_nuisance``, ``n_top_genes``, ``n_hvgs_actual``).

    Parameters
    ----------
    results_df
        Output of :func:`run_hvg_sweep`.
    """
    skip = {"n_nuisance", "n_top_genes", "n_hvgs_actual",
            "overall_n_in_data", "overall_n_in_hvg"}
    n_cols = [c for c in results_df.columns
              if c.startswith("n_") and c not in skip]
    cat_names = [c.replace("n_", "", 1) for c in n_cols]

    _BUILTIN_LABELS = {
        "mito": "Mitochondrial",
        "ribo": "Ribosomal",
        "hb": "Hemoglobin",
    }
    label_map = {c: _BUILTIN_LABELS.get(c, c) for c in cat_names}

    melted = results_df.melt(
        id_vars=["condition", "n_top_genes"],
        value_vars=n_cols,
        var_name="category_col",
        value_name="count",
    )
    melted["category"] = (
        melted["category_col"]
        .str.replace("n_", "", 1, regex=False)
        .map(label_map)
    )

    n_conditions = results_df["condition"].nunique()
    n_nvhg = results_df["n_top_genes"].nunique()
    if n_conditions * n_nvhg > 30:
        figsize = (max(figsize[0], n_conditions * n_nvhg * 0.25), figsize[1])

    fig, ax = plt.subplots(figsize=figsize)
    melted["x_label"] = (
        melted["condition"] + "\nn=" + melted["n_top_genes"].astype(str)
    )

    _BUILTIN_COLORS = {
        "Mitochondrial": "#e74c3c",
        "Ribosomal": "#f39c12",
        "Hemoglobin": "#8e44ad",
    }
    unique_labels = melted["category"].unique()
    tab10 = plt.cm.tab10.colors
    color_idx = 0
    palette = {}
    for lab in unique_labels:
        if lab in _BUILTIN_COLORS:
            palette[lab] = _BUILTIN_COLORS[lab]
        else:
            palette[lab] = tab10[color_idx % len(tab10)]
            color_idx += 1

    conditions_sorted = sorted(results_df["condition"].unique())
    n_top_sorted = sorted(results_df["n_top_genes"].unique())
    x_order = [
        f"{c}\nn={n}" for c in conditions_sorted for n in n_top_sorted
    ]

    sns.barplot(
        data=melted,
        x="x_label",
        y="count",
        hue="category",
        order=x_order,
        palette=palette,
        ax=ax,
    )
    ax.set_xlabel("")
    ax.set_ylabel("Number of nuisance genes in HVGs")
    ax.set_title("Nuisance gene breakdown in HVG sets")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="center",
                       fontsize=6)
    ax.legend(title="Category")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_nuisance_lines(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (10, 6),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Line plot of nuisance gene percentage vs n_top_genes, one line per condition."""
    fig, ax = plt.subplots(figsize=figsize)

    for cond, grp in results_df.groupby("condition"):
        grp_sorted = grp.sort_values("n_top_genes")
        ax.plot(
            grp_sorted["n_top_genes"],
            grp_sorted["pct_nuisance"],
            marker="o",
            label=cond,
        )

    ax.set_xlabel("n_top_genes")
    ax.set_ylabel("Nuisance genes (% of HVGs)")
    ax.set_title("Nuisance gene contamination vs HVG count")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def plot_quality_summary(
    results_df: pd.DataFrame,
    figsize: Tuple[float, float] = (12, 7),
    save_path: Optional[str] = None,
) -> plt.Figure:
    """Two-panel summary: marker coverage vs nuisance contamination.

    Left panel shows overall marker fraction (higher is better); right panel
    shows nuisance percentage (lower is better).  Enables quick visual
    comparison of which parameter configurations strike the best balance.

    Parameters
    ----------
    results_df
        Output of :func:`run_hvg_sweep`.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, sharey=True)

    for cond, grp in results_df.groupby("condition"):
        grp_sorted = grp.sort_values("n_top_genes")
        ax1.plot(grp_sorted["n_top_genes"], grp_sorted["overall_fraction"],
                 marker="o", label=cond)
        ax2.plot(grp_sorted["n_top_genes"], grp_sorted["pct_nuisance"],
                 marker="s", label=cond)

    ax1.set_xlabel("n_top_genes")
    ax1.set_ylabel("Fraction / Percentage")
    ax1.set_title("Marker coverage (higher = better)")
    ax1.set_ylim(bottom=0)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel("n_top_genes")
    ax2.set_title("Nuisance genes % (lower = better)")
    ax2.grid(True, alpha=0.3)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=7)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
    return fig


def recommend_optimal(
    results_df: pd.DataFrame,
    marker_weight: float = 0.6,
    nuisance_weight: float = 0.4,
    min_marker_coverage: float = 0.0,
    max_nuisance_pct: float = 100.0,
    verbose: bool = True,
) -> pd.Series:
    """Score every parameter combination and return the best row.

    Computes a composite quality score for each configuration::

        score = marker_weight * overall_fraction
              - nuisance_weight * (pct_nuisance / 100)

    Higher is better (more marker coverage, less nuisance contamination).
    Optionally filter rows that fail hard thresholds before scoring.

    Parameters
    ----------
    results_df
        Output of :func:`run_hvg_sweep`.
    marker_weight
        Weight given to marker coverage fraction (0–1 scale) in the score.
    nuisance_weight
        Weight given to nuisance percentage (0–1 scale) in the score.
    min_marker_coverage
        Exclude configurations with ``overall_fraction`` below this value
        (e.g. 0.7 = require at least 70 % marker coverage).
    max_nuisance_pct
        Exclude configurations with ``pct_nuisance`` above this value
        (e.g. 5.0 = reject if >5 % nuisance).
    verbose
        Print a summary of the recommendation.

    Returns
    -------
    pd.Series
        The row from *results_df* with the highest composite score,
        augmented with an ``hvg_quality_score`` field.
    """
    df = results_df.copy()
    df = df[df["overall_fraction"] >= min_marker_coverage]
    df = df[df["pct_nuisance"] <= max_nuisance_pct]
    if df.empty:
        raise ValueError(
            "No configurations survive the thresholds "
            f"(min_marker_coverage={min_marker_coverage}, "
            f"max_nuisance_pct={max_nuisance_pct}). "
            "Try relaxing the constraints."
        )

    df["hvg_quality_score"] = (
        marker_weight * df["overall_fraction"]
        - nuisance_weight * (df["pct_nuisance"] / 100)
    )

    best_idx = df["hvg_quality_score"].idxmax()
    best = df.loc[best_idx].copy()

    if verbose:
        print("=" * 60)
        print("  RECOMMENDED HVG PARAMETERS")
        print("=" * 60)
        print(f"  Layer:         {best['layer']}")
        print(f"  Flavor:        {best['flavor']}")
        bk = best["batch_key"]
        print(f"  Batch key:     {bk if bk != 'none' else 'None (batch-naive)'}")
        print(f"  n_top_genes:   {int(best['n_top_genes'])}")
        print(f"  ──────────────────────────────────")
        print(f"  Marker coverage:  {best['overall_fraction']:.1%} "
              f"({int(best['overall_n_in_hvg'])}/{int(best['overall_n_in_data'])})")
        print(f"  Nuisance genes:   {best['pct_nuisance']:.1f}% "
              f"({int(best['n_nuisance'])} of {int(best['n_hvgs_actual'])} HVGs)")
        print(f"  Quality score:    {best['hvg_quality_score']:.4f}  "
              f"(marker_w={marker_weight}, nuisance_w={nuisance_weight})")
        print("=" * 60)

    return best


def apply_optimal_hvgs(
    adata: AnnData,
    best: pd.Series,
) -> Set[str]:
    """Apply HVG selection to *adata* using the recommended parameters.

    Calls :func:`select_hvgs` with the layer / flavor / n_top_genes /
    batch_key from the recommendation, then writes the result into
    ``adata.var['highly_variable']``.

    Parameters
    ----------
    adata
        Annotated data matrix (**modified in place**).
    best
        Output of :func:`recommend_optimal`.

    Returns
    -------
    set of str
        The selected HVG names.
    """
    layer = best["layer"]
    flavor = best["flavor"]
    n_top = int(best["n_top_genes"])
    bk = best["batch_key"]
    batch_key = bk if bk != "none" else None

    hvgs = select_hvgs(adata, layer, flavor, n_top, batch_key=batch_key)

    adata.var["highly_variable"] = adata.var_names.isin(hvgs)
    print(f"Applied {len(hvgs)} HVGs to adata.var['highly_variable']  "
          f"(layer={layer}, flavor={flavor}, "
          f"n_top_genes={n_top}, batch_key={batch_key})")
    return hvgs


def plot_missing_markers_table(
    adata: AnnData,
    hvg_set: Set[str],
    marker_dict: Dict[str, List[str]] | None = None,
) -> pd.DataFrame:
    """Return a DataFrame showing which markers are missing from a given HVG set.

    Useful for inspecting a single chosen configuration in detail.

    Parameters
    ----------
    adata
        Annotated data matrix.
    hvg_set
        Set of HVG gene names.
    marker_dict
        Defaults to ``COARSE_MARKERS``.

    Returns
    -------
    pd.DataFrame
        Columns: cell_type, n_total, n_in_data, n_in_hvg, fraction, missing_genes.
    """
    if marker_dict is None:
        marker_dict = COARSE_MARKERS

    var_set = set(adata.var_names)
    rows = []
    for ct, genes in marker_dict.items():
        in_data = [g for g in genes if g in var_set]
        in_hvg = [g for g in in_data if g in hvg_set]
        missing = sorted(set(in_data) - hvg_set)
        rows.append({
            "cell_type": ct,
            "n_total": len(genes),
            "n_in_data": len(in_data),
            "n_in_hvg": len(in_hvg),
            "fraction": len(in_hvg) / len(in_data) if in_data else np.nan,
            "missing_genes": ", ".join(missing) if missing else "",
        })
    return pd.DataFrame(rows)
