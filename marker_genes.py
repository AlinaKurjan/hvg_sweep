"""Marker gene templates for HVG evaluation.

This module provides **example** marker dictionaries for use with
``hvg_evaluation.run_hvg_sweep``.  The defaults are deliberately simple and
fairly generic (broad immune / stromal / epithelial compartments) so that
they are:

- Easy to understand
- Easy to replace with your own biology

You will get the **best results** if you customise both:

1. ``COARSE_MARKERS`` – cell-type–specific markers for your tissue / organism
2. The nuisance gene sets passed as ``custom_nuisance`` to
   :func:`hvg_evaluation.run_hvg_sweep`

Practical tips for **choosing positive markers** (for ``COARSE_MARKERS``)
-----------------------------------------------------------------------
- Prefer genes that are **strongly enriched in one cell type** and low in
  most others (e.g. ``CD3D`` for T cells, ``PECAM1`` for endothelial cells).
- Use markers that are **robust across datasets** (classical, well-validated
  lineage markers from papers / atlases / CellMarker, etc.).
- Aim for **8–25 genes per coarse type** – enough to be stable, but not so
  many that every gene is only weakly specific.
- Avoid:
  - Housekeeping genes (e.g. ``ACTB``, ``GAPDH``, ribosomal ``RPL*``/``RPS*``)
  - Ubiquitous stress / immediate-early genes (e.g. ``FOS``, ``JUN``, ``HSPA1A``)
  - Cell-cycle genes (e.g. ``MKI67``, ``TOP2A``) – treat these as nuisance
  - Extremely tissue-specific structural genes that dominate counts but are
    not helpful for distinguishing *within-tissue* cell types (e.g. myofibre
    contractile genes in mixed muscle samples – good candidates for a
    nuisance category instead).

Practical tips for **choosing nuisance markers**
-----------------------------------------------
- The evaluation code already classifies **mitochondrial, ribosomal and
  hemoglobin genes automatically**, using prefixes (``MT-``, ``RPL*``,
  ``RPS*``, ``HB*``) or boolean flags in ``adata.var``.
- You can add **extra nuisance categories** via the ``custom_nuisance`` dict
  argument to :func:`run_hvg_sweep`, for example:

  >>> custom_nuisance = {
  ...     "tissue_structural": ["TTN", "NEB", "MYH7"],  # e.g., if you don't expect muscle cells to be present
  ...     "ambient_epithelium": ["KRT5", "KRT14", "KRT17"],  # e.g., skin / mucosa carry-over into other tissues
  ... }

- Good nuisance candidates are **highly abundant but biologically unhelpful
  for your question**, e.g.:
  - Ambient epithelial keratins when your tissue is mostly stroma
  - Strong lineage-structural genes you do not want to drive clustering
  - Immunoglobulin genes in non-B-cell–focused analyses
  - etc.

Gene symbols below use **human HGNC** convention; if you are working with
mouse or another organism, adjust to the naming used in ``adata.var_names``.
"""

from typing import Dict, List


# ===========================================================================
# COARSE MARKERS (EXAMPLE, EDIT FOR YOUR DATA)
#
# These broad compartments are intentionally generic and should be treated as
# a starting template.  For serious analyses, copy this file into your own
# project and replace / extend the entries with markers that fit your tissue.
# ===========================================================================

COARSE_MARKERS: Dict[str, List[str]] = {
    # Epithelial & barrier compartments
    "Epithelial": [
        "EPCAM",
        "KRT8",
        "KRT18",
        "KRT19",
        "KRT17",
        "CDH1",
        "CLDN4",
        "CLDN7",
        "MUC1",
    ],

    # Generic fibroblasts / stromal cells
    "Fibroblasts": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "PDGFRA",
        "FAP",
        "THY1",
        "VIM",
        "COL6A1",
        "COL6A3",
    ],

    # Vascular endothelium
    "Endothelial": [
        "PECAM1",
        "CDH5",
        "VWF",
        "ERG",
        "FLT1",
        "KDR",
        "CLDN5",
        "EMCN",
        "ESAM",
    ],

    # Pericytes / vascular smooth muscle
    "Pericytes_vSMC": [
        "PDGFRB",
        "RGS5",
        "ABCC9",
        "KCNJ8",
        "NOTCH3",
        "MCAM",
        "ACTA2",
        "TAGLN",
        "MYH11",
    ],

    # Myeloid lineage (monocytes, macrophages, dendritic cells)
    "Myeloid": [
        "LYZ",
        "CSF1R",
        "CD68",
        "AIF1",
        "ITGAM",
        "LST1",
        "C1QA",
        "C1QB",
        "C1QC",
        "FCGR3A",
        "S100A8",
        "S100A9",
    ],

    # T cells (CD4+, CD8+, Tregs)
    "T_Cells": [
        "CD3D",
        "CD3E",
        "CD3G",
        "TRAC",
        "IL7R",
        "CD4",
        "CD8A",
        "CD8B",
        "CCR7",
        "GZMK",
        "CCL5",
        "FOXP3",
    ],

    # NK cells and cytotoxic lymphocytes
    "NK_Cells": [
        "NKG7",
        "GNLY",
        "KLRD1",
        "NCAM1",
        "PRF1",
        "GZMB",
        "KLRF1",
    ],

    # B cells and plasma cells
    "B_Cells": [
        "MS4A1",
        "CD79A",
        "CD79B",
        "CD19",
        "CD74",
        "BANK1",
        "PAX5",
        "TCL1A",
    ],
    "Plasma_Cells": [
        "MZB1",
        "JCHAIN",
        "XBP1",
        "SDC1",
        "IGHG1",
        "IGKC",
    ],

    # Mast cells
    "Mast_Cells": [
        "KIT",
        "TPSAB1",
        "TPSB2",
        "CPA3",
        "HPGDS",
        "MS4A2",
    ],

    # Adipocytes
    "Adipocytes": [
        "PLIN1",
        "ADIPOQ",
        "LEP",
        "FABP4",
        "LPL",
        "CIDEC",
        "PPARG",
    ],

    # Neural / glial (peripheral or central, depending on context)
    "Neural_Glia": [
        "PLP1",
        "MBP",
        "SOX10",
        "S100B",
        "GFAP",
        "SNAP25",
        "SYT1",
        "TUBB3",
    ],

    # Erythroid cells and platelets
    "Erythrocytes": [
        "HBA1",
        "HBA2",
        "HBB",
        "ALAS2",
        "SLC4A1",
    ],
    "Platelets": [
        "PF4",
        "PPBP",
        "ITGA2B",
        "GP9",
        "TUBB1",
        "ITGB3",
    ],
}

# ===========================================================================
# Helper functions
# ===========================================================================

def filter_markers_to_adata(
    marker_dict: Dict[str, List[str]],
    var_names,
) -> Dict[str, List[str]]:
    """Filter marker dictionary to only include genes present in the dataset.

    Parameters
    ----------
    marker_dict : dict
        Cell type -> gene list mapping.
    var_names : Index or list
        Gene names present in the AnnData object.

    Returns
    -------
    dict
        Filtered marker dictionary (types with <2 genes removed).
    """
    var_set = set(var_names)
    filtered: Dict[str, List[str]] = {}
    for cell_type, genes in marker_dict.items():
        valid = [g for g in genes if g in var_set]
        if len(valid) >= 2:
            filtered[cell_type] = valid
    return filtered
