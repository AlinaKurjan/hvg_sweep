"""Canonical marker gene dictionary for human skeletal muscle cell types.

Single source of truth used across all within-study and cross-study analyses.

Annotation levels:
  - COARSE_MARKERS: ~18 major cell compartments
  - FINE_MARKERS: ~58 subtypes within compartments
  - NEGATIVE_MARKERS: genes that should NOT be expressed (for disambiguation)

Design principles for score_genes() compatibility:
  - Each list has 4-25 genes: enough for stable scoring, not so many that
    specificity is diluted.
  - Lists at the same level (coarse or fine) minimise overlap — shared genes
    weaken the ability of sc.tl.score_genes to discriminate.
  - Fine-level lists are designed to discriminate *within* the parent coarse
    compartment, not across all cell types.
  - Shared markers are flagged in comments rather than silently removed,
    because biological co-expression is real; disambiguation relies on the
    full gene-list context, not individual genes.

Sources
-------
Muscle literature:
  Rubenstein et al. 2020 Cell Reports — human VL scRNA-seq
  Perez et al. 2022 Aging Cell — human VL snRNA-seq aging
  De Micheli et al. 2020 Skeletal Muscle — human multi-tissue scRNA-seq
  Dos Santos et al. 2020 Nature Cell Biology — satellite cell states
  Charville et al. 2015 Stem Cell Reports — human satellite cell markers
  McKellar et al. 2021 eLife — muscle regeneration atlas
  Murgia et al. 2015 PNAS — fiber type proteomics
  Petrany & Bhatt et al. 2020 Nature Communications — mouse snRNA-seq fiber types
  Tabula Sapiens 2022 Science — skeletal muscle tissue
  Chakarov et al. 2019 Science — tissue macrophage subsets
  Kalucka et al. 2020 Cell Reports — organ-resolved EC atlas
  Human Protein Atlas (HPA v23) — tissue-enriched genes and specificity ratios
  CellMarker 2.0 — MSK-specific cell type marker database (635 entries)

Tendon / connective tissue literature:
  Harvey et al. 2019 Nature Communications — human tendon scRNA-seq
  Kendal et al. 2020 Nature Scientific Reports — human tendon multi-omic CITE-seq
  Zhang et al. 2023 eLife — human enthesis (bone-tendon junction) scRNA-seq
  Mimpen et al. 2024 FASEB J (PMID: 38742770) — human hamstring tendon snRNA-seq
    (10,533 nuclei; 4 healthy donors; spatial transcriptomics co-registered)
  Mimpen et al. 2025 J Physiol (PMID: 40232153) — ruptured human quadriceps
    tendon scRNA-seq (12,808 nuclei; healthy vs ruptured comparison)
  Mimpen et al. 2025 bioRxiv 2025.06.19.660575 — human tendon across the
    lifespan (spatial + snRNA-seq; Achilles and quadriceps; embryonic to adult)
"""

from typing import Dict, List


# ===========================================================================
# COARSE MARKERS — ~18 broad compartments
#
# 8-25 genes per type.  Minimal cross-type overlap.
# Used as the first pass in hierarchical scoring.
# ===========================================================================

COARSE_MARKERS: Dict[str, List[str]] = {

    "Myonuclei_Specialised": [
        # NMJ (neuromuscular junction) myonuclei
        "CHRNE", "CHRNA1", "COLQ", "COL25A1", "MUSK", "RAPSN", "LRP4",
        "CPEB3", "CUX1",  # scANVI-validated (CPEB3hi CUX1hi NMJ Myocytes)
        # COL13A1: top NMJ marker upstream of COLQ (Petrany 2020 mouse snRNA-seq)
        # DOK7: NMJ signaling scaffold
        "COL13A1", "DOK7",
        # MTJ (myotendinous junction) myonuclei
        "COL22A1", "PIEZO2", "COL24A1",
        # Immature / regenerating
        # NRAP, DYSF: transitional / unidentified myonuclei (Petrany 2020)
        "MYH3", "MYH8", "XIRP1", "MYMX", "NRAP", "DYSF",
        # Note: FLNC removed from coarse (broadly expressed in general myofibers,
        # not specific to specialised nuclei; retained in Regenerative_Myonuclei fine)
    ],

    "FAPs_Fibroblasts": [
        # Core FAP identity
        "PDGFRA", "THY1", "FAP", "HIC1", "NT5E",
        # Universal fibroblast ECM
        "DCN", "LUM", "COL1A1", "COL1A2", "COL3A1",
        # FAP-enriched signalling
        "PI16", "DPP4", "SFRP2", "SFRP4", "MFAP5", "GSN",
        "CFD", "FBLN1",
        # CD34: co-expressed with PDGFRA on FAPs (validated in multiple studies)
        # CD248 (Endosialin): tendon FAP marker (Harvey 2019)
        "CD34", "CD248",
        # Note: NGFR removed (belongs with Schwann_Cells_NonMyelinating)
        # Note: SMOC2, ANGPTL7 removed (less well-validated across datasets)
    ],

    "Endothelial_Cells": [
        # Pan-endothelial
        "PECAM1", "CDH5", "VWF", "ERG", "FLT1", "KDR",
        "TIE1", "ENG", "CLDN5", "EMCN",
        "ESAM", "ICAM2",
        # CLEC14A: validated in CellMarker seq data for muscle endothelial cells
        "CLEC14A",
        # Note: CD34 not added here despite EC expression — shared with FAPs
    ],

    "Mural_Cells": [
        # Pan-mural / pericyte
        "PDGFRB", "RGS5", "ABCC9", "KCNJ8", "NOTCH3", "HIGD1B",
        # CSPG4 (NG2): classical pericyte marker validated in CellMarker
        # MCAM (CD146): validated pericyte/mural marker in CellMarker
        "CSPG4", "MCAM",
        # vSMC-enriched
        "MYH11", "ACTA2", "CNN1", "TAGLN",
    ],

    "Macrophages_Monocytes": [
        # Pan-macrophage / myeloid
        "CD68", "CSF1R",
        # AIF1 (Iba1): pan-myeloid marker confirmed in CellMarker muscle data
        "AIF1",
        # Tissue-resident macrophage
        "CD163", "MRC1", "C1QA", "C1QB", "C1QC",
        "F13A1", "FOLR2", "MAF",
        # Monocyte
        "CD14", "FCGR3A", "LYZ",
        # Inflammatory (shared)
        "S100A8", "S100A9",
    ],

    "T_Cells": [
        # Pan-T (highest specificity markers first)
        "CD3D", "CD3E", "CD3G", "TRAC",
        # CD4 / CD8 lineage
        "IL7R", "CD4", "CD8A", "CD8B",
        # Effector
        "GZMA", "GZMK", "CCL5",
        # Regulatory
        "FOXP3",
        # Tissue residency
        "CD69",
    ],

    "NK_Cells": [
        "NKG7", "GNLY", "KLRD1", "NCAM1",
        "PRF1", "GZMB",
        "KLRF1", "KLRB1",
        "XCL1", "XCL2",
        "FCER1G",
    ],

    "Mast_Cells": [
        "KIT", "CPA3", "TPSAB1", "TPSB2", "MS4A2",
        "HDC", "GATA2", "HPGDS",
        # Note: HPGD removed (expressed in fibroblasts and other types;
        # HPGDS is the more specific mast cell marker)
    ],

    "B_Cells": [
        "CD79A", "CD79B", "MS4A1", "CD19",
        "TCL1A", "PAX5", "BANK1", "BLK",
    ],

    "Neural": [
        # Schwann cells (myelinating)
        "PLP1", "MPZ", "MBP", "PRX",
        # Schwann (general / non-myelinating)
        # NGFR: moved here from FAPs_Fibroblasts — more specific to Schwann cells
        "S100B", "CDH19", "SOX10", "NFASC", "NGFR",
        # Neurons / nerve fibres — note: rarely captured in standard muscle
        # scRNA-seq (low sensitivity in these tissue types)
        "SNAP25", "NRXN1", "SYT1", "NEFL",
    ],

    "Tendon_Ligament_Cells": [
        # Tenogenic transcription factors and structural ECM
        "SCX", "TNMD", "MKX", "FMOD", "THBS4", "KERA",
        "COL12A1", "TNC", "BGN",
        # Sheath / lubricating (Kendal 2020 Tenocyte E, Mimpen 2025 lifespan)
        "PRG4",
        # Fibrillar lineage (Mimpen 2025 lifespan bioRxiv: ABI3BP GAS2 fibroblasts)
        # ABI3BP: collagen cross-linking/fibrillar organisation
        # SPARC: matrix organisation; highly enriched in KERA/TNMD tenocytes
        "ABI3BP", "SPARC",
        # ECM remodelling (Mimpen 2025 J Physiol ruptured tendon)
        "CILP",
        # Note: TPPP3 not included at coarse (shared with Fibroblasts_PTGDShi)
        # Note: COMP not included at coarse (too broadly expressed in MSK ECM)
    ],

    "Chondrocytes": [
        # Definitive chondrocyte markers
        "SOX9", "COL2A1", "ACAN", "COL9A1", "COL11A1",
        "COL9A3", "CLEC3A",
        # COMP also in tendon ECM; context disambiguates
        "COMP",
    ],

    "Adipocytes": [
        "PLIN1", "ADIPOQ", "LEP", "FABP4",
        "GPAM", "PDE3B", "CIDEC", "PCK1",
        "LPL", "PPARG",
    ],

    "Mesothelial_Cells": [
        "MSLN", "WT1", "UPK3B", "KRT19",
        "CALB2", "ITLN1", "LRRN4",
    ],
}


# ===========================================================================
# FINE MARKERS — ~58 subtypes
#
# 4-12 genes per subtype.  Focused on discriminating *within* the parent
# coarse compartment rather than across all cell types.
# ===========================================================================

FINE_MARKERS: Dict[str, List[str]] = {

    # --- Myofiber subtypes (parent: Myofibers) -------------------------

    "Myofibers_type_I": [
        # Definitive slow myosin
        "MYH7", "MYH7B",
        # Slow troponins
        "TNNT1", "TNNC1", "TNNI1",
        # Slow calcium handling
        "ATP2A2", "CASQ2",
        # Type I-enriched
        "MYL3", "MYOZ2",
        "CA3",  # scANVI-validated (TNNI2low CA3hi Type I)
        # MYL2 (slow myosin light chain, Petrany 2020 top Type I marker)
        # TPM3 (slow tropomyosin, Petrany 2020)
        "MYL2", "TPM3",
    ],

    "Myofibers_type_IIA": [
        # Definitive fast-oxidative myosin
        "MYH2",
        # Fast troponins (shared with IIX)
        "TNNT3", "TNNC2", "TNNI2",
        # IIA-enriched (more oxidative)
        "ANKRD2", "MYOM3", "CASQ1",
        # MYBPC1: enriched in oxidative fibers (Petrany 2020)
        "MYBPC1",
        # Note: CD36 also enriched in oxidative fibers but shared with ECs and
        # Macrophages_LAM; not included here
    ],

    "Myofibers_type_IIX": [
        # Definitive fast-glycolytic myosin
        "MYH1",
        # Fast troponins (shared with IIA)
        "TNNT3", "TNNC2", "TNNI2",
        # IIX-enriched (more glycolytic)
        # MYBPC2: moved from coarse to here (IIX-biased)
        "MYLK2", "ACTN3", "MYBPC2", "ATP2A1",
        # MYL1: fast myosin light chain enriched in IIB/IIX clusters (Petrany 2020)
        "MYL1",
    ],

    # --- Specialised myonuclei (parent: Myonuclei_Specialised) ---------

    "NMJ_Myonuclei": [
        "CHRNE", "CHRNA1", "COLQ", "COL25A1",
        "MUSK", "RAPSN", "LRP4", "DOK7",
        "PRKAR1A", "ABLIM2", "UTRN",
        "CPEB3", "CUX1",  # scANVI-validated
        # COL13A1 (top NMJ marker, Petrany 2020)
        # KCNQ5, GRAMD1B (Petrany 2020 NMJ cluster)
        "COL13A1", "KCNQ5", "GRAMD1B",
        "ETV5",
    ],

    "MTJ_Myonuclei": [
        "COL22A1", "PIEZO2", "COL24A1",
        "COL6A1", "COL6A3", "FSTL1",
        # LAMA2, LIMCH1, PRKG1, FRAS1 (Petrany 2020 MTJ cluster)
        "LAMA2", "LIMCH1", "PRKG1", "FRAS1",
    ],

    "Regenerative_Myonuclei": [
        "MYH3", "MYH8", "XIRP1", "FLNC",
        # NRAP, DYSF: transitional / unidentified myonuclei (Petrany 2020)
        "NRAP", "DYSF",
        # Note: NCAM1 removed (strong NK cell marker — shared at fine level
        # causes NK ↔ regenerating myofiber confusion; see NK_Cells)
    ],

    "Immature_Myonuclei": [
        "MYMX", "MYOG", "MYF6", "PCNA", "CDK1",
    ],

    # --- Satellite cell subtypes (parent: Satellite_Cells) -------------

    "Satellite_Cells_Quiescent": [
        "PAX7", "MYF5", "SPRY1", "CALCR",
        "CHODL", "HEYL",
        "CD82", "FGFR4", "ITGA7", "DLK1",
        # De Micheli 2020 MuSC1 surface markers
        "EGFR", "CD151", "CD81", "SDC2",
        # Note: HEY1 removed at fine level (also a top Endothelial_Arterial
        # Notch marker; retained in Satellite_Cells coarse where disambiguation
        # from ECs works via the broader gene list)
    ],

    "Satellite_Cells_Activated": [
        # Core activation TFs
        "MYOD1", "MYOG", "MYF6",
        # Proliferation
        "CDK1", "PCNA", "MKI67",
        # Surface
        "CXCR4", "SDC4",
        # De Micheli 2020 MuSC2 activated markers
        "CD44", "TNFRSF12A",
    ],

    # --- FAP / Fibroblast subtypes (parent: FAPs_Fibroblasts) ----------

    "FAPs": [
        "PDGFRA", "THY1", "FAP", "HIC1", "NT5E",
        "CD34", "PI16", "DPP4",
        # CD248 (Endosialin): also marks tendon FAPs (Harvey 2019)
        "CD248",
    ],

    "FAPs_Adipogenic": [
        "EBF1", "TSHZ2", "LAMA2", "ABCA9", "RBMS3",
        "ABCA10", "ABCA6", "FBN1", "ZFPM2",
        # Note: ABCA10 also marks a homeostatic tendon fibroblast subtype
        # in Mimpen 2025 J Physiol (ABCA10hi); context distinguishes them
    ],

    "FAPs_Myofibroblast": [
        # Fibrogenic / contractile FAP state, encompasses scANVI TAGLNhi ADIRFhi
        # Note: TAGLN, ACTA2, MYL9 shared with vSMCs AND tendon SMMCs_Tendon
        # SERPINE1 and POSTN emphasis helps distinguish from vSMC context
        "POSTN", "FN1", "TGFBI", "TAGLN", "ADIRF",
        "ACTA2", "MYL9", "SERPINE1",
    ],

    "Fibroblasts_COL3A1hi": [
        # scANVI-validated: COL3A1hi LUMhi (structural matrix fibroblasts)
        "COL3A1", "LUM", "DCN", "FBLN1", "COL6A1",
    ],

    "Fibroblasts_PTGDShi": [
        # scANVI-validated: PTGDShi APODhi (quiescent interstitial FAPs)
        # Note: APOD+ subtype overlaps with tendon Tenocytes_APOD_FAPlike
        # (Kendal 2020 Tenocyte D) in tendon-containing samples
        "PTGDS", "APOD", "CFD", "CD55", "TPPP3",
    ],

    # --- Tenocyte / tendon-ligament subtypes (parent: Tendon_Ligament_Cells) ---
    #
    # Tenocyte subtype nomenclature follows Kendal et al. 2020 (Tenocyte A-E),
    # supplemented by Harvey 2019, Zhang 2023, and Mimpen 2024/2025.

    "Tenocytes_SCX": [
        # Core SCX+ midsubstance tenocytes (Kendal 2020, Harvey 2019)
        "SCX", "TNMD", "MKX", "FMOD", "THBS4", "KERA",
        "COL12A1", "TNC", "BGN",
        # Fibrillar lineage (Mimpen 2025 lifespan bioRxiv):
        # ABI3BP: collagen cross-linking regulator; GAS2: cytoskeletal regulator
        # SPARC: ECM organisation; highly enriched in KERA/TNMD tenocytes
        "ABI3BP", "SPARC", "GAS2",
    ],

    "Tenocytes_PTX3_Inflammatory": [
        # PTX3+ pro-inflammatory tenocytes (Kendal 2020 Tenocyte A)
        "PTX3", "KRT19", "CHI3L1", "CLDN11", "PENK",
        "SERPINE2", "CXCL1", "CXCL6",
    ],

    "Tenocytes_KRT7_ECM": [
        # KRT7+ ECM/microfibril tenocytes (Kendal 2020 Tenocyte B)
        "KRT7", "COL4A1", "POSTN", "THBS1",
        "TIMP3", "IGFBP5", "IGFBP7", "LTBP2",
        # FBLN1hi / repair markers (Mimpen 2025 FBLN1hi subtype, J Physiol):
        # FBLN1 (fibrillin-associated ECM), NOX4 (oxidative stress in repair),
        # CILP (cartilage intermediate layer protein, ECM remodelling)
        "FBLN1", "NOX4", "CILP",
    ],

    "Tenocytes_APOD_FAPlike": [
        # APOD+ FAP-like tenocytes (Kendal 2020 Tenocyte D)
        # Note: overlaps with Fibroblasts_PTGDShi in muscle datasets — these
        # likely represent the same quiescent mesenchymal state in different
        # tissue contexts
        "APOD", "LY6E", "CXCL14", "PLA2G2A",
        "FBLN2", "COL6A3",
        # PI16+ LCT overlap (Mimpen 2025 lifespan bioRxiv COL3A1 PI16 lineage)
        "PI16",
    ],

    "Tenocytes_TPPP3_Sheath": [
        # TPPP3/PRG4+ chondrogenic tenocytes / tendon sheath
        # (Kendal 2020 Tenocyte E, Harvey 2019 Sheath)
        "TPPP3", "PRG4", "CLU", "PRELP",
        "PCOLCE2", "CRTAC1", "CILP", "CILP2",
        # CREB5: lubricating fibroblast marker (Mimpen 2025 lifespan bioRxiv
        # CREB5 PRG4 fibroblasts at tendon sheath interface)
        "CREB5",
    ],

    "Tenocytes_NR4A1_Homeostatic": [
        # Healthy / homeostatic tendon fibroblasts (Mimpen 2025 J Physiol,
        # NR4A1hi subtype enriched in uninjured tissue; stress-response
        # and anti-apoptotic signalling via NR4A family nuclear receptors)
        "NR4A1", "NR4A3", "NAMPT", "SEMA4A",
    ],

    "Tenocytes_ADAM12_Repair": [
        # Repair / stress-response tenocytes (Mimpen 2025 J Physiol,
        # ADAM12hi subtype enriched in ruptured tendons; also relevant in
        # muscle injury and healing contexts)
        "ADAM12", "TNC", "POSTN", "THBS2",
    ],

    "Tendon_Stem_Progenitor_Cells": [
        # TSPCs — tendon-resident progenitors (Harvey 2019 Sheath cluster)
        "TPPP3", "PDGFRA", "PRG4", "CD44", "CD34",
        "S100A4", "AQP1",
    ],

    "Tendon_FAPs": [
        # Tendon-resident FAPs (Harvey 2019 T-FAPs; Mimpen 2025 lifespan bioRxiv
        # COL3A1 PI16 loose connective tissue fibroblast lineage)
        "PI16", "CD248", "COL14A1", "NID1",
        "DPT", "COL5A3", "CCL11",
    ],

    "Tendon_MTJ_Fibroblasts": [
        # Tendon-muscle junction fibroblasts (Mimpen 2025 lifespan bioRxiv
        # FGF14 THBS4 fibroblasts enriched at the myotendinous junction)
        "FGF14", "THBS4", "COL12A1", "FNDC1",
    ],

    "SMMCs_Tendon": [
        # Smooth muscle-mesenchymal cells in tendon (Kendal 2020 Tenocyte C)
        # Note: TAGLN, ACTA2, MYL9 shared with vSMCs and FAPs_Myofibroblast;
        # ITGA7, RGS5, VCAM1 help contextualise within tendon
        "ITGA7", "TAGLN", "MYL9", "ACTA2",
        "RGS5", "COL4A1", "VCAM1",
    ],

    "Enthesis_Fibrocartilage": [
        # Bone-tendon junction chondrocytes / fibrocartilage (Zhang 2023 eLife)
        "SOX9", "COL2A1", "ACAN", "COL9A1",
        "COL11A1", "CLEC3A", "WWP2", "MFGE8", "TNN",
    ],

    "Enthesoblasts": [
        # Transitional enthesis cells bridging tendon and fibrocartilage
        # (Zhang 2023 eLife)
        "COL12A1", "POSTN", "COMP", "TNMD", "THBS2",
    ],

    "Ligament_Fibroblasts": [
        # ACL / ligament fibroblasts (CellMarker)
        "COL1A1", "COL3A1", "DCN", "LUM",
        "ELN", "FBN1", "COMP", "PRG4",
    ],

    "Synovial_Fibroblasts_Lining": [
        # Synovial lining layer (CellMarker, HPA)
        "PRG4", "HAS1", "PDPN", "CD55",
    ],

    "Synovial_Fibroblasts_Sublining": [
        # Synovial sublining (CellMarker seq data)
        "THY1", "COL1A1", "PDGFRA", "ISLR", "GGT5",
    ],

    # --- Endothelial subtypes (parent: Endothelial_Cells) --------------

    "Endothelial_Arterial": [
        "GJA5", "HEY1", "DLL4", "EFNB2", "SOX17",
        "SEMA3G", "BMX", "DKK2",
        "CLU", "FBLN5",
        # Note: HEY1 also present in Satellite_Cells coarse (Notch signaling
        # shared); retained here as it is well-validated for arterial ECs
    ],

    "Endothelial_Capillary": [
        "CA4", "RGCC", "GPIHBP1", "BTNL9",
        "FCN3", "SLC9A3R2", "MFSD2A",
        # CD36: metabolic, enriched in muscle capillaries;
        # also in Macrophages_LAM — flagged as shared
        "CD36",
    ],

    "Endothelial_Venous": [
        "ACKR1", "SELE", "SELP", "EPHB4", "NR2F2",
        "PLVAP",
        "CD74",  # MHC-II, scANVI: ACKR1hi CD74hi
    ],

    "Endothelial_Lymphatic": [
        "PROX1", "FLT4", "PDPN",
        "CCL21", "TFF3", "MMRN1",
        # LYVE1 moved here from Macrophages_Resident: it is expressed on
        # tissue-resident macrophages but more definitive for lymphatic ECs;
        # removing it from Macrophages_Resident improves scoring accuracy
        "LYVE1",
    ],

    # --- Mural subtypes (parent: Mural_Cells) -------------------------

    "Pericytes": [
        "RGS5", "ABCC9", "KCNJ8", "PDGFRB", "NOTCH3",
        "HIGD1B", "NDUFA4L2",
        # CSPG4 (NG2) and MCAM (CD146): classical pericyte markers
        "CSPG4", "MCAM",
    ],

    "vSMCs": [
        "MYH11", "ACTA2", "CNN1", "TAGLN",
        "SMTN", "PLN", "RERGL",
    ],

    # --- Myeloid subtypes (parent: Macrophages_Monocytes) --------------
    # Note: DCs and neutrophils are included under this coarse compartment.

    "Macrophages_Resident": [
        "CD163", "MRC1", "C1QA", "C1QB", "C1QC",
        "F13A1", "FOLR2", "MAF", "STAB1",
        "SIGLEC1",
        # Note: LYVE1 removed (moved to Endothelial_Lymphatic where it is more
        # definitive; improves fine-level scoring accuracy)
    ],

    "Macrophages_LAM": [
        # Lipid-associated macrophages (relevant in aging/obese muscle)
        "TREM2", "LIPA", "GPNMB", "SPP1",
        "FABP5", "LGALS3", "APOE", "CTSD",
        # CD36: also in Endothelial_Capillary — flagged as shared
        "CD36",
    ],

    "Macrophages_Inflammatory": [
        "CD68", "S100A8", "S100A9", "S100A12",
        "IL1B", "TNF", "CCL3", "CCL4",
    ],

    "Monocytes_Classical": [
        "CD14", "S100A8", "S100A9", "LYZ", "VCAN", "FCN1",
    ],

    "Monocytes_Nonclassical": [
        "FCGR3A", "CDKN1C", "LST1", "MS4A7",
    ],

    "cDC1": [
        "XCR1", "CLEC9A", "BATF3", "IRF8", "CADM1",
    ],

    "cDC2": [
        "CD1C", "CLEC10A", "FCER1A", "HLA-DQA1", "ITGAX",
    ],

    "pDC": [
        "LILRA4", "IL3RA", "IRF7", "TCF4", "CLEC4C",
    ],

    "Neutrophils": [
        "FCGR3B", "CSF3R", "CXCR2", "S100A8", "S100A9",
    ],

    # --- Lymphoid T-cell subtypes (parent: T_Cells) --------------------

    "T_Cells_CD4": [
        "CD3D", "CD3E", "TRAC", "CD4", "IL7R",
        "TCF7", "LEF1",
    ],

    "T_Cells_CD8": [
        "CD3D", "CD3E", "CD8A", "CD8B",
        "GZMK", "GZMA",
    ],

    "T_Cells_Treg": [
        "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
    ],

    "T_Cells_TRM": [
        # Tissue-resident memory T cells
        "CD69", "ITGAE", "ZNF683", "CXCR6",
    ],

    # --- NK (parent: NK_Cells) -----------------------------------------

    "NK_Cells": [
        "NKG7", "GNLY", "KLRD1", "NCAM1",
        "PRF1", "GZMB", "KLRF1", "KLRB1",
        "XCL1", "XCL2", "FCER1G",
        # Note: NCAM1 (CD56) is the defining NK marker; removed from
        # Regenerative_Myonuclei fine to prevent NK ↔ myonuclei confusion
    ],

    # --- B lineage (parent: B_Cells) -----------------------------------

    "B_Cells": [
        "CD79A", "CD79B", "MS4A1", "CD19",
        "TCL1A", "PAX5", "BANK1", "BLK",
    ],

    "Plasma_Cells": [
        "MZB1", "JCHAIN", "XBP1", "IRF4",
        "PRDM1", "SDC1", "IGHG1",
    ],

    # --- Mast (parent: Mast_Cells) -------------------------------------

    "Mast_Cells": [
        "KIT", "CPA3", "TPSAB1", "TPSB2", "MS4A2",
        "HDC", "GATA2", "HPGDS",
    ],

    # --- Neural subtypes (parent: Neural) ------------------------------

    "Schwann_Cells_Myelinating": [
        "PLP1", "MPZ", "MBP", "PRX", "MAG",
    ],

    "Schwann_Cells_NonMyelinating": [
        # NGFR moved here from FAPs_Fibroblasts (more specific to Schwann cells)
        "CDH19", "SOX10", "S100B", "NGFR", "NFASC",
    ],

    "Neurons": [
        "SNAP25", "NRXN1", "SYT1", "NEFL", "RBFOX3",
    ],

    # --- Other (1:1 with coarse) ---------------------------------------

    "Adipocytes": [
        "PLIN1", "ADIPOQ", "LEP", "FABP4",
        "GPAM", "PDE3B", "CIDEC", "PCK1",
        "LPL", "PPARG",
    ],

    "Mesothelial_Cells": [
        "MSLN", "WT1", "UPK3B", "KRT19",
        "CALB2", "ITLN1", "LRRN4",
    ],

    "Erythrocytes": [
        "HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1",
    ],

    "Platelets": [
        "PF4", "PPBP", "ITGA2B", "GP9", "TUBB1", "ITGB3",
    ],
}


# ===========================================================================
# FINE -> COARSE mapping
# ===========================================================================

FINE_TO_COARSE: Dict[str, str] = {
    # Myofibers
    "Myofibers_type_I": "Myofibers",
    "Myofibers_type_IIA": "Myofibers",
    "Myofibers_type_IIX": "Myofibers",
    # Specialised myonuclei
    "NMJ_Myonuclei": "Myonuclei_Specialised",
    "MTJ_Myonuclei": "Myonuclei_Specialised",
    "Regenerative_Myonuclei": "Myonuclei_Specialised",
    "Immature_Myonuclei": "Myonuclei_Specialised",
    # Satellite cells
    "Satellite_Cells_Quiescent": "Satellite_Cells",
    "Satellite_Cells_Activated": "Satellite_Cells",
    # FAPs / Fibroblasts
    "FAPs": "FAPs_Fibroblasts",
    "FAPs_Adipogenic": "FAPs_Fibroblasts",
    "FAPs_Myofibroblast": "FAPs_Fibroblasts",
    "Fibroblasts_COL3A1hi": "FAPs_Fibroblasts",
    "Fibroblasts_PTGDShi": "FAPs_Fibroblasts",
    # Tendon / ligament (separate coarse compartment — NOT FAPs_Fibroblasts)
    # Removed: "Tenocytes": "FAPs_Fibroblasts" — tenocytes have a distinct
    # SCX+ developmental lineage and should not be grouped with FAPs
    "Tenocytes_SCX": "Tendon_Ligament_Cells",
    "Tenocytes_PTX3_Inflammatory": "Tendon_Ligament_Cells",
    "Tenocytes_KRT7_ECM": "Tendon_Ligament_Cells",
    "Tenocytes_APOD_FAPlike": "Tendon_Ligament_Cells",
    "Tenocytes_TPPP3_Sheath": "Tendon_Ligament_Cells",
    "Tenocytes_NR4A1_Homeostatic": "Tendon_Ligament_Cells",
    "Tenocytes_ADAM12_Repair": "Tendon_Ligament_Cells",
    "Tendon_Stem_Progenitor_Cells": "Tendon_Ligament_Cells",
    "Tendon_FAPs": "Tendon_Ligament_Cells",
    "Tendon_MTJ_Fibroblasts": "Tendon_Ligament_Cells",
    "SMMCs_Tendon": "Tendon_Ligament_Cells",
    "Enthesoblasts": "Tendon_Ligament_Cells",
    "Ligament_Fibroblasts": "Tendon_Ligament_Cells",
    "Synovial_Fibroblasts_Lining": "Tendon_Ligament_Cells",
    "Synovial_Fibroblasts_Sublining": "Tendon_Ligament_Cells",
    # Chondrocytes (new coarse compartment)
    "Enthesis_Fibrocartilage": "Chondrocytes",
    # Endothelial
    "Endothelial_Arterial": "Endothelial_Cells",
    "Endothelial_Capillary": "Endothelial_Cells",
    "Endothelial_Venous": "Endothelial_Cells",
    "Endothelial_Lymphatic": "Endothelial_Cells",
    # Mural
    "Pericytes": "Mural_Cells",
    "vSMCs": "Mural_Cells",
    # Myeloid (all map to Macrophages_Monocytes coarse)
    "Macrophages_Resident": "Macrophages_Monocytes",
    "Macrophages_LAM": "Macrophages_Monocytes",
    "Macrophages_Inflammatory": "Macrophages_Monocytes",
    "Monocytes_Classical": "Macrophages_Monocytes",
    "Monocytes_Nonclassical": "Macrophages_Monocytes",
    "cDC1": "Macrophages_Monocytes",
    "cDC2": "Macrophages_Monocytes",
    "pDC": "Macrophages_Monocytes",
    "Neutrophils": "Macrophages_Monocytes",
    # T cells
    "T_Cells_CD4": "T_Cells",
    "T_Cells_CD8": "T_Cells",
    "T_Cells_Treg": "T_Cells",
    "T_Cells_TRM": "T_Cells",
    # NK
    "NK_Cells": "NK_Cells",
    # B lineage
    "B_Cells": "B_Cells",
    "Plasma_Cells": "B_Cells",
    # Mast
    "Mast_Cells": "Mast_Cells",
    # Neural
    "Schwann_Cells_Myelinating": "Neural",
    "Schwann_Cells_NonMyelinating": "Neural",
    "Neurons": "Neural",
    # Other (1:1)
    "Adipocytes": "Adipocytes",
    "Mesothelial_Cells": "Mesothelial_Cells",
    "Erythrocytes": "Erythrocytes",
    "Platelets": "Platelets",
}


# ===========================================================================
# NEGATIVE MARKERS — genes that should NOT be expressed in a given type
#
# Used by hierarchical_score_and_assign() to penalise misassignment.
# Only defined for the coarse level; the most important disambiguation pairs.
# ===========================================================================

NEGATIVE_MARKERS: Dict[str, List[str]] = {
    # ---- Non-immune cell types: should not express PTPRC (CD45) -------
    "Myofibers": [
        "PTPRC", "PDGFRA", "PECAM1", "CDH5", "PAX7",
    ],
    "Myonuclei_Specialised": [
        "PTPRC", "PDGFRA", "PECAM1", "CDH5",
    ],
    "Satellite_Cells": [
        # Satellite cells are PAX7+ but should NOT express mature myofiber
        # genes or mesenchymal lineage markers at high levels
        "MYH7", "MYH1", "MYH2", "TTN", "NEB",
        "PDGFRA", "PECAM1", "PTPRC",
    ],
    "FAPs_Fibroblasts": [
        "PECAM1", "CDH5", "PAX7", "PTPRC", "TTN", "NEB",
    ],
    "Endothelial_Cells": [
        "PDGFRA", "COL1A1", "PAX7", "PTPRC", "TTN",
    ],
    "Mural_Cells": [
        # vSMCs share ACTA2/TAGLN with myofibroblasts — disambiguate via
        # negative markers rather than positive ones
        "PECAM1", "CDH5", "TTN", "NEB", "PTPRC", "PAX7",
    ],
    "Neural": [
        "PECAM1", "CDH5", "PDGFRA", "PTPRC", "TTN",
    ],
    "Adipocytes": [
        "PECAM1", "CDH5", "PDGFRA", "PTPRC", "PAX7", "TTN",
    ],
    "Mesothelial_Cells": [
        "PECAM1", "PTPRC", "TTN", "PAX7",
    ],
    "Tendon_Ligament_Cells": [
        "PTPRC", "PECAM1", "CDH5", "PAX7",
        # Should not express mature muscle contractile genes
        "TTN", "NEB", "MYH7", "MYH1", "MYH2",
    ],
    "Chondrocytes": [
        "PTPRC", "PECAM1", "CDH5", "PAX7",
        "TTN", "NEB",
        # PDGFRA: distinguishes FAPs from chondrocytes in enthesis context
        "PDGFRA",
    ],
    # ---- Immune cell types: should express PTPRC but not stromal ------
    "Macrophages_Monocytes": [
        "CD3D", "CD3E", "CD79A", "MS4A1", "PAX7", "PDGFRA", "PECAM1",
    ],
    "T_Cells": [
        "CD79A", "MS4A1", "CD68", "CSF1R", "PECAM1", "PDGFRA", "KIT",
    ],
    "NK_Cells": [
        "CD3D", "CD3E", "CD79A", "MS4A1", "CD68", "PECAM1",
    ],
    "B_Cells": [
        "CD3D", "CD3E", "CD68", "CSF1R", "PECAM1", "PAX7",
    ],
    "Mast_Cells": [
        "CD3D", "CD3E", "CD79A", "CD68", "PECAM1",
    ],
}

NUISANCE_MARKERS: Dict[str, List[str]] = {
    "myofibre": ["TTN", "NEB", "DES", "DMD", "ACTA1", "RYR1",
        "MYH7", "MYH1", "MYH2", "TNNT1", "TNNT3", "TNNC1", "TNNC2", "TNNI1", "TNNI2", "MYBPC1",  "MYL1", "MYL2",
        "CKM","MB"],
    "histones": [
        "HIST1H1C","HIST1H1D","HIST1H2AC","HIST1H2AD","HIST1H2AE","HIST1H2AG",
        "HIST1H2AJ","HIST1H2AK","HIST1H2AL","HIST1H2AM","HIST1H2BC","HIST1H2BD",
        "HIST1H2BE","HIST1H2BF","HIST1H2BG","HIST1H2BH","HIST1H2BI","HIST1H2BJ",
        "HIST1H3A","HIST1H3B","HIST1H3C","HIST1H3D","HIST1H3E","HIST1H3F",
        "HIST2H2AA3","HIST2H2AC","HIST2H2BE","HIST2H3A","HIST2H3C","HIST2H4A"
    ],
    "cell_cycle": [
        "MKI67","TOP2A","UBE2C","CDC20","CCNB1","CCNB2","CCNA2","CDK1",
        "CENPF","CENPA","AURKB","BIRC5","PCNA","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","TYMS"
    ],
    "ambient_epithelium_keratins": [
        "KRT1","KRT5","KRT10","KRT14","KRT15","KRT17","KRT18","KRT19"
    ],
    "ambient_immunoglobulin": ["IGKC","IGHG1","IGHG2","IGHG3","IGHG4","IGLC1","IGLC2","IGLC3","IGLC7"],
    "ambient_plasma_proteins": ["ALB","APOA1","APOA2","APOE","AHSG","FGA","FGB","FGG"],
    "erythrocytes": [
        "HBA1", "HBA2", "HBB", "ALAS2", "SLC4A1",
    ],

    "platelets": [
        "PF4", "PPBP", "ITGA2B", "GP9", "TUBB1", "ITGB3",
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
    filtered = {}
    for cell_type, genes in marker_dict.items():
        valid = [g for g in genes if g in var_set]
        if len(valid) >= 2:
            filtered[cell_type] = valid
    return filtered


def get_dotplot_markers() -> Dict[str, List[str]]:
    """Return a curated subset of markers suitable for a summary dotplot.

    These are 2-4 highly specific genes per major cell type, chosen for
    their discriminative power in skeletal muscle and tendon.

    Returns
    -------
    dict
        Cell type -> short gene list for dotplot display.
    """
    return {
        # Myofibers
        "Type I Myocytes": ["MYH7", "TNNT1", "MYL2"],
        "Type IIA Myocytes": ["MYH2", "ANKRD2"],
        "Type IIX Myocytes": ["MYH1", "MYLK2", "MYL1"],
        "NMJ Myonuclei": ["CHRNE", "COLQ", "COL13A1"],
        "MTJ Myonuclei": ["COL22A1", "PIEZO2", "LAMA2"],
        # Stem / mesenchymal
        "Satellite Cells": ["PAX7", "MYF5", "CALCR", "EGFR"],
        "FAPs": ["PDGFRA", "THY1", "DCN"],
        "Tenocytes": ["TNMD", "SCX", "KERA"],
        "Tendon Sheath": ["PRG4", "TPPP3", "CREB5"],
        # Vascular
        "EC Arterial": ["GJA5", "DLL4"],
        "EC Capillary": ["CA4", "RGCC"],
        "EC Venous": ["ACKR1", "SELP"],
        "EC Lymphatic": ["PROX1", "LYVE1"],
        "Pericytes": ["RGS5", "KCNJ8", "CSPG4"],
        "vSMCs": ["MYH11", "CNN1"],
        # Immune — myeloid
        "Resident Mac": ["CD163", "F13A1", "FOLR2"],
        "LAM": ["TREM2", "GPNMB", "SPP1"],
        "Inflammatory Mac": ["S100A8", "IL1B"],
        "Monocytes": ["CD14", "LYZ", "FCN1"],
        "cDC2": ["CD1C", "CLEC10A"],
        # Immune — lymphoid
        "T Cells": ["CD3D", "CD3E"],
        "NK Cells": ["NKG7", "GNLY"],
        "B Cells": ["CD79A", "MS4A1"],
        "Plasma Cells": ["MZB1", "JCHAIN"],
        "Mast Cells": ["KIT", "TPSB2"],
        # Other
        "Schwann Cells": ["PLP1", "MPZ"],
        "Adipocytes": ["PLIN1", "ADIPOQ"],
        "Chondrocytes": ["SOX9", "COL2A1", "ACAN"],
        "Mesothelial": ["MSLN", "WT1"],
        "Erythrocytes": ["HBA1", "HBB"],
        "Platelets": ["PF4", "PPBP"],
    }


def get_diagnostic_markers() -> Dict[str, List[str]]:
    """Return the 1-3 most pathognomonic genes per cell type.

    These are suitable for quick validation on feature plots —
    if a cluster expresses these genes, the identity is near-certain.

    Returns
    -------
    dict
        Cell type -> 1-3 most diagnostic genes.
    """
    return {
        "Myofibers": ["TTN", "DES"],
        "Type_I": ["MYH7", "MYL2"],
        "Type_IIA": ["MYH2"],
        "Type_IIX": ["MYH1", "MYL1"],
        "NMJ": ["CHRNE", "COL13A1"],
        "MTJ": ["COL22A1"],
        "Satellite_Cells": ["PAX7", "EGFR"],
        "FAPs": ["PDGFRA"],
        "Tenocytes": ["TNMD", "SCX"],
        "Chondrocytes": ["SOX9", "COL2A1"],
        "Endothelial": ["PECAM1"],
        "Pericytes": ["RGS5", "KCNJ8"],
        "vSMCs": ["MYH11"],
        "Macrophages": ["CD68", "AIF1"],
        "T_Cells": ["CD3D"],
        "NK_Cells": ["NKG7", "NCAM1"],
        "B_Cells": ["CD79A"],
        "Mast_Cells": ["KIT"],
        "Schwann": ["PLP1"],
        "Adipocytes": ["PLIN1"],
        "Erythrocytes": ["HBA1"],
    }
