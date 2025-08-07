# -*- coding: utf-8 -*-
"""
Constants and configuration for EspressoPro package.
"""

from __future__ import annotations

# Reference datasets for best-localised tracks
_REF_SIMPLIFIED = ("BestSimplified",)   # use *only* best‑localised tracks
_REF_DETAILED   = ("BestDetailed",)

# Detailed cell type labels
_DETAILED_LABELS = [
    "CD14_Mono", "CD16_Mono", "cDC1", "cDC2", "ILC", "Macrophage", "MkP",
    "HSC_MPP", "MEP", "ErP", "GMP", "LMPP", "Pre-Pro-B", "EoBaMaP", "Plasma",
    "B_Naive", "B_Memory", "Immature_B", "CD4_T_Naive", "CD4_T_Memory", "CD4_CTL",
    "Treg", "CD8_T_Naive", "CD8_T_Memory", "MAIT", "NK_CD56_dim", "NK_CD56_bright",
    "Erythroblast", "Stroma", "pDC", "gdT", "dnT", "Myeloid_progenitor",
    "Pre-B", "Pro-B",
]

# Simplified cell type classes mapping
SIMPLIFIED_CLASSES = {
    lbl: [f"{ref}.{lbl}.predscore" for ref in _REF_SIMPLIFIED] for lbl in (
        "NK", "HSPC", "Erythroid", "pDC", "Monocyte", "Myeloid",
        "CD4_T", "CD8_T", "B", "cDC", "Other_T", "ILC", "Macrophage", "Stroma",
    )
}

# Detailed cell type classes mapping  
DETAILED_CLASSES = {
    lbl: [f"{ref}.{lbl}.predscore" for ref in _REF_DETAILED] for lbl in _DETAILED_LABELS
}

# Parent-child mappings for hierarchical annotation
SIMPLIFIED_PARENT_MAP = {
    "Immature":         ["CommonSimplified.HSPC.predscore"],
    "Mature":           [f"CommonSimplified.{l}.predscore" for l in SIMPLIFIED_CLASSES if l != "HSPC"],
    "Immature/Mature":  [f"CommonSimplified.{l}.predscore" for l in SIMPLIFIED_CLASSES],
}

# Mapping from Broad/Simplified "parent" to allowed Detailed classes
DETAILED_PARENT_MAP = {
    "HSPC": [
        "CommonDetailed.HSC_MPP.predscore",
        "CommonDetailed.ErP.predscore", "CommonDetailed.GMP.predscore",
        "CommonDetailed.LMPP.predscore", "CommonDetailed.EoBaMaP.predscore",
        "CommonDetailed.Pre_Pro_B.predscore", "CommonDetailed.MkP.predscore",
        "CommonDetailed.MEP.predscore", "CommonDetailed.Pro_B.predscore",
        "CommonDetailed.Myeloid_progenitor.predscore",
    ],
    "Erythroid": [
        "CommonDetailed.ErP.predscore", "CommonDetailed.Erythroblast.predscore",
    ],
    "pDC": [
        "CommonDetailed.EoBaMaP.predscore", "CommonDetailed.pDC.predscore",
    ],
    "Monocyte": [
        "CommonDetailed.CD14_Mono.predscore", "CommonDetailed.CD16_Mono.predscore",
        "CommonDetailed.cDC.predscore",       "CommonDetailed.Macrophage.predscore",
    ],
    "Myeloid": ["CommonDetailed.Myeloid_progenitor.predscore"],
    "cDC":     ["CommonDetailed.cDC1.predscore", "CommonDetailed.cDC2.predscore"],
    "ILC":     ["CommonDetailed.ILC.predscore"],
    "Other_T": ["CommonDetailed.dnT.predscore",  "CommonDetailed.gdT.predscore"],
    "NK":      ["CommonDetailed.NK_CD56_dim.predscore",
                "CommonDetailed.NK_CD56_bright.predscore"],
    "CD4_T": [
        "CommonDetailed.CD4_T_Naive.predscore",  "CommonDetailed.CD4_T_Memory.predscore",
        "CommonDetailed.CD4_CTL.predscore",      "CommonDetailed.Treg.predscore",
    ],
    "CD8_T": [
        "CommonDetailed.CD8_T_Naive.predscore",  "CommonDetailed.CD8_T_Memory.predscore",
        "CommonDetailed.MAIT.predscore",         "CommonDetailed.gdT.predscore",
    ],
    "B": [
        "CommonDetailed.Plasma.predscore",       "CommonDetailed.B_Naive.predscore",
        "CommonDetailed.B_Memory.predscore",     "CommonDetailed.Immature_B.predscore",
        "CommonDetailed.Pre_B.predscore",        "CommonDetailed.Pro_B.predscore",
    ],
    "Stroma": ["CommonDetailed.Stroma.predscore"],
}

# Mast cell marker signatures
MAST_POS = ['FcεRIα', 'CD117', 'CD62L']          # ↑ mast
MAST_NEG = [                                     # ↓ mast
    'CD303', 'CD304', 'CD123', 'CD34', 'CD8', 'CD4',
    'CD138', 'CD7', 'CD10', 'CD11b', 'CD5', 'CD141', 'CD1c',
]
