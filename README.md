# Single-cell transcriptomic analysis of granulosa cells in a mouse PCOS model

Analysis code accompanying our manuscript on single-cell dissection of granulosa
cell (GC) heterogeneity in polycystic ovary syndrome (PCOS), identifying
***Rora*** as a key regulator of GC state transitions.

## Overview

This repository contains the analysis pipeline for our single-cell RNA-seq
dataset (mouse ovary / PCOS model), together with external validation on the
public datasets **GSE240688** and **GSE98595**. Scripts are organized by figure
(`Figures/FigureN/`); each folder holds the code used to generate the
corresponding main and supplementary panels.

### Analysis workflow

1. **Data preprocessing and cell-type annotation** — quality control,
   doublet removal, Harmony batch correction, clustering and annotation of
   major ovarian cell types; spatial transcriptomic analysis.
2. **Granulosa cell subclustering** — GC subsetting, removal of
   hemoglobin-high (erythrocyte-contaminated) clusters, re-annotation,
   subgroup differential expression and GO enrichment, and Augur-based
   cell-state prioritization.
3. **Trajectory and pseudotime analysis** — Monocle3, PAGA/DPT, CytoTRACE2,
   Palantir and CellRank2, including diffusion maps, root/terminal-state
   selection, branch detection and fate-probability projection.
4. **Gene regulatory network analysis** — hdWGCNA co-expression modules and
   pySCENIC transcription-factor regulons, integrated by intersection
   analysis (Venn) together with the pseudotime and prioritization results
   to nominate candidate regulators.
5. **In silico perturbation** — CellOracle-based virtual knockout of
   candidate genes.
6. **Cell–cell communication** — CellChat and NicheNet ligand–receptor
   inference with sender/receiver ranking.
7. **Drug sensitivity screening** — drug2cell drug-sensitivity scoring on
   cell populations of interest.
8. **External validation** — differential expression and ssGSEA of candidate
   signatures on GSE240688 and GSE98595.

## Repository structure

```text
Figures/
├── Figure2/  # QC, cell-type annotation, composition and spatial analysis
│   ├── FigureS1_UMAP_after_doublet_removal.ipynb
│   ├── cell_annotation_proportion.ipynb
│   ├── FigureS2_celltype_group_DEA.ipynb
│   ├── FigureS2_celltype_DEG_GSEA_clean.R
│   └── spatial_analysis.R
├── Figure3/  # Granulosa cell subclustering and subgroup characterization
│   ├── GC_subclustering_RBCremoval.ipynb
│   ├── GC_annotation_DEA.ipynb
│   ├── GC_GO_enrichment.R
│   └── GC_augur.R
├── Figure4/  # Trajectory, pseudotime and co-expression analysis
│   ├── monocle3.R
│   ├── palantir.ipynb
│   ├── dpt_cytotrace2_palantir.ipynb
│   ├── FigureS3_cellrank2_visualization.ipynb
│   ├── FigureS4_pseudotime_palantir.R
│   └── Figure4_S5-7_hdWGCNA.R
├── Figure5/  # Regulon analysis, candidate nomination and external validation
│   ├── pySCENIC_TF_analysis.R
│   ├── pySCENIC_post_analysis.R
│   ├── Rora_venn.R
│   ├── GSE240688_external_validation.R
│   └── GSE98595_ssGSEA_validation.R
├── Figure6/  # In silico perturbation
│   ├── celloracle_virtual_knockout.ipynb
│   └── virtual_knockout_clean.R
├── Figure7/  # Cell–cell communication
│   ├── cellchat.R
│   └── NicheNet.R
└── Figure8/  # Drug sensitivity screening
    ├── drug2cell.R
    └── drug2cell.ipynb
```

## Requirements

- **R**: Seurat, monocle3, WGCNA (hdWGCNA), AUCell/SCENIC downstream
  tooling, CellChat, NicheNet, Augur, clusterProfiler, ggplot2.
- **Python**: scanpy, anndata, palantir, cytotrace2-py, cellrank, scvelo,
  celloracle, drug2cell, pandas, matplotlib, seaborn.

Conda environments used during the study are noted as comments in the
individual scripts. All file paths are relative (`./`) and can be adapted
to the local directory layout.

## Usage

Run the scripts in the order implied by the figure numbering (Figure2 ->
Figure8). Intermediate `.h5ad`/`.csv` objects are read and written through
relative paths; the pipeline assumes the output folder of each stage is
present (or created by `os.mkdir` as in the scripts).

## Data availability

- Raw and processed scRNA-seq data of this study: [GSE268919](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268919).
- External validation datasets: [GSE240688](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE240688)
  and [GSE98595](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98595).

## Citation

If you use this code, please cite our manuscript: [DOI / link to be added upon
publication].

## License

[To be chosen before release, e.g. MIT for code.]
