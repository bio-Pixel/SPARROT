
# SPARROT <img src="man/figures/SPARROT_logo.png" align="right" height="140"/>

**SPARROT** (*Spatial Proximity Analysis of Regional Relationships, Overlap, and Transcriptome*)  
A toolkit for analyzing and visualizing spatial colocalization and interactions in spatial transcriptomics data.

---

## 🐦 Installation

```r
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("bio-Pixel/SPARROT")
```

---

## 📦 Features

- Seamlessly convert Seurat spatial objects to `SparrotObj`
- Binarize cell type probability using GMM or main-valley threshold
- Compute spatial overlap (Dice, Jaccard, MCC) between genes and cell types
- Quantify spatial proximity using nearest-neighbor distance metrics
- Visualize multi-celltype spatial maps with transparency encoding
- Interface with [Nebulosa](https://github.com/constantAmateur/Nebulosa) for smooth feature density plots
- Permutation tests for colocalization significance
- Concave hull outlining for spatial domains

---

## 📘 Vignette

You can find a step-by-step usage example in the [SPARROT Workflow Vignette](vignettes/SPARROT_Workflow.Rmd).  
*(Or replace with HTML if using `pkgdown`)*

---

## 🧬 Example

```r
library(SPARROT)

# Convert Seurat to SPARROT
obj <- convertSeuratToSparrot(seurat_obj, cell_prob = pred_matrix)

# Plot spatial probabilities
plotMultiCellTypeProb(obj, celltype = c("B", "Plasma"))

# Compute gene-celltype overlap
computeGeneCelltypeOverlap(obj, gene = "MZB1", celltype = "Plasma")
```
