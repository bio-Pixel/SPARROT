
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
- Compute spatial overlap (Dice-Sørensen coefficient, Jaccard index, Matthews correlation coefficient) between genes and cell types
- Permutation tests for colocalization significance
- Quantify spatial proximity using nearest-neighbor distance metrics
- Visualize multi-celltype spatial maps with transparency encoding

---

## 📘 Vignette

You can find a step-by-step usage example in the
1. [Myocardial Infarction](vignettes/SPARROT_Workflow1.md)
2. [Myocardial Infarction](vignettes/SPARROT_Workflow2.md)
3. [Myocardial Infarction](vignettes/SPARROT_Workflow3.md)

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
