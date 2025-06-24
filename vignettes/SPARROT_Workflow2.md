
# 🟢 Spatial Transcriptomics Showcase: Human Lymph Node

10x Genomics obtained fresh frozen human lymph node tissue from BioIVT Asterand Human Tissue Specimens. The tissue was embedded and cryosectioned as described in Visium Spatial Protocols - Tissue Preparation Guide (Demonstrated Protocol CG000240). The Seurat-formatted spatial object can be downloaded directly from the following link:

**🔗 [lymph_node_seuobj.rds (.rds)](https://drive.google.com/file/d/1Lwx0_M8dgNUNJ646UFMdbe9kT3wgp3CN/view?usp=drive_link)**

The corresponding cell type composition, estimated using *cell2location*, is available at:

**🔗 [cell2location annotations (.rds)](https://drive.google.com/file/d/1F45TjbJHtta5cPPth8YIZDIw5M7uNyk7/view?usp=drive_link)**

---

## 📦 Workflow Summary

### 1. Load data and create `SparrotObj`

```r
library(SPARROT)
library(Seurat)
library(ggplot2)
library(igraph)
library(tidygraph)
library(ggraph)


# Load Seurat spatial transcriptomics object 
seu <- readRDS("lymph_node_seuobj.rds")
```
```r
# Load cell-type deconvolution matrix 
cpm <- readRDS("LM_CellProb_cell2location.rds")
view(cpm)
```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/LM_cpm.png?raw=true" width="1000"/>


```r
identical(colnames(seu), rownames(cpm))
```
```r
#> [1] TRUE
```

Create SPARROT object from Seurat object:
```r
cc <- convertSeuratToSparrot(seu, cell_prob = cpm)
```
You can also create the SPARROT object by:
```r
cc = createSparrotObj(coords = GetTissueCoordinates(seu),
                      cell_prob = cpm,
                      expr = GetAssayData(seu, layer = 'data'))
cc
```
```r
#> An object of class 'SparrotObj'

#> Number of spots/cells:  4039 
#> Number of genes:        33538 
#> Number of cell types:   34 
#> Meta data columns:
#>   • coords:     row, col 
#>   • binarized:  bin_B_Cycling, bin_B_GC_DZ, bin_B_GC_LZ, bin_B_GC_prePB, bin_B_IFN, bin_B_activated, bin_B_mem, bin_B_naive, bin_B_plasma, bin_B_preGC, bin_DC_CCR7., bin_DC_cDC1, bin_DC_cDC2, bin_DC_pDC, bin_Endo, bin_FDC, bin_ILC, bin_Macrophages_M1, bin_Macrophages_M2, bin_Mast, bin_Monocytes, bin_NK, bin_NKT, bin_T_CD4., bin_T_CD4._TfH, bin_T_CD4._TfH_GC, bin_T_CD4._naive, bin_T_CD8._CD161., bin_T_CD8._cytotoxic, bin_T_CD8._naive, bin_T_TIM3., bin_T_TfR, bin_T_Treg, bin_VSMC 
#>   • others:     orig.ident, nCount_Spatial, nFeature_Spatial 
#> 
#> Use @meta.data, @expr, @cell_prob, or accessor methods to explore.
```

---
### 2. Visualize Multi-celltype Probabilities

```r
plotMultiCellTypeProb(cc, celltype =c('T_CD4._naive', 'B_naive', 'FDC'),
                      color = c(`T_CD4._naive` = "#2D81FF", B_naive = "#00B37F", FDC = '#FF6A00'),
                      outline = F, coord.fixed = T) 
```

<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/LM_Miltiprop.png?raw=true" width="500"/>

---

### 3. Pairwise Spatial Overlap Significance Analysis for Specified Cell Types
The *computePairwiseCelltypeOverlap()* function is designed to efficiently calculate spatial overlap (**metric = c("dice", "jaccard", "mcc")**) for **all cell type pairs in a given SparrotObj**. When analyzing datasets with many cell types, the number of pairwise comparisons grows rapidly, making parallel computing highly beneficial. If **ncore > 1** is specified, parallel computation is enabled. 
```r
# Run overlap computation using 10 cores
pm <- computePairwiseCelltypeOverlap(cc, metric = "dice", ncore = 10)
```
If ncore = NULL or ncore = 1, the function defaults to sequential execution.
This function returns a **data frame** with each cell type pair, their **overlap score**, and **corresponding p-value**. 
For **downstream analysis**, we selected statistically significant edges:
```r
# Filter pairs with significant spatial overlap
net_filtered <- subset(pm, pvalue < 0.05)
net_filtered
```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/LM_net_data.png?raw=true" width="500"/>

## 🌐 Visualizing Cell-Cell Spatial Co-localization Networks
To visualize how different cell types co-localize in the tissue space, we built a cell-type network where:
Nodes represent cell types
Edges represent significant spatial co-localization
Edge weight corresponds to the Dice score
Communities are detected using Louvain clustering to reveal spatially interacting groups
```r
library(igraph)
library(tidygraph)
library(ggraph)

# Create igraph object from filtered results
g <- graph_from_data_frame(net_filtered[, c("celltype1", "celltype2", "dice")], directed = FALSE)

# Assign edge weight
E(g)$weight <- net_filtered$dice

# Convert to tidygraph for community detection
tg <- as_tbl_graph(g) %>%
  mutate(community = as.factor(group_louvain()))

# Plot using ggraph
set.seed(123)
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray60", alpha = 0.7) +
  geom_node_point(aes(color = community), size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void() +
  theme(legend.position = "none", aspect.ratio = 0.8)
```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/LM_net.png?raw=true" width="600"/>

---

## 📋 Session Info


```r
sessionInfo()
```

```
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26100)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8
[2] LC_CTYPE=Chinese (Simplified)_China.utf8
[3] LC_MONETARY=Chinese (Simplified)_China.utf8
[4] LC_NUMERIC=C
[5] LC_TIME=Chinese (Simplified)_China.utf8

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
[1] ggraph_2.2.1       tidygraph_1.3.1    igraph_2.0.3
[4] ggplot2_3.5.0      Seurat_5.0.3       SeuratObject_5.0.1
[7] sp_2.1-3           SPARROT_0.5.0
```


---
**SPARROT v0.5.0** — Spatial Proximity Analysis of Regional Relationships, Overlap, and Transcriptome  
Created by Xu Pan | [GitHub](https://github.com/bio-Pixel/SPARROT)
