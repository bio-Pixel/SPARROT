
# 🧠 SPARROT Case Study  

# Spatial Transcriptomics Showcase: Myocardial Infarction Sample (P9)

For this showcase, we utilize spatial transcriptomics data from [Kuppe et al., 2022](https://www.nature.com/articles/s41586-022-05060-x), who generated a comprehensive spatial multi-omic atlas of human myocardial infarction. Specifically, we analyze a 10X Visium spatial slide (P9) representing an ischemic zone of the human heart post-infarction. The Seurat-formatted spatial object can be downloaded directly from the following link:

**🔗 [10X_Visium_ACH0012.tar.gz (Zenodo)](https://zenodo.org/records/6580069/files/10X_Visium_ACH0012.tar.gz?download=1)**

The corresponding cell type composition, estimated using *cell2location*, is available at:

**🔗 [cell2location annotations (.rds)](https://drive.google.com/file/d/1YWocGsNZ929NKrZP0Jbfi-iBG9c2JR4j/view?usp=drive_link)**

---

## 📦 Workflow Summary

### 1. Load data and create `SparrotObj`

```r
library(SPARROT)
library(Seurat)

# Load Seurat spatial transcriptomics object and cell-type deconvolution matrix 
seu <- readRDS("ACH0012.rds")
cpm <- readRDS("P9_CellProb_cell2location.rds")

# Create SPARROT object from Seurat object
cc <- convertSeuratToSparrot(seu, cell_prob = cpm)
```
You can also create the SPARROT object by:
```r
cc = createSparrotObj(coords = GetTissueCoordinates(seu),
                      cell_prob = cpm,
                      expr = GetAssayData(seu, layer = 'data'))
cc
#> An object of class 'SparrotObj'
#> 
#> Number of spots/cells:  4361 
#> Number of genes:        15972 
#> Number of cell types:   11 
#> Meta data columns:
#>   • coords:     row, col 
#>   • binarized:  bin_Adipocyte, bin_Cardiomyocyte, bin_Endothelial, bin_Fibroblast, bin_Lymphoid, bin_Mast, bin_Myeloid, bin_Neuronal, bin_Pericyte, bin_Cycling.cells, bin_vSMCs 
#> 
#> Use @meta.data, @expr, @cell_prob, or accessor methods to explore.
```

---

### 2. Visualize Multi-celltype Probabilities

```r
plotMultiCellTypeProb(cc, celltype = c("Endothelial", "Fibroblast", "Cardiomyocyte"),
                      color = c(Cardiomyocyte = "#2D81FF", Fibroblast = "#00B37F", Endothelial = "#FF6A00"),
                      outline = FALSE, coord.fixed = TRUE)
```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/P9_cardio_prop.png?raw=true" width="500"/>

---

### 3. Evaluate Pairwise Celltype Overlap

```r
evaluate_overlap_metrics(bin1 = as.logical(cc@meta.data[, "bin_Cardiomyocyte"]),
                         bin2 = as.logical(cc@meta.data[, "bin_Endothelial"]),
                         coords = cc@coords)
```

```r
pm <- computePairwiseCelltypeOverlap(cc, metric = "dice", ncore = 10)
net <- subset(pm, pvalue < 0.05)
saveRDS(net, "pairwise_dice_sig.rds")
```

---

### 4. Build Colocalization Network

```r
library(igraph)
library(tidygraph)
library(ggraph)

g <- graph_from_data_frame(net[, c("celltype1", "celltype2", "dice")], directed = FALSE)
E(g)$weight <- net$dice
tg <- as_tbl_graph(g) %>% mutate(community = as.factor(group_louvain()))
```

```r
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray60", alpha = 0.7) +
  geom_node_point(aes(color = community), size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void()
```

📷 **Output:**

![](vignettes/figs/ACH0012_cellcolocal_map.png)

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
