
---
title: "Cardiomyocyte Spatial Proximity Analysis with SPARROT"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Cardiomyocyte Spatial Proximity Analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 5
)
library(SPARROT)
library(Seurat)
library(ggplot2)
library(ggraph)
library(tidygraph)
library(igraph)
```

### 1. Data Preparation

```{r load-data}
setwd("C:/Users/Pixel/Desktop/Project/PDAC_ecotype/test_ST/cardio_P9")
seu <- readRDS("ACH0012.rds")
cpm <- read.csv("adata_obs.csv", row.names = 1)
cpm <- cpm[colnames(seu), 7:17]
saveRDS(cpm, "CellProb_cell2location.rds")

cc <- convertSeuratToSparrot(seu, cell_prob = cpm)
```

### 2. Multi-celltype Visualization

```{r multi-prob-plot, fig.cap="Multi-celltype probability map of Endothelial, Fibroblast, and Cardiomyocyte"}
plotMultiCellTypeProb(cc, 
                      celltype = c('Endothelial','Fibroblast','Cardiomyocyte'),
                      color = c(Cardiomyocyte = "#2D81FF", Fibroblast = "#00B37F", Endothelial = '#FF6A00'),
                      outline = FALSE, coord.fixed = TRUE)
```

### 3. Pairwise Overlap Metrics

```{r pairwise-metrics}
evaluate_overlap_metrics(
  bin1 = as.logical(cc@meta.data[,"bin_Cardiomyocyte"]),
  bin2 = as.logical(cc@meta.data[,"bin_Endothelial"]),
  coords = cc@coords
)
```

```{r run-pairwise-overlap, eval=FALSE}
pm <- computePairwiseCelltypeOverlap(cc, "dice", ncore = 10)
saveRDS(pm, "pairwise_dice_sig.rds")
```

### 4. Colocalization Network Visualization

```{r load-network}
net <- readRDS("pairwise_dice_sig.rds")
net_filtered <- subset(net, pvalue < 0.05)
g <- graph_from_data_frame(net_filtered[, c("celltype1", "celltype2", "dice")], directed = FALSE)
E(g)$weight <- net_filtered$dice
tg <- as_tbl_graph(g) %>% mutate(community = as.factor(group_louvain()))
```

```{r colocalization-network, fig.cap="Cell-type colocalization network based on Dice coefficient"}
ggraph(tg, layout = "fr") +
  geom_edge_link(aes(width = weight), color = "gray60", alpha = 0.7) +
  geom_node_point(aes(color = community), size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  scale_edge_width(range = c(0.5, 3)) +
  theme_void() +
  theme(legend.position = "none", aspect.ratio = 0.8)
```

### 5. Session Info

```{r session-info}
sessionInfo()
R version 4.3.3 (2024-02-29 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] ggraph_2.2.1       tidygraph_1.3.1    igraph_2.0.3       ggplot2_3.5.0      Seurat_5.0.3       SeuratObject_5.0.1
[7] sp_2.1-3           SPARROT_0.5.0     

```
