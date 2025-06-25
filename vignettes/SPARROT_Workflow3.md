# 🩺 Spatial Transcriptomics Showcase: Human Pancreatic Ductal Adenocarcinoma, Stereo-seq

The pathological tumor size was 21 mm, with evidence of neural, vesicle, and para-pancreatic infiltration.The patient is a 58-year-old male diagnosed with moderately differentiated pancreatic ductal adenocarcinoma (PDAC) located in the head of the pancreas. Genetic testing revealed a KRAS G12R mutation. To investigate the spatial architecture of the tumor microenvironment, we performed spatial transcriptomic profiling using the **Stereo-seq** on this patient's resected tumor tissue. 100 × 100 DNA Nanoballs (DNB) (*bin100*, 50-μm diameter) were combined as a pseudo spot, with transcripts of the same gene aggregated within each spot. The SPARROT-formatted object with cell type composition, estimated by *Tangram*, can be downloaded directly from the following link:

**🔗 [PUMCH_D2_SparrotObj.rds (.rds)](https://drive.google.com/file/d/18OnrCikbKFvYnGHvR0oBA5EKH94DrORf/view?usp=sharing)**

---

## 📦 Workflow Summary

### 1. Load data `SparrotObj`

```r
library(SPARROT)
library(Seurat)
library(ggplot2)

# Load SparrotObj 
cc <- readRDS("PUMCH_D2_SparrotObj.rds")
```
```r
#> An object of class 'SparrotObj'

#> Number of spots/cells:  14036 
#> Number of genes:        19409 
#> Number of cell types:   31 
#> Meta data columns:
#>   • coords:     row, col 
#>   • binarized:  bin_KRT19., bin_CD8..T, bin_Fibroblast, bin_ADM, bin_Acinar, bin_Macrophage, bin_EC, bin_CD4..Th, bin_Ductal, bin_SMC, bin_Endocrine, bin_Treg, bin_Bm, bin_Neutrophil, bin_PC, bin_Tn, bin_Bn, bin_Stellate, bin_Mast, bin_NK, #>#>  bin_Monocyte, bin_CD4..Tfh, bin_DC, bin_Schwann, bin_NKT, bin_MDSC, bin_ILC, bin_Tuft, bin_GC.B, bin_Proliferating, bin_Adipocyte 

#> Use @meta.data, @expr, @cell_prob, or accessor methods to explore.
```

---

### 2. Visualize Multi-celltype Probabilities

```r
plotMultiCellTypeProb(cc, celltype = c( "Bm","PC","KRT19."  ), 
                      color = c( Bm = "#2D81FF", PC = "#00B37F", KRT19. = "#A259FF"),
                      outline = TRUE, concavity = 2, coord.fixed = TRUE)
```

<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/PUMCH_D2_multiProb-01.png?raw=true" width="500"/>

---

### 4. Visualize Gene expression density 

```r
spFeatureDensityPlot(cc, features = c("PDGFRA","RYR2","PECAM1"), outline = F)
```

<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/P9_cardio_expr.png?raw=true" width="1000"/>

---

### 5. Evaluation of Spatial Co-localization of Cardiomyocytes with Fibroblasts/Endothelial cells
*evaluate_overlap_metrics()* evaluates the spatial overlap between two binary spatial patterns — for instance, between two cell types or between a gene expression region and a cell type — in spatial transcriptomics data. It computes both the **overlap scores (Dice-Sørensen coefficient, Jaccard index, Matthews correlation coefficient)** and **permutation-based p-values** to assess statistical significance.

```r
evaluate_overlap_metrics(bin1 = as.logical(cc@meta.data[, "bin_Cardiomyocyte"]),
                         bin2 = as.logical(cc@meta.data[, "bin_Fibroblast"]),
                         coords = cc@coords)
```
```r
#>       dice p_dice   jaccard p_jaccard        mcc p_mcc
#>1 0.3724632      1 0.2288509         1 -0.4687531     1
```

> **Interpretation**: Although *Cardiomyocytes* and *Fibroblasts* show moderate spatial overlap (*Dice* = 0.37),  
> the lack of statistical significance (*p* = 1) and a negative *MCC* suggest their distributions are likely  
> **mutually exclusive** in this tissue region and **not spatially co-localized beyond chance**.

```r
evaluate_overlap_metrics(bin1 = as.logical(cc@meta.data[, "bin_Cardiomyocyte"]),
                         bin2 = as.logical(cc@meta.data[, "bin_Endothelial"]),
                         coords = cc@coords)
```
```r
#>       dice p_dice   jaccard p_jaccard       mcc p_mcc
#>1 0.6410665      0 0.4717423         0 0.1813916     0
```

> **Interpretation**: *Cardiomyocytes* and *Endothelial cells* exhibit a **strong spatial overlap**  
> (*Dice* = 0.64, *Jaccard* = 0.47), with all permutation-based p-values < 0.001.  
> This indicates a **statistically significant co-localization**, suggesting these two cell types  
> may occupy shared niches or interact closely in the infarct zone.

---

### 6. Evaluation of Spatial Co-localization of Cardiomyocytes with Gene Expression

```r
computeGeneCelltypeOverlap(cc, gene = "RYR2", celltype = "Fibroblast")
```
```r
#>fitting ...
#>  |===================================================================================================================| 100%
#>       dice p_dice   jaccard p_jaccard        mcc p_mcc
#>1 0.3614404      1 0.2205842         1 -0.3097191     1
```

> **Interpretation**: The spatial expression of *RYR2* shows only modest overlap with *Fibroblasts*  
> (*Dice* = 0.36), and all p-values are non-significant (*p* = 1), suggesting that *RYR2* is likely expressed  
> in regions that are **spatially distinct** from *Fibroblast*-enriched areas.

```r
computeGeneCelltypeOverlap(cc, gene = "PDGFRA", celltype = "Fibroblast")
```
```r
#>fitting ...
#>  |===================================================================================================================| 100%
#>       dice p_dice   jaccard p_jaccard      mcc p_mcc
#>1 0.5524716      0 0.3816654         0 0.231789     0
```

> **Interpretation**: The expression of *PDGFRA* exhibits a **significant spatial colocalization**  
> with *Fibroblasts* (*Dice* = 0.55, *p* < 0.001), indicating that *PDGFRA* may mark or functionally associate  
> with *Fibroblast* populations in this tissue context.


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
