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

### 2. Visualize Probabilities of Memory B Cells (Bm), Plasma Cells (PC), and KRT19+ Cells (Tumor)

```r
plotMultiCellTypeProb(cc, celltype = c( "Bm","PC","KRT19."  ), 
                      color = c( Bm = "#2D81FF", PC = "#00B37F", KRT19. = "#A259FF"),
                      outline = TRUE, concavity = 2, coord.fixed = TRUE)
```

<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/PUMCH_D2_multiProb-01.png?raw=true" width="500"/>

---

### 3. Evaluation of Spatial Co-localization of Plasma Cells with *ADA2* Expression and Memory B Cells with *CD27* Expression

```r
computeGeneCelltypeOverlap(cc, "ADA2", "PC")
```
```r
#> fitting ...
#>  |================================================================================================| 100%
#>       dice p_dice   jaccard p_jaccard       mcc p_mcc
#> 1 0.3823684      0 0.2363755         0 0.1575453     0
```
```r
computeGeneCelltypeOverlap(cc, "CD27", "Bm")
```
```r
#> fitting ...
#> |================================================================================================| 100%
#> dice p_dice   jaccard p_jaccard       mcc p_mcc
#>1 0.2966482      0 0.1741556         0 0.1898258     0
```
> **Interpretation**:  *ADA2* shows spatial overlap with Plasma Cells.
> All permutation p-values are statistically significant (p < 0.001), indicating that
> *ADA2* expression is non-randomly enriched in spatial locations occupied by Plasma Cells.
> The positive MCC further supports a true spatial association, suggesting that *ADA2*⁺ PCs may represent a
> distinct spatial or functional subpopulation within the tissue.
> 
> **Interpretation**:  *CD27* exhibits spatial co-localization with Memory B cells (Bm).
> All overlap metrics are statistically significant (p < 0.001), supporting a non-random spatial relationship.
> The positive MCC further indicates that *CD27*⁺ expression is preferentially enriched in Bm-dominated regions,
> consistent with *CD27* being a canonical marker for Memory B cells.
---
### 4. Joint Binarized Spatial Domains of *ADA2*⁺ Plasma and *CD27*⁺ Memory B Cells
To investigate spatial interactions between *ADA2*⁺ plasma cells and *CD27*⁺ memory B cells, we first binarized gene expression using a Gaussian Mixture Model (GMM)-based thresholding. We then defined four spatially restricted subpopulations: *ADA2*⁺ and *ADA2*⁻ plasma cells, and *CD27*⁺ and *CD27*⁻ memory B cells. Spatial distributions of these subpopulations were visualized, and their proximity was quantified using a bidirectional nearest-neighbor distance metric to assess potential lineage relationships or spatial associations.
```r
# Binarize ADA2 expression using Gaussian Mixture Model (GMM)
ADA2_bin = as.logical(binarizeByGMM(cc@expr["ADA2",]))

# Define ADA2-positive plasma cells (ADA2⁺ PC) by intersecting ADA2⁺ and plasma cell annotations
ADA2pos_PC = as.logical(cc@meta.data$bin_PC)
ADA2pos_PC[!ADA2_bin] = FALSE

# Define ADA2-negative plasma cells (ADA2⁻ PC) by excluding ADA2⁺ from plasma cells
ADA2neg_PC = as.logical(cc@meta.data$bin_PC)
ADA2neg_PC[ADA2_bin] = FALSE

# Store binary indicators in metadata
cc@meta.data$bin_ADA2pos_PC = as.numeric(ADA2pos_PC)
cc@meta.data$bin_ADA2neg_PC = as.numeric(ADA2neg_PC)

# Repeat the same procedure for CD27 expression and memory B cells
CD27_bin = as.logical(binarizeByGMM(cc@expr["CD27",]))

# Define CD27-positive memory B cells (CD27⁺ Bm)
CD27pos_Bm = as.logical(cc@meta.data$bin_Bm)
CD27pos_Bm[!CD27_bin] = FALSE

# Define CD27-negative memory B cells (CD27⁻ Bm)
CD27neg_Bm = as.logical(cc@meta.data$bin_Bm)
CD27neg_Bm[CD27_bin] = FALSE

# Store binary indicators in metadata
cc@meta.data$bin_CD27pos_Bm = as.numeric(CD27pos_Bm)
cc@meta.data$bin_CD27neg_Bm = as.numeric(CD27neg_Bm)

# Visualize the spatial distribution of each subpopulation
plotCellType(cc, celltype = "ADA2pos_PC") + coord_fixed()
plotCellType(cc, celltype = "ADA2neg_PC") + coord_fixed()
plotCellType(cc, celltype = "CD27pos_Bm") + coord_fixed()
plotCellType(cc, celltype = "CD27neg_Bm") + coord_fixed()

```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/PUMCH_D2_ADA2_CD27cell.png?raw=true" width="1000"/>

We compared the spatial overlap between CD27⁺ memory B cells and ADA2⁺ or ADA2⁻ plasma cells.

```r
evaluate_overlap_metrics(bin1 = cc@meta.data$bin_CD27pos_Bm,
                         bin2 = cc@meta.data$bin_ADA2pos_PC,
                         coords = cc@coords)
```
```r
#>       dice p_dice  jaccard p_jaccard      mcc p_mcc
#>1 0.4753053      0 0.311738         0 0.315965     0
```
```r
evaluate_overlap_metrics(bin1 = cc@meta.data$bin_CD27pos_Bm,
                         bin2 = cc@meta.data$bin_ADA2neg_PC,
                         coords = cc@coords)
```
```r
#>       dice p_dice   jaccard p_jaccard        mcc p_mcc
#>1 0.2381575      1 0.1351753         1 -0.2510521     1
```
> **Interpretation**:  Spatial overlap metrics revealed that *CD27*⁺ memory B cells are significantly co-localized with *ADA2*⁺ plasma cells (Dice = 0.48, MCC = 0.32, p < 0.001), while showing no overlap and even spatial exclusion with *ADA2*⁻ plasma cells (MCC = –0.25, p = 1), supporting a local differentiation origin of *ADA2*⁺ plasma cells within the immune microenvironment.
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
