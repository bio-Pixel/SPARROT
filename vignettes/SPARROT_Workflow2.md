
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
library(ggplot2)

# Load Seurat spatial transcriptomics object 
seu <- readRDS("ACH0012.rds")
```
```r
# Load cell-type deconvolution matrix 
cpm <- readRDS("P9_CellProb_cell2location.rds")
view(cpm)
```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/P9_cpm.png?raw=true" width="1000"/>


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
### 2. Visualize celltype Probability

```r
p1 = plotCellTypeProb(cc, celltype = "Fibroblast")+
      coord_fixed()+
      ggtitle("Fibroblast")+ xlab(NULL)+ ylab(NULL)+
      theme_bw()+
      theme(axis.ticks = element_blank(),
            axis.text = element_blank())

p3 = plotCellTypeProb(cc, celltype = "Endothelial")+
  coord_fixed()+
  ggtitle("Endothelial")+ xlab(NULL)+ ylab(NULL)+
  theme_bw()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
library(patchwork)
p1+p2+p3

```
<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/P9_cardio_FCardioprop.png?raw=true" width="1000"/>

### 3. Visualize Multi-celltype Probabilities

```r
plotMultiCellTypeProb(cc, celltype = c("Endothelial", "Fibroblast", "Cardiomyocyte"),
                      color = c(Cardiomyocyte = "#2D81FF", Fibroblast = "#00B37F", Endothelial = "#FF6A00"),
                      outline = FALSE, coord.fixed = TRUE)
```

<img src="https://github.com/bio-Pixel/SPARROT/blob/main/vignettes/P9_cardio_prop.png?raw=true" width="500"/>

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
