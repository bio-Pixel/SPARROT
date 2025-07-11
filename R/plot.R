#' Compute concave hull of coordinate set
#'
#' @param coords Matrix of 2D coordinates
#' @param concavity Concavity parameter (lower = smoother)
#'
#' @return A data.frame of polygon coordinates
#' @export
concave_dat <- function(coords, concavity = 2) {
  df <- coords
  colnames(df) <- c("x", "y")
  hull_df <- as.data.frame(concaveman::concaveman(as.matrix(df), concavity = concavity))
  colnames(hull_df) <- c("V1", "V2")
  return(hull_df)
}

#' Plot probability of a single cell type
#'
#' @param object SparrotObj
#' @param celltype Cell type to plot
#' @param outline Whether to add concave outline
#' @param pt.size Point size
#'
#' @return A ggplot object
#' @export
plotCellTypeProb <- function(object, celltype, outline = FALSE, pt.size = 1) {
  if (!celltype %in% object@celltypes) {
    stop(paste("Cell type", celltype, "not found in object@celltypes"))
  }
  df <- as.data.frame(object@meta.data)
  df$prob <- object@cell_prob[, celltype]

  p <- ggplot(df, aes(x = row, y = col, color = prob)) +
    geom_point(size = pt.size) +
    scale_color_gradientn(colours = c("gray90","#D9AAD7","#A765B1","#A765B1")) +
    theme_void() +
    ggtitle(paste("Cell type probability:", celltype)) +
    theme(legend.position = "right")

  if (outline) {
    hull_df <- concave_dat(coords = object@coords)
    p <- p + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }
  return(p)
}

#' Plot binarized presence of a cell type
#'
#' @param object SparrotObj
#' @param celltype Cell type to plot
#' @param outline Whether to add outline
#' @param pt.size Point size
#' @param color Color for present cells
#'
#' @return A ggplot object
#' @export
plotCellType <- function(object, celltype, outline = FALSE, pt.size = 1, color = "red") {
  bin_col <- paste0("bin_", celltype)
  if (!bin_col %in% colnames(object@meta.data)) {
    stop(paste("Binarized column", bin_col, "not found in meta.data"))
  }
  df <- as.data.frame(object@meta.data)
  df$binary <- factor(df[[bin_col]], levels = c(0, 1), labels = c("absent", "present"))

  p <- ggplot(df, aes(x = row, y = col, color = binary)) +
    geom_point(size = pt.size) +
    scale_color_manual(values = c("grey80", color)) +
    theme_void() +
    ggtitle(paste("Cell type presence:", celltype)) +
    theme(legend.position = "right")

  if (outline) {
    hull_df <- concave_dat(coords = object@coords)
    p <- p + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }
  return(p)
}

#' Plot multiple cell type probabilities with normalized alpha+color blending
#'
#' @param object SparrotObj
#' @param celltype Cell types to plot
#' @param pt.size Point size
#' @param outline Whether to add outline
#' @param color Optional named color map
#' @param coord.fixed Fix aspect ratio
#' @param legend Whether to add legend
#' @param normalized Normalize probability values to 0-1
#'
#' @return A ggplot object or cowplot
#' @export
plotMultiCellTypeProb <- function(object, celltype = NULL, pt.size = 1, outline = TRUE, concavity = 2,
                                  color = NULL, coord.fixed = TRUE, legend = TRUE, normalized = TRUE) {
  prob <- object@cell_prob
  meta <- as.data.frame(object@meta.data)
  coords <- as.data.frame(object@coords)
  
  if (is.null(celltype)) {
    celltype <- colnames(prob)
  }
  if (!all(celltype %in% colnames(prob))) {
    stop("Some celltypes not found in cell_prob matrix.")
  }
  
  color_map <- if (!is.null(color)) color else object@params$celltype_colors
  base <- ggplot()  + geom_point(data = coords, aes(x = row, y = col), color = "gray99")
  
  for (ct in celltype) {
    prob_vec <- prob[, ct]
    
    if(normalized){
      rng <- range(prob_vec, na.rm = TRUE)
      if (diff(rng) == 0) {
        prob_vec <- rep(0, length(prob_vec))
      } else {
        prob_vec <- (prob_vec - rng[1]) / diff(rng)
      }
      
      threshold <- find_main_valley_threshold(prob_vec)
      valid_idx <- which(prob_vec > threshold)
    }else{
      threshold <- find_main_valley_threshold(prob_vec)
      valid_idx <- which(prob_vec > threshold)      
    }
    
    
    if (length(valid_idx) > 0) {
      prob_sel <- prob_vec[valid_idx]
      rows <- coords[valid_idx, , drop = FALSE]
      df <- data.frame(row = rows$row, col = rows$col, value = prob_sel)
      df = rbind(df, c(df$row[1], df$col[1], threshold))
      base <- base +
        geom_point(data = df, aes(x = row, y = col, alpha = value , color = value), size = pt.size) +
        scale_color_gradientn(colours = c("gray96", colorspace::lighten(color_map[ct], amount = 0.8), color_map[ct]), 
                              values = c(0,threshold,1))+
        ggnewscale::new_scale_color()
    }
  }
  
  base <- base +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
  
  if (coord.fixed) {
    base <- base + coord_fixed()
  }
  
  if (outline) {
    hull_df <- concave_dat(coords, concavity = concavity)
    base <- base + geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)
  }
  
  if (legend) {
    legend_df <- data.frame(
      row = seq_along(celltype),
      col = 1,
      celltype = factor(celltype, levels = celltype)
    )
    
    legend_layer <- ggplot(legend_df, aes(x = row, y = col, color = celltype)) +
      guides(color = guide_legend(override.aes = list(alpha = 1))) +
      geom_point(size = 3) +
      scale_color_manual(values = color_map[celltype]) +
      theme_void() +
      theme(legend.position = "right")
    
    legend_plot <- cowplot::get_legend(legend_layer)
    base <- cowplot::ggdraw() +
      cowplot::draw_plot(base, 0, 0, 0.85, 1) +
      cowplot::draw_plot(legend_plot, 0.85, 0, 0.15, 1)
  }
  return(base)
}

#' Plot spatial feature density using Nebulosa
#'
#' @param object A SparrotObj
#' @param features Vector of gene names
#' @param outline Add concave outline
#' @param raster Rasterize points
#' @param pt.size Point size
#' @param method Kernel method ("ks", "wkde")
#' @param joint Joint density
#' @param ncol Number of columns for multiple plots
#'
#' @return A ggplot object
#' @export
spFeatureDensityPlot <- function(object, features, outline = TRUE, raster = TRUE, pt.size = 1,
                         method = c("ks", "wkde"), joint = FALSE, ncol = NULL) {
  if (!requireNamespace("Nebulosa", quietly = TRUE)) {
    stop("Please install the Nebulosa package first.")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Please install the Seurat package.")
  }

  expr <- object@expr
  coords <- object@coords

  missing_genes <- setdiff(features, rownames(expr))
  if (length(missing_genes) > 0) {
    stop("These features are not found in the expression matrix: ", paste(missing_genes, collapse = ", "))
  }

  expr_mat <- as.matrix(expr[features, , drop = FALSE])
  seu <- Seurat::CreateSeuratObject(counts = CreateAssayObject(expr_mat))

  colnames(coords) = c("pixel_1","pixel_2")
  seu@reductions$pixel = Seurat::CreateDimReducObject(coords, key = "pixel")

  Seurat::DefaultAssay(seu) <- "RNA"

  if (joint){
    pp = Nebulosa::plot_density(seu, reduction = "pixel", features, 
                          size = pt.size,  raster = raster, method = method,
                          joint = joint)
    return(pp)
  }

  if (outline) {
    hull_df <- concave_dat(coords = object@coords)
    p = lapply(features,function(z){
              p = Nebulosa::plot_density(seu, reduction = "pixel", z, 
                          size = pt.size,  raster = raster, method = method)+ 
              ggplot2::coord_fixed()+
              ggplot2::scale_color_gradientn(colours = c( "white","#D9AAD7", "#A765B1"))+ 
              geom_polygon(data = hull_df, aes(x = V1, y = V2), color = "black", alpha = 0)+
              theme_bw() +
              theme(axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank())
              return(p)
          })
    pp = patchwork::wrap_plots(p, ncol = ncol)
  }else{
    p = lapply(features,function(z){
              p = Nebulosa::plot_density(seu, reduction = "pixel", z, 
                          size = pt.size,  raster = raster, method = method)+ 
              ggplot2::coord_fixed()+
              ggplot2::theme_void()+
              ggplot2::scale_color_gradientn(colours = c( "white", "#D9AAD7", "#A765B1"))+ 
              theme_bw() +
              theme(axis.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.text = element_blank())
              return(p)
          })
    pp = patchwork::wrap_plots(p, ncol = ncol)   
  }

  return(pp)
}

