# SPARROT R package documentation help file snippets
# You can include these in the corresponding .R files or roxygen-style headers for each function

#' Create a SparrotObj object
#'
#' @description Create a new SparrotObj containing coordinates, expression, cell type probabilities, and metadata.
#'
#' @param coords A matrix of spatial coordinates (row, col).
#' @param expr A gene expression matrix (dgCMatrix).
#' @param cell_prob A numeric matrix of cell type probabilities.
#' @param meta.data Optional metadata as a DataFrame.
#' @param params A list of additional parameters.
#'
#' @return An object of class SparrotObj
#' @export
createSparrotObj <- function(coords, expr, cell_prob, meta.data = NULL, params = list()) {
  coords <- as.matrix(coords)
  colnames(coords) <- c("row", "col")
  expr <- methods::as(expr, "dgCMatrix")
  cell_prob <- as.matrix(cell_prob)

  if (!all(rownames(coords) == colnames(expr))) {
    stop("Row names of coords and column names of expr must be identical and in the same order.")
  }
  if (!all(rownames(coords) == rownames(cell_prob))) {
    stop("Row names of coords and cell_prob must be identical and in the same order.")
  }

  if (is.null(meta.data)) {
    meta.data <- S4Vectors::DataFrame(row.names = rownames(coords))
  } else if (!methods::is(meta.data, "DataFrame")) {
    meta.data <- S4Vectors::DataFrame(meta.data)
  }

  binary_matrix <- apply(cell_prob, 2, classify_by_main_valley)
  colnames(binary_matrix) <- paste0("bin_", colnames(binary_matrix))

  meta.data <- S4Vectors::cbind(
    S4Vectors::DataFrame(coords),
    meta.data,
    S4Vectors::DataFrame(binary_matrix)
  )

  color_50 <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#4a6fe3",
                "#b5bbe3", "#bec1d4", "#d6bcc0", "#bb7784", "#8dd593", "#f0b98d", "#f6c4e1",
                "#023fa5", "#7D58B9", "#11c638", "#F0E442", "#ef9708", "#ead3c6", "#FEC260",
                "#8e063b", "#d33f6a", "#e6afb9", "#0d6c0d", "#0fcfc0", "#7382BC", "#e07b91",
                "#F7C394", "#ade87c", "#2D81FF", "#FF6A00", "#00B37F", "#e07b91", "#A259FF",
                "#0099C6", "#FF33CC", "#66CC66", "#FFB347", "#9933FF", "#CC3366", "#33CCCC",
                "#FF6666", "#66B2FF", "#228B22", "#999933", "#CC99FF", "#FF99CC", "#669999", "#FFCC99", "#0066CC")

  n_ct <- length(colnames(cell_prob))
  if (n_ct > length(color_50)) {
    extra_colors <- grDevices::rainbow(n_ct - length(color_50))
    all_colors <- c(color_50, extra_colors)
  } else {
    all_colors <- color_50
  }
  celltype_colors <- setNames(all_colors[seq_len(n_ct)], colnames(cell_prob))

  params <- c(params, list(celltype_colors = celltype_colors))

  new("SparrotObj",
      coords = coords,
      expr = expr,
      cell_prob = cell_prob,
      meta.data = meta.data,
      celltypes = colnames(cell_prob),
      params = params)
}

#' Convert a Seurat object to a SparrotObj
#'
#' @param seurat_obj A Seurat object with spatial coordinates
#' @param cell_prob A matrix of cell type probabilities
#' @param expr_slot Expression slot to extract (e.g., "data")
#' @param assay Assay name (defaults to DefaultAssay())
#' @param params A list of additional parameters
#'
#' @return A SparrotObj
#' @export
convertSeuratToSparrot <- function(seurat_obj, cell_prob, expr_slot = "data", assay = NULL, params = list()) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.")
  }

  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  expr <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = expr_slot)
  coords <- Seurat::GetTissueCoordinates(seurat_obj)
  coords <- as.matrix(coords)
  colnames(coords) <- c("row", "col")

  if (!all(rownames(cell_prob) %in% rownames(coords))) {
    stop("Row names of cell_prob must match those of Seurat spatial coordinates.")
  }

  cell_prob <- cell_prob[rownames(coords), , drop = FALSE]

  createSparrotObj(coords = coords,
                   expr = expr,
                   cell_prob = cell_prob,
                   meta.data = seurat_obj@meta.data,
                   params = params)
}
