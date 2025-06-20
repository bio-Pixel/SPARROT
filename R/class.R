
setClass("SparrotObj",
         slots = list(
           coords = "matrix",
           expr = "dgCMatrix",
           cell_prob = "matrix",
           meta.data = "DataFrame",
           celltypes = "character",
           params = "list"
         )
)

setGeneric("getCoords", function(object) standardGeneric("getCoords"))
setMethod("getCoords", "SparrotObj", function(object) object@coords)

setGeneric("getExpr", function(object) standardGeneric("getExpr"))
setMethod("getExpr", "SparrotObj", function(object) object@expr)

setGeneric("getCellProb", function(object) standardGeneric("getCellProb"))
setMethod("getCellProb", "SparrotObj", function(object) object@cell_prob)

setGeneric("getMeta", function(object) standardGeneric("getMeta"))
setMethod("getMeta", "SparrotObj", function(object) object@meta.data)

setGeneric("getCelltypes", function(object) standardGeneric("getCelltypes"))
setMethod("getCelltypes", "SparrotObj", function(object) object@celltypes)

setGeneric("getParams", function(object) standardGeneric("getParams"))
setMethod("getParams", "SparrotObj", function(object) object@params)

setMethod("show", "SparrotObj", function(object) {
  cat("An object of class 'SparrotObj'\n\n")
  cat("Number of spots/cells: ", nrow(object@coords), "\n")
  cat("Number of genes:       ", nrow(object@expr), "\n")
  cat("Number of cell types:  ", length(object@celltypes), "\n")

  meta_cols <- colnames(object@meta.data)
  coord_cols <- colnames(object@coords)
  bin_cols <- grep("^bin_", meta_cols, value = TRUE)
  other_cols <- setdiff(setdiff(meta_cols, coord_cols), bin_cols)

  cat("Meta data columns:\n")
  cat("  • coords:    ", paste(coord_cols, collapse = ", "), "\n")
  if (length(bin_cols)) {
    cat("  • binarized: ", paste(bin_cols, collapse = ", "), "\n")
  }
  if (length(other_cols)) {
    cat("  • others:    ", paste(other_cols, collapse = ", "), "\n")
  }

  cat("\nUse @meta.data, @expr, @cell_prob, or accessor methods to explore.\n")
})
