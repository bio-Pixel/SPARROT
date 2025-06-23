#' Compute pairwise overlap between cell types
#'
#' @param object A SparrotObj object
#' @param metric One of "dice", "jaccard", or "mcc"
#' @param ncore Number of cores to use (set to NULL or 1 for sequential execution)
#' @return A data.frame with celltype pairs, overlap values, and p-values
#' @export
computePairwiseCelltypeOverlap <- function(object, metric = c("dice", "jaccard", "mcc"), ncore = NULL) {
  metric <- match.arg(metric)
  ct <- object@celltypes
  combos <- t(combn(ct, 2))

  meta_bin <- as.data.frame(S4Vectors::as.data.frame(object@meta.data[, grep("^bin_", colnames(object@meta.data))]))
  coords <- object@coords

  pb <- txtProgressBar(min = 0, max = nrow(combos), style = 3)

  do_one <- function(pair) {
    bin1 <- as.logical(meta_bin[[paste0("bin_", pair[1])]])
    bin2 <- as.logical(meta_bin[[paste0("bin_", pair[2])]])
    res <- evaluate_overlap_metrics(
      bin1, bin2,
      coords = coords,
      expand_bin_dist = 1
    )
    setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
    data.frame(celltype1 = pair[1], celltype2 = pair[2],
               value = res[[metric]],
               pvalue = res[[paste0("p_", metric)]],
               stringsAsFactors = FALSE)
  }

  if (!is.null(ncore) && ncore > 1) {
    cl <- parallel::makeCluster(ncore)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, varlist = c("meta_bin", "coords", "metric", "evaluate_overlap_metrics", "expand_bin_chebyshev"), envir = environment())
    results <- parallel::parApply(cl, combos, 1, do_one)
  } else {
    results <- apply(combos, 1, do_one)
  }

  close(pb)
  df <- do.call(rbind, results)
  return(df)
}


#' Evaluate overlap metrics (Dice, Jaccard, MCC) with permutation test
#'
#' @param bin1 Logical vector
#' @param bin2 Logical vector
#' @param coords Optional coordinate matrix for Chebyshev expansion
#' @param bin_size_um Bin size in microns (optional, for future)
#' @param expand_bin_dist Expansion distance in bins
#' @param n_perm Number of permutations
#' @param seed Random seed
#' @param perm_index_mat Optional matrix of precomputed permutation indices (rows = elements, cols = permutations)
#' @return Data frame with Dice, Jaccard, MCC and permutation p-values
#' @export
evaluate_overlap_metrics <- function(bin1, bin2, coords = NULL, bin_size_um = 50,
                                     expand_bin_dist = 1, n_perm = 1000, seed = 123,
                                     perm_index_mat = NULL) {
  set.seed(seed)
  stopifnot(length(bin1) == length(bin2))
  bin1 <- as.logical(bin1)
  bin2 <- as.logical(bin2)

  if (!is.null(coords) && expand_bin_dist > 0) {
    bin1 <- expand_bin_chebyshev(bin1, coords, max_bin_dist = expand_bin_dist)
    bin2 <- expand_bin_chebyshev(bin2, coords, max_bin_dist = expand_bin_dist)
  }

  TP <- sum(bin1 & bin2)
  FP <- sum(!bin1 & bin2)
  FN <- sum(bin1 & !bin2)
  TN <- sum(!bin1 & !bin2)

  dice <- if ((2 * TP + FP + FN) > 0) 2 * TP / (2 * TP + FP + FN) else NA
  jaccard <- if ((TP + FP + FN) > 0) TP / (TP + FP + FN) else NA
  denom <- sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN)
  mcc <- if (denom > 0) ((TP * TN) - (FP * FN)) / denom else NA

  # Permutations
  if (is.null(perm_index_mat)) {
    perm_index_mat <- replicate(n_perm, sample(length(bin2)))
  }

  perm_stats <- apply(perm_index_mat, 2, function(idx) {
    bin2_perm <- bin2[idx]
    TP_p <- sum(bin1 & bin2_perm)
    FP_p <- sum(!bin1 & bin2_perm)
    FN_p <- sum(bin1 & !bin2_perm)
    TN_p <- sum(!bin1 & !bin2_perm)

    d <- if ((2 * TP_p + FP_p + FN_p) > 0) 2 * TP_p / (2 * TP_p + FP_p + FN_p) else NA
    j <- if ((TP_p + FP_p + FN_p) > 0) TP_p / (TP_p + FP_p + FN_p) else NA
    den <- sqrt(TP_p + FP_p) * sqrt(TP_p + FN_p) * sqrt(TN_p + FP_p) * sqrt(TN_p + FN_p)
    m <- if (den > 0) ((TP_p * TN_p) - (FP_p * FN_p)) / den else NA
    c(d, j, m)
  })

  p_dice <- mean(perm_stats[1, ] >= dice, na.rm = TRUE)
  p_jaccard <- mean(perm_stats[2, ] >= jaccard, na.rm = TRUE)
  p_mcc <- mean(perm_stats[3, ] >= mcc, na.rm = TRUE)

  data.frame(
    dice = dice, p_dice = p_dice,
    jaccard = jaccard, p_jaccard = p_jaccard,
    mcc = mcc, p_mcc = p_mcc
  )
}


evaluate_bidirectional_nn_distance <- function(bin1, bin2, coords, bin_size_um = 50) {
  stopifnot(length(bin1) == length(bin2))
  bin1 <- as.logical(bin1)
  bin2 <- as.logical(bin2)

  if (!is.matrix(coords)) coords <- as.matrix(coords)
  if (nrow(coords) != length(bin1)) stop("Coords row number must match bin vectors.")

  idx1 <- which(bin1)
  idx2 <- which(bin2)

  if (length(idx1) == 0 || length(idx2) == 0) {
    return(data.frame(
      mean_nn1to2_bin = NA, sd_nn1to2_bin = NA, se_nn1to2_bin = NA,
      mean_nn2to1_bin = NA, sd_nn2to1_bin = NA, se_nn2to1_bin = NA,
      mean_nn1to2_um = NA, sd_nn1to2_um = NA, se_nn1to2_um = NA,
      mean_nn2to1_um = NA, sd_nn2to1_um = NA, se_nn2to1_um = NA
    ))
  }

  dmat <- fields::rdist(coords[idx1, , drop = FALSE], coords[idx2, , drop = FALSE])

  nn1to2 <- apply(dmat, 1, min)
  nn2to1 <- apply(dmat, 2, min)

  mean_nn1to2 <- mean(nn1to2)
  sd_nn1to2 <- sd(nn1to2)
  se_nn1to2 <- sd_nn1to2 / sqrt(length(nn1to2))

  mean_nn2to1 <- mean(nn2to1)
  sd_nn2to1 <- sd(nn2to1)
  se_nn2to1 <- sd_nn2to1 / sqrt(length(nn2to1))

  data.frame(
    mean_nn1to2_bin = mean_nn1to2, sd_nn1to2_bin = sd_nn1to2, se_nn1to2_bin = se_nn1to2,
    mean_nn2to1_bin = mean_nn2to1, sd_nn2to1_bin = sd_nn2to1, se_nn2to1_bin = se_nn2to1,
    mean_nn1to2_um = mean_nn1to2 * bin_size_um, sd_nn1to2_um = sd_nn1to2 * bin_size_um, se_nn1to2_um = se_nn1to2 * bin_size_um,
    mean_nn2to1_um = mean_nn2to1 * bin_size_um, sd_nn2to1_um = sd_nn2to1 * bin_size_um, se_nn2to1_um = se_nn2to1 * bin_size_um
  )
}

#' Compute overlap between gene expression and cell type presence
#'
#' @param object A SparrotObj object
#' @param gene A gene symbol (character)
#' @param celltype A cell type name matching object@celltypes
#' @return A data frame with Dice, Jaccard, and MCC metrics
#' @export
computeGeneCelltypeOverlap <- function(object, gene, celltype) {
  if (!gene %in% rownames(object@expr)) {
    stop("Gene not found in expression matrix: ", gene)
  }
  if (!celltype %in% object@celltypes) {
    stop("Cell type not found in object: ", celltype)
  }

  # Extract gene expression and binarize
  gene_expr <- object@expr[gene, ]
  gene_bin <- binarizeByGMM(gene_expr)

  # Extract cell type binary vector
  bin_col <- paste0("bin_", celltype)
  if (!bin_col %in% colnames(object@meta.data)) {
    stop("Binarized cell type column not found in metadata: ", bin_col)
  }
  cell_bin <- as.logical(object@meta.data[[bin_col]])

  # Compute overlap metrics
  res <- evaluate_overlap_metrics(gene_bin, cell_bin, coords = object@coords, expand_bin_dist = 0)
  return(res)
}

  
