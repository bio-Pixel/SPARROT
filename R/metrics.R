
evaluate_overlap_metrics <- function(bin1, bin2, coords = NULL, bin_size_um = 50, expand_bin_dist = 1,
                                     n_perm = 1000, seed = 123) {
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

  perm_dice <- numeric(n_perm)
  perm_jaccard <- numeric(n_perm)
  perm_mcc <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    bin2_perm <- sample(bin2)
    TP_p <- sum(bin1 & bin2_perm)
    FP_p <- sum(!bin1 & bin2_perm)
    FN_p <- sum(bin1 & !bin2_perm)
    TN_p <- sum(!bin1 & !bin2_perm)

    perm_dice[i] <- if ((2 * TP_p + FP_p + FN_p) > 0) 2 * TP_p / (2 * TP_p + FP_p + FN_p) else NA
    perm_jaccard[i] <- if ((TP_p + FP_p + FN_p) > 0) TP_p / (TP_p + FP_p + FN_p) else NA
    denom_p <- sqrt(TP_p + FP_p) * sqrt(TP_p + FN_p) * sqrt(TN_p + FP_p) * sqrt(TN_p + FN_p)
    perm_mcc[i] <- if (denom_p > 0) ((TP_p * TN_p) - (FP_p * FN_p)) / denom_p else NA
  }

  p_dice <- mean(perm_dice >= dice, na.rm = TRUE)
  p_jaccard <- mean(perm_jaccard >= jaccard, na.rm = TRUE)
  p_mcc <- mean(perm_mcc >= mcc, na.rm = TRUE)

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
