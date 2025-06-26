#' Binarize numeric vector using GMM
#'
#' @param score_vector A numeric vector
#'
#' @return Integer vector (0/1)
#' @export
binarizeByGMM <- function(score_vector) {
  if (!is.numeric(score_vector)) {
    stop("Input must be a numeric vector")
  }

  model <- mclust::Mclust(score_vector, G = 2)
  labels <- model$classification
  group_means <- tapply(score_vector, labels, mean)
  positive_label <- which.max(group_means)
  binary_label <- as.integer(labels == positive_label)
  names(binary_label) <- names(score_vector)
  return(binary_label)
}

#' Find threshold between main valley in density
#'
#' @param score_vector Numeric vector
#'
#' @return Numeric threshold
#' @export
find_main_valley_threshold <- function(score_vector) {
  dens <- density(score_vector)
  y_vals <- dens$y
  x_vals <- dens$x

  peaks_idx <- which(diff(sign(diff(y_vals))) == -2) + 1
  valleys_idx <- which(diff(sign(diff(y_vals))) == 2) + 1

  if (length(peaks_idx) < 2) {
    warning("Less than two peaks detected, returning median as threshold")
    return(median(score_vector))
  }

  peak_heights <- y_vals[peaks_idx]
  top_peaks_idx <- sort(peaks_idx[order(peak_heights, decreasing = TRUE)[1:2]])

  valleys_between <- valleys_idx[valleys_idx > top_peaks_idx[1] & valleys_idx < top_peaks_idx[2]]

  if (length(valleys_between) == 0) {
    warning("No valley between main peaks, returning median as threshold")
    return(median(score_vector))
  }

  valley_heights <- y_vals[valleys_between]
  main_valley_idx <- valleys_between[which.min(valley_heights)]

  threshold <- x_vals[main_valley_idx]
  return(threshold)
}

#' Classify by main valley threshold
#'
#' @param score_vector Numeric vector
#' @param plot Whether to show density plot
#'
#' @return Integer vector (0/1)
#' @export
classify_by_main_valley <- function(score_vector, plot = FALSE) {
  threshold <- find_main_valley_threshold(score_vector)
  labels <- as.integer(score_vector > threshold)
  names(labels) <- names(score_vector)

  if (plot) {
    dens <- density(score_vector)
    plot(dens, main = "Density with main valley threshold", xlab = "Score")
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
  }

  return(labels)
}
