
expand_bin_chebyshev <- function(bin_vec, coords, max_bin_dist = 1) {
  stopifnot(length(bin_vec) == nrow(coords))

  true_idx <- which(bin_vec)
  if(length(true_idx) == 0) return(rep(FALSE, length(bin_vec)))

  true_coords <- coords[true_idx, , drop = FALSE]

  chebyshev_dist <- function(mat1, mat2) {
    n <- nrow(mat1)
    m <- nrow(mat2)
    dist_mat <- matrix(0, n, m)
    for(i in 1:n) {
      dist_mat[i, ] <- apply(abs(sweep(mat2, 2, mat1[i,], FUN = "-")), 1, max)
    }
    dist_mat
  }

  dist_mat <- chebyshev_dist(coords, true_coords)
  min_dist_to_true <- apply(dist_mat, 1, min)
  expanded_bin <- min_dist_to_true <= max_bin_dist
  return(expanded_bin)
}

plot_expand_chebyshev <- function(bin_vec, coords, max_bin_dist = 1, pt.size = 3, coord.fixed = TRUE) {
  expanded_bin <- expand_bin_chebyshev(bin_vec, coords, max_bin_dist)
  df <- data.frame(
    x = coords[,1],
    y = coords[,2],
    original = bin_vec,
    expanded = expanded_bin
  )
  df$category <- factor(
    ifelse(df$original, "Original TRUE",
           ifelse(df$expanded & !df$original, "Expanded TRUE", "FALSE")),
    levels = c("Original TRUE", "Expanded TRUE", "FALSE")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = category)) +
    ggplot2::geom_point(size = pt.size) +
    ggplot2::scale_color_manual(values = c("blue", "red", "grey80")) +
    ggplot2::labs(title = paste("Bin expansion by Chebyshev distance =", max_bin_dist),
                  color = "Category") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank()
    )

  if (coord.fixed) {
    p <- p + ggplot2::coord_fixed()
  }

  return(p)
}
