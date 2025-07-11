% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/expand.R
\name{plot_expand_chebyshev}
\alias{plot_expand_chebyshev}
\title{Plot expanded Chebyshev region}
\usage{
plot_expand_chebyshev(
  bin_vec,
  coords,
  max_bin_dist = 1,
  pt.size = 3,
  coord.fixed = TRUE
)
}
\arguments{
\item{bin_vec}{Logical vector of binary labels}

\item{coords}{Matrix of coordinates}

\item{max_bin_dist}{Expansion distance in bins}

\item{pt.size}{Point size}

\item{coord.fixed}{Whether to fix aspect ratio}
}
\value{
A ggplot object
}
\description{
Plot expanded Chebyshev region
}
