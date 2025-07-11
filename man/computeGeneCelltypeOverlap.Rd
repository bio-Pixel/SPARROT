% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics.R
\name{computeGeneCelltypeOverlap}
\alias{computeGeneCelltypeOverlap}
\title{Compute overlap between gene expression and cell type presence}
\usage{
computeGeneCelltypeOverlap(object, gene, celltype)
}
\arguments{
\item{object}{A SparrotObj}

\item{gene}{Gene name}

\item{celltype}{Cell type name}
}
\value{
Data frame with Dice, Jaccard, and MCC values
}
\description{
Compute overlap between gene expression and cell type presence
}
