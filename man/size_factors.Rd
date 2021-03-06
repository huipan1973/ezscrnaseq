% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/size_factors.R
\name{size_factors}
\alias{size_factors}
\title{Calculate deconvolving size factors from cell pools}
\usage{
size_factors(
  sce,
  min.size = 10,
  max.size = 3000,
  min.mean = 0.1,
  method = "igraph",
  group.col = NULL,
  seed = 100,
  ncores = 1,
  prefix = NULL,
  plot = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing expression values, usually counts.}

\item{min.size}{An integer scalar specifying the minimum size of each cluster.}

\item{max.size}{Maximum cluster size.}

\item{min.mean}{A numeric scalar specifying the filter to be applied on the average count for each filter prior to computing ranks.
Only used when \code{use.ranks=TRUE}, see \code{?\link[scran]{scaledColRanks}} for details.}

\item{method}{String specifying the clustering method to use.
\code{"hclust"} uses hierarchical clustering while \code{"igraph"} uses graph-based clustering.}

\item{group.col}{column name of the groups.}

\item{seed}{Random seed.}

\item{ncores}{Number of cores.}

\item{prefix}{Prefix for file name for the QC metrics histograms.}

\item{plot}{TRUE/FASLE for whether plot the QC metrics histograms.}

\item{verbose}{TRUE/FASLE for specifying whether diagnostics should be printed to screen.}
}
\value{
A SingleCellExperiment object with size factors.
}
\description{
Calculate deconvolving size factors from cell pools using \pkg{ scran} \code{computeSumFactors}
}
