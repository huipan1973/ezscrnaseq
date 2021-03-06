% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_metrics.r
\name{qc_metrics}
\alias{qc_metrics}
\title{Quality control on the cells}
\usage{
qc_metrics(
  sce,
  sym_col = "symbol",
  by_nmads = TRUE,
  thresholds = c(3, 3, 3),
  ncores = 1,
  prefix = NULL,
  plot = TRUE,
  write = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{sce}{A SingleCellExperiment object containing expression values, usually counts.}

\item{sym_col}{The column name for the gene symbols in \code{rowData(sce)}.}

\item{by_nmads}{TRUE/FASLE for using number of median absolute deviation as thresholds.}

\item{thresholds}{Numbers of median absolute deviation if \code{by_nmads} is TRUE, otherwise the actual counts or
percentages.}

\item{ncores}{Number of cores.}

\item{prefix}{Prefix for file name for the QC metrics histograms.}

\item{plot}{TRUE/FASLE for whether plot the QC metrics histograms.}

\item{write}{TRUE/FASLE for whether write the table of filtered cells.}

\item{verbose}{TRUE/FASLE for specifying whether diagnostics should be printed to screen.}
}
\value{
A SingleCellExperiment object.
}
\description{
Quality control on the cells i.e. filter cells by library sizes, number of expressed genes and mitochondrial gene
proportion, visualize the QC metrics by histograms.
}
