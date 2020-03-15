#' Calculate deconvolving size factors from cell pools
#'
#' Calculate deconvolving size factors from cell pools using \pkg{ scran} \code{computeSumFactors}
#'
#' @param group.col column name of the grouping variable on the samples in colData(sce).
#' @param max.size Maximum cluster size.
#' @param seed Random seed.
#' @inheritParams qc_metrics
#' @inheritParams scran::quickCluster
#' @inheritParams scran::computeSumFactors
#' @return A SingleCellExperiment object with size factors.
#' @export

size_factors <- function(sce, min.size=10, max.size=3000, min.mean=0.1, method="igraph",  group.col=NULL,
  				seed=100, ncores=1, prefix=NULL, plot=TRUE, verbose=TRUE){

  #method <- match.arg(method)
  stopifnot(min.size <= max.size , ncores > 0, min.mean >=0, is.numeric(seed), is.logical(verbose), is.logical(plot))
  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))

  suppressWarnings(set.seed(seed = seed, sample.kind = "Rounding"))
  suppressWarnings(clusters <- scran::quickCluster(sce, min.size=min.size, min.mean=min.mean, method=method, 
			BPPARAM=bp))
  #warning message from quickCluster
  #Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE, :
  #You're computing too large a percentage of total singular values, use a standard svd instead.

  if (verbose){
    message("\nNumber of cells in clusters:\n")
    print(kable(table(clusters)))
  }

  sce <- scran::computeSumFactors(sce, cluster=clusters, max.cluster.size=max.size, min.mean=min.mean, positive=TRUE, 
					BPPARAM=bp)
  BiocParallel::bpstop(bp)

  if (verbose){
    message("\nSummary of size factors:\n")
    print(summary(SingleCellExperiment::sizeFactors(sce)))
  }


  #rm zero size factors
  if (any(SingleCellExperiment::sizeFactors(sce)==0)){
    if (verbose) message("\nRemoving", sum(SingleCellExperiment::sizeFactors(sce)==0), "cells that have size factor of zero.\n")
    sce <- sce[, SingleCellExperiment::sizeFactors(sce) > 0, drop=FALSE]
  }

  # scatter not working sce$total_counts is always NULL may be decrypted
  if (plot){
    if (is.null(sce$total_counts)) stop("Total counts are not calculated yet.")
   # @param group.col column name of the grouping variable on the samples in colData(sce). Not working
    grDevices::pdf(paste(c(prefix, "size_factor_scatter.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    if (!is.null(group.col)){
      cols <- as.numeric(as.factor(colData(sce)[, group.col]))+1
      if(max(cols) > 8) stop("Group number can't be more than 8.")
      graphics::plot(SingleCellExperiment::sizeFactors(sce), sce$total_counts/1e3, log="xy", ylab="Library size (thousands)", 
			xlab="Size factor", main="Size factors from deconvolution", 
			col=scales::alpha(grDevices::palette()[cols], 0.3), pch=16)
      legd <- sort(unique(colData(sce)[, group.col]))
      graphics::legend("topleft", col=2:(length(legd)+1), pch=16, cex=1.2, legend=legd)
    } else {
      graphics::plot(SingleCellExperiment::sizeFactors(sce), sce$total_counts/1e3, log="xy", ylab="Library size (thousands)", 
			xlab="Size factor", main="Size factors from deconvolution", 
			col= scales::alpha(grDevices::palette()[1], 0.3), pch=16)
    }
  }
  return(sce)
}
