#' Quality control on the cells
#'
#' Quality control on the cells i.e. filter cells by library sizes, number of expressed genes and mitochondrial gene
#' proportion, visualize the QC metrics by histograms.
#'
#' @param sce A SingleCellExperiment object containing expression values, usually counts.
#' @param sym_col The column name for the gene symbols in \code{rowData(sce)}.
#' @param by_nmads TRUE/FALSE; should number of median absolute deviation be used as thresholds?
#' @param thresholds Numeric vector of length 3. If \code{by_nmads} is `FALSE`, the elements represent thresholds for
#' total number of reads, number of detected genes, percentage of mitochondrial genes. If \code{by_nmads} is `TRUE`,
#' thresholds for all 3 elements is calculated in terms of number of median absolute deviations, which should not exceed 5.
#' @param ncores Number of cores.
#' @param prefix Prefix for file name for the QC metrics histograms.
#' @param plot TRUE/FALSE for whether to plot the QC metrics histograms.
#' @param write TRUE/FALSE for whether to write the table of filtered cells.
#' @param verbose TRUE/FALSE for specifying if diagnostics should be printed to screen.
#' @import knitr
#' @import SummarizedExperiment
#' @import SingleCellExperiment
#' @import scater
#' @import scran
#' @import Hmisc
#' @import BiocParallel
#' @import pheatmap
#' @return A SingleCellExperiment object.
#' @export

qc_metrics <- function(sce, sym_col="symbol", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, prefix=NULL, plot=TRUE, write=TRUE,
			verbose=TRUE){

  stopifnot(ncores > 0, is.logical(verbose), is.logical(plot), is.logical(write))
  # specifying the mitochodial genes
  is.mito <- grepl("^(M|m)(T|t)-", rowData(sce)[, sym_col])
  n_mito <- sum(is.mito)
  if (verbose) message("Number of mitochondrial genes:", n_mito, "\n")

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))
  #sce <- calculateQCMetrics(sceMock, feature_controls=list(Mt=is.mito), BPPARAM=bp)
  #sce1 <- scater::perFeatureQCMetrics(sce, subsets=list(Mt=is.mito), BPPARAM=bp)
  sce1 <- scater::perCellQCMetrics(sce, subsets=list(Mt=is.mito), BPPARAM=bp)
  BiocParallel::bpstop(bp)

  # qc calculation
  if (by_nmads) {
    if (any(thresholds > 5)) stop("Thresholds are too big for using NMADS")
    #libsize.drop <- isOutlier(sce$total_counts, nmads=thresholds[1], type="lower", log=TRUE)
    #feature.drop <- isOutlier(sce$total_features_by_counts, nmads=thresholds[2], type="lower", log=TRUE)
    libsize.drop <- isOutlier(sce1$sum, nmads=thresholds[1], type="lower", log=TRUE)
    feature.drop <- isOutlier(sce1$detected, nmads=thresholds[2], type="lower", log=TRUE)
    if (n_mito >0) mito.drop <- suppressWarnings(isOutlier(sce1$subsets_Mt_sum, nmads=thresholds[3], type="higher"))
  } else {
    if(any(thresholds < 10)) stop("Thresholds are too small for using actual counts or percentages")
    libsize.drop <- sce1$sum < thresholds[1]
    feature.drop <- sce1$detected < thresholds[2]
    if (n_mito > 0) mito.drop <- sce1$subsets_Mt_sum > thresholds[3]
  }

  # qc tab
  if (n_mito > 0){
    qc <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop),
                     Remaining=sum(!(libsize.drop | feature.drop | mito.drop)))
    #cutoff <- data.frame(LibSize=max(sce$total_counts[libsize.drop]), Feature=max(sce$total_features_by_counts[feature.drop]),
    #                    Mito=min(sce$pct_counts_Mt[mito.drop]))
    cutoff <- data.frame(LibSize=max(sce1$subsets_Mt_sum[libsize.drop]), Feature=max(sce1$subsets_Mt_detected[feature.drop]),
                        Mito=min(sce1$subsets_Mt_sum[mito.drop]))
  } else {
    qc <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=0,
                     Remaining=sum(!(libsize.drop | feature.drop)))
    cutoff <- data.frame(LibSize=max(sce1$subsets_Mt_sum[libsize.drop]), Feature=max(sce1$subsets_Mt_detected[feature.drop]))
  }

  if (verbose){
    message("Filtered cells\n")
    print(kable(qc))
    message("\nCutoff\n")
    print(kable(cutoff))
  }

  if (write) utils::write.csv(qc, paste(c(prefix, "filtered_cells.csv"), collapse="_"), row.names=FALSE)

  # qc metric
  if (plot){
    if (n_mito > 0){
      grDevices::pdf(paste(c(prefix,"qc_metrics_hist.pdf"), collapse="_"), 9, 3)
      graphics::par(mfrow=c(1,3), mar=c(5.1, 4.1, 1.1, 0.1), oma=c(0, 0, 2, 0))

    } else{
      grDevices::pdf(paste(c(prefix,"qc_metrics_hist.pdf"), collapse="_"), 6, 3)
      graphics::par(mfrow=c(1,2), mar=c(5.1, 4.1, 1.1, 0.1), oma=c(0, 0, 2, 0))
    }
    on.exit(grDevices::dev.off())

    graphics::hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", breaks=30, col="grey80", ylab="Number of cells")
    graphics::abline(v=cutoff$LibSize[1]/1e3, lty=2)

    graphics::hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", breaks=30, col="grey80",
                   ylab="Number of cells")
    graphics::abline(v=cutoff$Feature[1], lty=2)

    if (n_mito > 0){
      graphics::hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", main="", breaks=30, col="grey80",  ylab="Number of cells")
      graphics::abline(v=cutoff$Mito[1], lty=2)
    }
    mtitle("Histograms of QC Metrics", cex.m=1.2)
  }

  if (n_mito > 0){
    keep <- !(libsize.drop | feature.drop | mito.drop)
  } else{
    keep <- !(libsize.drop | feature.drop)
  }

  sce <- sce[keep,, drop=FALSE]
  return(sce)
}
