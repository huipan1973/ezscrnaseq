#' Quality control on the cells
#'
#' Quality control on the cells i.e. filter cells by library sizes, number of expressed genes and mitochondrial gene
#' proportion, visualize the QC metrics by histograms.
#'
#' @param sce A SingleCellExperiment object containing expression values, usually counts.
#' @param sym_col The column name for the gene symbols in \code{rowData(sce)}.
#' @param by_nmads TRUE/FASLE for using number of median absolute deviation as thresholds.
#' @param thresholds Numbers of median absolute deviation if \code{by_nmads} is TRUE, otherwise the actual counts or
#' percentages.
#' @param ncores Number of cores.
#' @param prefix Prefix for file name for the QC metrics histograms.
#' @param plot TRUE/FASLE for whether plot the QC metrics histograms.
#' @param write TRUE/FASLE for whether write the table of filtered cells.
#' @param verbose TRUE/FASLE for specifying whether diagnostics should be printed to screen.
#' @return A SingleCellExperiment object.
#' @export

qc_metrics <- function(sce, sym_col="symbol", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, prefix=NULL, plot=TRUE, write=TRUE,
			verbose=TRUE){

  stopifnot(ncores > 0, is.logical(verbose), is.logical(plot), is.logical(write))
  # specifying the mitochodial genes
  is.mito <- grepl("^(M|m)(T|t)-", SummarizedExperiment::rowData(sce)[, sym_col])
  n_mito <- sum(is.mito)
  if (verbose) message("Number of mitochondrial genes:", n_mito, "\n")

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(BiocParallel::bpstart(bp))

  stats <- as.data.frame(scater::perCellQCMetrics(sce, subsets=list(Mt=is.mito), BPPARAM=bp))
  BiocParallel::bpstop(bp)

  SummarizedExperiment::colData(sce) <- cbind(SummarizedExperiment::colData(sce), stats)

  # qc calculation
  if (by_nmads) {
    if (any(thresholds > 5)) stop("Thresholds are too big for unsing NMADS")
    libsize.drop <- scater::isOutlier(stats$sum, nmads=thresholds[1], type="lower", log=TRUE)
    feature.drop <- scater::isOutlier(stats$detected, nmads=thresholds[2], type="lower", log=TRUE)
    if (n_mito >0) mito.drop <- suppressWarnings(scater::isOutlier(stats$subsets_Mt_percent, nmads=thresholds[3], type="higher"))

  } else{
    if(any(thresholds < 10)) stop("Thresholds are too small for unsing actual counts or percentages")

    libsize.drop <- stats$sum < thresholds[1]
    feature.drop <- stats$detected < thresholds[2]
    if (n_mito > 0) mito.drop <- stats$subsets_Mt_percent > thresholds[3]
  }

  # qc tab
  if (n_mito > 0){
    qc <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop),
                     Remaining=sum(!(libsize.drop | feature.drop | mito.drop)))

    cutoff <- data.frame(LibSize=max(stats$sum[libsize.drop]), Feature=max(stats$detected[feature.drop]),
                        Mito=min(stats$subsets_Mt_percent[mito.drop]))

  } else{
    qc <- data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=0,
                     Remaining=sum(!(libsize.drop | feature.drop)))
    cutoff <- data.frame(LibSize=max(stats$sum[libsize.drop]), Feature=max(stats$detected[feature.drop]))
  }

  if (verbose){
    message("Filtered cells\n")
    print(knitr::kable(qc))
    message("\nCutoff\n")
    print(knitr::kable(cutoff))
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

    graphics::hist(stats$sum/1e3, xlab="Library sizes (thousands)", main="", breaks=30, col="grey80", ylab="Number of cells")
    graphics::abline(v=cutoff$LibSize[1]/1e3, lty=2, col="red")

    graphics::hist(stats$detected, xlab="Number of expressed genes", main="", breaks=30, col="grey80", ylab="Number of cells")
    graphics::abline(v=cutoff$Feature[1], lty=2, col="red")

    if (n_mito >0){
      graphics::hist(stats$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", main="", breaks=30, col="grey80",  ylab="Number of cells")
      graphics::abline(v=cutoff$Mito[1], lty=2, col="red")
    }
    Hmisc::mtitle("Histograms of QC Metrics", cex.m=1.2)
  }

  if (n_mito > 0){
    keep <- !(libsize.drop | feature.drop | mito.drop)
  } else{
    keep <- !(libsize.drop | feature.drop)
  }

  sce <- sce[, keep]
  return(sce)
}
