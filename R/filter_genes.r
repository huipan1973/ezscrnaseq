#' Filter genes base on avarage counts
#'
#' Filter genes that have avarage counts fewer than \code{cutoff)}
#'
#' @param cutoff The cutoff of avarage counts for filtering. `cutoff` should be >= 0.
#' @inheritParams qc_metrics
#' @inheritParams scater::calcAverage
#' @return A SingleCellExperiment object.
#' @export

filter_genes <- function(sce, cutoff=0, ncores=1, prefix=NULL, plot=TRUE, write=TRUE, verbose=TRUE){

  stopifnot(cutoff >= 0, ncores > 0, is.logical(verbose), is.logical(plot), is.logical(write))
  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))
  rowData(sce)$ave.count <- scater::calculateAverage(sce, BPPARAM=bp)
  BiocParallel::bpstop(bp)

  keep <- rowData(sce)$ave.count > cutoff
  n_keep <- sum(keep)
  if (verbose) message("\nNumber of genes kept:", n_keep, "\n")

  if (plot){
    grDevices::pdf(paste(prefix, "average_counts.pdf", sep="_"))
    on.exit(grDevices::dev.off())

    graphics::hist(log10(rowData(sce)$ave.count), breaks=100, main="", col="grey", xlab=expression(Log[10]~"average count"))
    if (cutoff > min(rowData(sce)$ave.count)) graphics::abline(v=log10(cutoff), lty=2)
  }

  g2k <- as.data.frame(table(keep))
  if (write) utils::write.csv(g2k, paste(prefix, "genes_to_keep.csv", sep="_"), row.names=FALSE)

  sce <- sce[keep,, drop=FALSE]
  return(sce)
}
