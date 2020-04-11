#' Make a technical trend & decompose the gene-level variance
#'
#' Make a technical trend using \pkg{scran} \code{makeTechTrend} and decompose the gene-level variance using \pkg{scran} \code{modelGeneVar}.
#'
#' @param assay_type A string specifying which assay values to use, e.g., "counts" or "logcounts".
#' @inheritParams qc_metrics
#' @inheritParams scran::modelGeneVarByPoisson
#' @return A list of a dataframe for the variance components and a function accepting a mean log-expression as input and returning the variance of tecnical noise.
#' @export

tech_trend <- function(sce, dispersion=0, assay_type="logcounts", block=NULL, design=NULL, ncores=1, prefix=NULL, plot=TRUE){

  stopifnot(ncores > 0, dispersion >=0, is.logical(plot))

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(BiocParallel::bpstart(bp))
  var_tech <- scran::modelGeneVarByPoisson(x=sce, dispersion=dispersion, assay.type=assay_type, block=block, design=design, BPPARAM=bp)
  BiocParallel::bpstop(bp)

  var_tech_trend <- S4Vectors::metadata(var_tech)$trend

  if (plot){
    grDevices::pdf(paste(c(prefix, "mean_variance_trend.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    graphics::plot(var_tech$mean, var_tech$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
    graphics::curve(var_tech_trend(x), col="dodgerblue", add=TRUE, lwd=2)
    graphics::legend("topright", legend="Technical noise", lty=1, lwd=2, col="dodgerblue")
  }

  return(list(dec=as.data.frame(var_tech), trend=var_tech_trend))
}
