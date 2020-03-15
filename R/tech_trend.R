#' Make a technical trend & decompose the gene-level variance
#'
#' Make a technical trend using \pkg{scran} \code{makeTechTrend} and decompose the gene-level variance using
#' \pkg{scran} \code{modelGeneVar}.
#'
#' @param assay_type A string specifying which assay values to use, e.g., "counts" or "logcounts".
#' @inheritParams qc_metrics
#' @inheritParams scran::makeTechTrend
#' @inheritParams scran::modelGeneVar
#' @inheritParams stats::loess
#' @return A function accepting a mean log-expression as input and returning the variance of the log-expression as the output.
#' @export

tech_trend <- function(sce, dispersion=0, assay_type="logcounts", ncores=1, size.factors=1, prefix=NULL, plot=TRUE){

  stopifnot(ncores > 0, dispersion >= 0, size.factors >= 0, is.logical(plot))
  # no spike ref: https://github.com/MarioniLab/scran/issues/7
  # according to scran::multiBlockVar(), use logcounts for tech trend
  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))
  var_fit_trend <- scran::makeTechTrend(dispersion=dispersion, x=sce, BPPARAM=bp)

  var_fit <- scran::modelGeneVar(sce, assay.type=assay_type)
  #deprecated function
  #var_fit <- trendVar(sce, parametric=FALSE, loess.args=list(span=span), use.spikes=FALSE, assay.type=assay_type)
  #var_out <- decomposeVar(sce, fit=var_fit, block=block, design=design, assay.type=assay_type, BPPARAM=bp)

  BiocParallel::bpstop(bp)

  if (plot){
    grDevices::pdf(paste(c(prefix, "mean_variance_trend.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())

    graphics::plot(var_fit$mean, var_fit$total, pch=16, cex=0.6, xlab="Mean log-expression", ylab="Variance of log-expression")
    #graphics::curve(var_fit_trend(var_fit$mean), col="dodgerblue", add=TRUE, lwd=2)
    #graphics::curve(var_fit$trend(var_fit$mean), col="red", add=TRUE, lwd=2)
    #graphics::legend("topright", legend=c("Technical noise", "All variance"), lty=1, lwd=2, col=c("dodgerblue", "red"))

   #p<- ggplot() +  xlab("Mean log-expression")+ ylab("Variance of log-expression") +
   #              geom_point(data=plot, aes(x=var_fit.mean , y=var_fit.total)) +
   #              stat_function(data=data.frame(var_fit_trend(var_fit$mean)), aes(x), fun=eq, colour="red") +
   #              stat_function(data=data.frame(var_fit1$trend(var_fit$mean)), aes(x), fun=eq, colour="dodgerblue")
  }

  return(var_fit_trend)
}
