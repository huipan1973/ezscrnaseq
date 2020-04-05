#' Find marker genes for cell clusters
#'
#' Find marker genes for cell clusters using \pkg{scran} \code{findMarkers}
#'
#' @param clusters Vector specify cell clusters.
#' @param annot gene annotation.
#' @param block Vector specify blocking.
#' @param lfc Log fold-change.
#' @param direction Direction of change.
#' @param assay_type A string specifying which assay values to use, e.g., "counts" or "logcounts".
#' @param fdr_cutoff FDR cutoff for top marker genes.
#' @inheritParams qc_metrics
#' @inheritParams scran::findMarkers
#' @return A data.frame for the statistics and annotation of top marker geness
#' @export

find_markers <- function(sce, clusters, annot, block=NULL, lfc=1,  test.type=c("t", "wilcox", "binom"),
                          direction=c("up", "down", "any"),  pval.type=c("any", "all"), assay_type=c("logcounts", "counts", "corrected"),
						              ncores=1, fdr_cutoff=0.25, prefix=NULL, write=TRUE){

  direction <- match.arg(direction)
  test.type <- match.arg(test.type)
  pval.type <- match.arg(pval.type)
  assay_type <- match.arg(assay_type)
  stopifnot(ncores  > 0, is.numeric(lfc), is.logical(write), is.numeric(fdr_cutoff))

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(BiocParallel::bpstart(bp))
  markers <- scran::findMarkers(sce, groups=clusters, block=block, pval.type=pval.type, lfc=lfc,
                         test.type=test.type, direction=direction, assay.type=assay_type, BPPARAM=bp)
  BiocParallel::bpstop(bp)

  clus <- names(markers)
  if (pval.type=="any"){
    marker_sets <- lapply(seq_along(markers), function(i) tryCatch(data.frame(markers[[i]][markers[[i]][,3] <
			      fdr_cutoff, 1:3], Cluster=clus[i]), error=function(e) NULL))
    suffix <- "simes"
  } else if(pval.type=="all"){
    marker_sets <- lapply(seq_along(markers), function(i) tryCatch(data.frame(markers[[i]][markers[[i]][,2] <
			     fdr_cutoff, 1:2], Cluster=clus[i]), error=function(e) NULL))
    suffix <- "iut"
  }

  if (all(sapply(marker_sets, is.null))) stop (paste0("No markers have FDR < ", fdr_cutoff))

  marker_sets <- lapply(marker_sets, function(m) {m$ID <- rownames(m); m})
  marker_sets <- Reduce(rbind, marker_sets)
  marker_sets <- ezlimma::df_signif(marker_sets, 3)
  marker_sets <- data.frame(marker_sets, annot[marker_sets$ID, ])

  if (write) utils::write.csv(marker_sets, paste0(paste(c(prefix, "clusters_top_markers", suffix), collapse="_"), ".csv"),
	                       na="", row.names=FALSE)

  return(marker_sets)
}
