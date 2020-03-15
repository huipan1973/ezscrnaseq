#' Cell cycle phase classification
#'
#' Cell cycle phase classification using \pkg{scran} \code{cyclone}.
#'
#' @param organism `hsa` or `mmu` for human or mouse genes pre-trained marker sets.
#' @param gene.names Ensembl gene IDs for genes in \code{sce}. \code{length(gene.names)} should match with
#' \code{nrow(sce)}.
#' @param pairs Pair information for training marker sets. If set to `NULL`, uses pre-trained marker sets. Some genes
#'  in pairs should overlap with \code{row.names(sce)}.
#' @param seed Random seed.
#' @inheritParams qc_metrics
#' @inheritParams scran::cyclone
#' @return A list containing phases, scores, and normalized.scores.
#' @export

ezcyclone <- function(sce, organism=c("hsa", "mmu"), gene.names=NULL, pairs=NULL, ncores=1, seed=100, iter=1000,
                       min.iter=100, min.pairs=50, verbose=TRUE){

  organism <- match.arg(organism)

  stopifnot(nrow(pairs$G1) > 0, nrow(pairs$S) > 0, nrow(pairs$G2M) > 0, ncores > 0, iter >= min.iter, min.pairs >0,
		is.logical(verbose), is.numeric(seed))

  if (is.null(gene.names)){
      gene.names=rownames(sce)
   } else {
      stopifnot(length(intersect(rownames(sce), gene.names)) == length(gene.names), nrow(sce) == length(gene.names))
   }

  if (is.null(pairs)){
    if (organism=="hsa"){
      pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    } else if(organism=="mmu"){
      pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
    }
  }

  suppressMessages(gene.pairs <- reshape2::melt(pairs))
  genes <- union(gene.pairs$first, gene.pairs$second)
  stopifnot(length(intersect(genes, row.names(sce))) > 0)

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))
  suppressWarnings(set.seed(seed = seed, sample.kind = "Rounding"))

  assignments <- scran::cyclone(sce, pairs=pairs, gene.names=gene.names, iter=iter, min.iter=min.iter, min.pairs=min.pairs,
                         BPPARAM=bp, verbose=verbose)
  BiocParallel::bpstop(bp)
  stopifnot(any(!is.na(assignments$normalized.scores)), any(!is.na(assignments$phases)), any(!is.na(assignments$scores)))
  return(assignments)
}
