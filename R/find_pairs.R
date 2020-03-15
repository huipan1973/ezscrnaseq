#' Cell cycle training step of the pair-based prediction
#'
#' Cell cycle pairs creation using \pkg{scran} \code{sandbag}.
#'
#' @param genes Ensembl gene IDs for genes in \code{sce}. If \code{genes} are set to `NULL` , genes will be taken
#' from \code{sce} .
#' @param G1 Vector of column indices of G1 genes.
#' @param S Vector of column indices of S genes.
#' @param G2M Vector of column indices of G2M genes.
#' @inheritParams qc_metrics
#' @inheritParams scran::sandbag
#' @details Genes in all three phases are required.
#'
#' Minimum 1000 genes should overlap between \code{sce} and \code{genes}.
#'
#' Function throws error if no pairs are found.
#' @return A list containing training gene pairs for \code{G1}, \code{S}, and \code{G2M} phase.
#' @export

find_pairs <- function(sce, G1, S, G2M, genes=NULL){
  if (is.null(genes)) genes <- rownames(sce)
  stopifnot(is.numeric(G1), is.numeric(G2M), is.numeric(S), !is.null(rownames(sce)),
            length(intersect(rownames(sce), genes)) > 1000)

  sce2 <- sce[which(rownames(sce) %in% genes),, drop=FALSE]
  pairs <- scran::sandbag(sce2, list(G1=G1, S=S, G2M=G2M))
  stopifnot(nrow(pairs$G1) > 0, nrow(pairs$S) > 0, nrow(pairs$G2M) > 0)
  return(pairs)
}
