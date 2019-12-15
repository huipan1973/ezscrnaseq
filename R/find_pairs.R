#' Cell cycle training step of the pair-based prediction  
#'
#' Cell cycle pairs creation using \pkg{scran} \code{sandbag}
#'
#' @param genes.list Ensembl gene IDs for genes in \code{sce}.
#' @param G1 vector of column number having G1 genes expression
#' @param S vector of column number having S genes expression
#' @param G2M vector of column number having G2M genes expression
#' @inheritParams scran::sandbag
#' @return A list containing training gene pairs for G1, S and G2M phase.
#' @export

find_pairs <- function(sce, genes.list=genes.list ,G1=G1, S=S, G2M=G2M){
  
  stopifnot(!is.null(G1), !is.null(S), !is.null(G2M), !is.null(genes.list))
  sce2<-sce[row.names(sce) %in% genes.list,]
  pairs <- sandbag(sce2, list(G1=G1, S=S, G2M=G2M))
  return(pairs)
}
