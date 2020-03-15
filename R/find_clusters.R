#' Find cell clusters
#'
#' Find cell clusters using \pkg{igraph}.
#'
#' @param use_dimred A string specifying whether existing values in \code{reducedDims(sce)} should be used.
#' @param seed Random seed.
#' @param snn_k The number of nearest neighbors to consider during graph construction.
#' @param method "walktrap" or "spinglass" for finding communities in graphs via short random walks or a spin-glass model
#' and simulated annealing.
#' @param min_member Minimal number of cluster members.
#' @inheritParams qc_metrics
#' @inheritParams scran::buildSNNGraph
#' @inheritParams igraph::cluster_walktrap
#' @inheritParams igraph::cluster_spinglass
#' @return A SingleCellExperiment object with cell cluster information.
#' @export

find_clusters <- function(sce, use_dimred="PCA", seed=100, snn_k=10, ncores=1, method=c("walktrap", "spinglass"), steps=4,
                          spins=25, min_member=20, prefix=NULL, plot=TRUE, verbose=TRUE){

  method <- match.arg(method)
  stopifnot(is.logical(verbose), is.logical(plot), ncores > 0, is.numeric(seed), snn_k > 1, steps > 1, spins > 1 , min_member > 1)

  suppressWarnings(set.seed(seed = seed, sample.kind = "Rounding"))

  cl_type <- ifelse(.Platform$OS.type=="windows", "SOCK", "FORK")
  bp <- BiocParallel::SnowParam(workers=ncores, type=cl_type)
  BiocParallel::register(bpstart(bp))
  # snn
  snn_gr <- scran::buildSNNGraph(sce, use.dimred=use_dimred, k=snn_k, BPPARAM=bp)
  BiocParallel::bpstop(bp)

  # cluster
  if(method=="walktrap"){
    cluster_out <- igraph::cluster_walktrap(snn_gr, steps=steps)
  } else if(method=="spinglass"){
    cluster_out <- igraph::cluster_spinglass(snn_gr, spins=spins)
  }

  sce$Cluster <- factor(cluster_out$membership)
  sce$Cluster <- factor(paste0("clus_", sce$Cluster), levels=paste0("clus_", levels(sce$Cluster)))
  nc <- table(sce$Cluster)

  if (verbose){
    message("Clusters found:\n")
    print(kable(nc))
  }

  # modularity score
  ms <- igraph::modularity(cluster_out)
  if (verbose) message("\nModularity score: ", ms, "\n")

  # total weight between nodes
  #mod_out <- clusterModularity(snn_gr, sce$Cluster, get.values=TRUE)
  mod_out <- scran::clusterModularity(snn_gr, sce$Cluster, get.weights=TRUE)
  ratio <- log2(mod_out$observed / mod_out$expected + 1)

  if (plot){
    grDevices::pdf(paste(c(prefix, "clusters_total_weights.pdf"), collapse="_"))
    on.exit(grDevices::dev.off())
    pheatmap::pheatmap(ratio, scale="none", cluster_rows=FALSE, cluster_cols=FALSE, color=grDevices::colorRampPalette(c("white",
			"blue"))(100))
  }

  # rm small clusters
  num_small <- sum(nc < min_member)
  if(num_small > 0){
    if (verbose) message("\nRemoving ", num_small, " clusters that have cells less than ", min_member, "\n")
    sce <- sce[, sce$Cluster %in% names(nc)[nc >= min_member], drop=FALSE]
    sce$Cluster <- droplevels(sce$Cluster)
  }

  return(sce)
}
