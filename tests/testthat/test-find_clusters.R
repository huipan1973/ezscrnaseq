context("find_clusters")

test_that("ncores=2", {
  sce1 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)
  sce1 <- scater::logNormCounts(sce1)
  metadata(sce1)$log.exprs.offset <- 1
  trend <- tech_trend(sce1, ncores=2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend$trend, assay.type="logcounts", max.rank=100))

  # method=NULL
  sceclust1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE)
  sceclust2 <- find_clusters(sce1, snn_k=5, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sceclust1, sceclust2)

  # method="spinglass"
  sceclust1 <- find_clusters(sce1, snn_k=5, ncores=2, plot=FALSE, verbose=FALSE, method="spinglass")
  sceclust2 <- find_clusters(sce1, snn_k=5, ncores=1, plot=FALSE, verbose=FALSE, method="spinglass")
  expect_equal(sceclust1, sceclust2)

  # default vs user provided
  sceclust1 <- find_clusters(sce1, plot=FALSE, verbose=FALSE)
  sceclust2 <- find_clusters(sce1, snn_k=10, ncores=1, plot=FALSE, verbose=FALSE, method="walktrap", use_dimred="PCA",
                            seed=100, steps=4, min_member=20)
  expect_equal(sceclust1, sceclust2)

  # defaul spins vs user provided spins for method spinglass
  sceclust1 <- find_clusters(sce1, ncores=1, plot=FALSE, verbose=FALSE, method="spinglass")
  sceclust2 <- find_clusters(sce1, spins=25, ncores=1, plot=FALSE, verbose=FALSE, method="spinglass")
  expect_equal(sceclust1, sceclust2)


  # negative tests
  expect_error(find_clusters(sce1, ncores=1, verbose=1))
  expect_error(find_clusters(sce1, ncores=1, plot=1))
  expect_error(find_clusters(sce1, ncores=0, plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, seed="abc", plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, snn_k=0, plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, steps=0, plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, spins=0, plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, min_member=0, plot=FALSE, verbose=FALSE))
  expect_error(find_clusters(sce1, ncores=1, method="myCluster", plot=FALSE, verbose=FALSE))
})


test_that("truth table", {
  sc <- matrix(0, nrow = 200, ncol = 10)
  colnames(sc) <- c(paste0("Cell_",1:10))
  sc[sc > 0] <- 0
  sc[1:100, 1:5] <- 3
  sc[101:200, 6:10] <- 7
  scVar <- sc
  itr <- round(runif(1) * 100)
  for(n in 1:itr)
  {
    i <- sample(1:200, 1, replace=T) # row number
    j <- sample(1:10, 1, replace=T) #col number
    v <- sample(1:10, 1, replace=T) # value to instert
    scVar[i,j] <- v
  }
  scVar1 <- SingleCellExperiment(assays = list(counts = scVar))
  sce1 <- scater::logNormCounts(scVar1)
  metadata(sce1)$log.exprs.offset <- 1
  trend <- tech_trend(sce1, ncores=2, plot=FALSE)
  seed <- sample(1:100,1,replace=T)
  expect_warning(set.seed(seed = seed, sample.kind = "Rounding"))
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend$trend, assay.type="logcounts", max.rank=100))
  sceclust1 <- find_clusters(sce1, snn_k=2, ncores=2, plot=FALSE , min_member=2, verbose=FALSE)
  expect_true(all(as.character(sceclust1$Cluster) == c(rep("clus_1", times=5), rep("clus_2", times=5))) ||
  all(as.character(sceclust1$Cluster) == c(rep("clus_2", times=5), rep("clus_1", times=5))))
})
