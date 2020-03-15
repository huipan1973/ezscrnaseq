context("find_markers")

test_that("find_markers", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE)
  sce2 <- scater::logNormCounts(sce1)
  metadata(sce2)$log.exprs.offset <- 1
  trend <- tech_trend(sce2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce3 <- denoisePCA(sce2, technical=trend, assay.type="logcounts", max.rank=100))
  sceCluster3 <- find_clusters(sce3, snn_k=5, plot=FALSE, verbose=FALSE)

  # ncores 1 vs 2
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=1,
  		annot=data.frame(rowData(sce1),row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)

  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=2,
		 annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)
  expect_equal(clus_df1, clus_df2)

  # test.type ="wilcox"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, test.type ="wilcox")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, test.type ="wilcox")
  expect_equal(clus_df1, clus_df2)

  # test.type ="binom"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, test.type ="binom")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, test.type ="binom")
  expect_equal(clus_df1, clus_df2)

  # direction ="down"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, direction="down" )
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, direction="down")
  expect_equal(clus_df1, clus_df2)

  # pval.type ="all"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, pval.type ="all")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, pval.type ="all")
  expect_equal(clus_df1, clus_df2)

  # assay_type="counts",
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, assay_type="counts")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		 row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, assay_type="counts")
  expect_equal(clus_df1, clus_df2)

  # default vs user provided test
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=1,
  		annot=data.frame(rowData(sce1),row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)

  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=2,
		 annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, lfc =1,
		 direction ="up", pval.type ="any", assay_type="logcounts", test.type ="t")
  expect_equal(clus_df1, clus_df2)
  # test.type ="t"

  # fdr_cutoff=0.25
  expect_error(find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), write=FALSE, ncores=1, assay_type="counts"))
  expect_error(find_markers(sceCluster3, clusters=sceCluster3$Cluster, fdr_cutoff=0.25, annot=data.frame(rowData(sce1),
		 row.names=rownames(sce1)), write=FALSE, ncores=1, assay_type="counts"))

  #negative test ncores=0
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		 row.names=rownames(sce1)), fdr_cutoff=Inf, ncores=0, write=FALSE))
  #negative test write not logical
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=1))
  #negative test lfc not numeric
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, lfc="up"))
  #negative test fdr_cutoff not numeric
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff="Inf"))

  #negative test  test.type other
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, test.type="abc", write=FALSE))
  #negative test  pval.type other
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, pval.type=0.05, write=FALSE))
  #negative test direction other
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, direction="any", write=FALSE))
  #negative test  assay_type other
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, assay_type="MyAssay", write=FALSE))
})

test_that("truth table", {
  sc <- matrix(0, nrow = 200, ncol = 10)
  colnames(sc) <- c(paste0("Cell_", 1:10))
  sc[sc > 0] <- 0
  sc[1:100, 1:5] <- 3
  sc[101:200, 6:10] <- 7
  scVar <- sc
  itr <- round(runif(1) * 100)
  for(n in 1:itr){
    i <- sample(1:200, 1, replace=T) # row number
    j <- sample(1:10, 1, replace=T) #col number
    v <- sample(1:10, 1, replace=T) # value to instert
    scVar[i,j] <- v
  }
  scVar1 <- SingleCellExperiment(assays = list(counts = scVar))
  rowDataSc <- rowData(sce)[1:200,]
  rowDataSc$BaseGeneMean <- apply(sc, 1, mean)
  rowDataSc$GeneMean <- apply(sc, 1, mean)
  rowData(scVar1) <- rowDataSc
  sce1 <- scater::logNormCounts(scVar1)
  metadata(sce1)$log.exprs.offset <- 1
  expect_warning(trend <- tech_trend(sce1, ncores=2, plot=FALSE))

  seed <- sample(1:100, size=1, replace=TRUE)
  expect_warning(set.seed(seed = seed, sample.kind = "Rounding"))
  expect_warning(sce1 <- denoisePCA(sce1, technical=trend, assay.type="logcounts", max.rank=100))
  sceclust1 <- find_clusters(sce1, snn_k=2, ncores=2, plot=FALSE , min_member=2, verbose=FALSE)

  clus_df1 <- find_markers(sceclust1, clusters=sceclust1$Cluster, ncores=1,
               annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)

  expect_equal(levels(factor(clus_df1$Cluster)), levels(sceclust1$Cluster))
  expect_equal(max(clus_df1$GeneMean), max(rowDataSc$GeneMean))
  expect_equal(max(clus_df1$BaseGeneMean), max(rowDataSc$BaseGeneMean))
})
