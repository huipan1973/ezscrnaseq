context("find_markers")

test_that("find_markers", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, plot=FALSE, verbose=FALSE)
  sce2 <- scater::logNormCounts(sce1)
  metadata(sce2)$log.exprs.offset <- 1
  trend <- tech_trend(sce2, plot=FALSE)
  expect_warning(set.seed(seed = 100, sample.kind = "Rounding"))
  expect_warning(sce3 <- denoisePCA(sce2, technical=trend$trend, assay.type="logcounts", max.rank=100))
  sceCluster3 <- find_clusters(sce3, snn_k=5, plot=FALSE, verbose=FALSE)

  # ncores 1 vs 2
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=1,
  		annot=data.frame(rowData(sce1),row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)

  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, ncores=2,
		 annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)
  expect_equal(clus_df1, clus_df2)
  # test.type ="t"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster,
  		annot=data.frame(rowData(sce1), row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE)
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, test.type ="t",
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
  # direction ="up"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, ncores=1, write=FALSE )
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, direction="up")
  expect_equal(clus_df1, clus_df2)
  # direction ="down"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, direction="down" )
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, direction="down")
  expect_equal(clus_df1, clus_df2)
  # pval.type ="any"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
  		row.names=rownames(sce1)), fdr_cutoff=Inf, ncores=1, write=FALSE)
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, pval.type ="any")
  expect_equal(clus_df1, clus_df2)
  # pval.type ="all"
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, pval.type ="all")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, pval.type ="all")
  expect_equal(clus_df1, clus_df2)

  # assay_type="logcounts",
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, ncores=1, write=FALSE)
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, ncores=1, write=FALSE, assay_type="logcounts")
  expect_equal(clus_df1, clus_df2)

  # assay_type="counts",
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, assay_type="counts")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		 row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=2, assay_type="counts")
  expect_equal(clus_df1, clus_df2)
  # lfc =1
  clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, assay_type="counts")
  clus_df2 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, lfc =1, annot=data.frame(rowData(sce1),
		 row.names=rownames(sce1)), fdr_cutoff=Inf, write=FALSE, ncores=1, assay_type="counts")
  expect_equal(clus_df1, clus_df2)

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
  #negative test  assay_type other
  expect_error(clus_df1 <- find_markers(sceCluster3, clusters=sceCluster3$Cluster, annot=data.frame(rowData(sce1),
		row.names=rownames(sce1)), fdr_cutoff=Inf, assay_type="MyAssay", write=FALSE))

})
