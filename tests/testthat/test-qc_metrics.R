context("qc_metrics")

test_that("by_nmads", {
  sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE)
  sce2 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE)
  expect_equal(sce1, sce2)
})

test_that("not_by_nmads", {
  expect_warning(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(90000,1200,50), plot=FALSE, write=FALSE,
		ncores=1, verbose=FALSE))
  expect_warning(sce2 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(90000,1200,50), plot=FALSE, write=FALSE,
		ncores=2, verbose=FALSE))
  expect_equal(sce1, sce2)
})

test_that("large_nmads", {
  expect_error(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(10,3,3), plot=FALSE, write=FALSE,
		 verbose=FALSE))
})

test_that("small counts", {
  expect_error(sce1 <- qc_metrics(sce, sym_col="Gene", by_nmads=FALSE, thresholds=c(9,600,20), plot=FALSE, write=FALSE,
		verbose=FALSE))
})

test_that("non mito genes", {
  sceNew <- sce
  rowData(sceNew)[, "Gene"] <- row.names(sce)
  sce1 <- qc_metrics(sceNew, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE)
  sce2 <- qc_metrics(sceNew, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE)
  expect_equal(sce1, sce2)
})


test_that("negative test", {
  #ncore =0
  expect_error(qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=0, plot=FALSE, write=FALSE,
		verbose=FALSE))
  # logical tests for argument
  expect_error(qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=1))
  expect_error(qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, write=1))
  expect_error(qc_metrics(sce, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, verbose=1))
})

test_that("truth table", {
  #by_nmads
  sc <- matrix(sample(0:200, 20000, replace=T), nrow = 2000, ncol = 10)
  colnames(sc) <- c(paste0("Cell_",1:10))
  row.names(sc) <- row.names(sce)
  #sc[2,]<-220
  sc[,1]<-0
  sc[,9]<-0
  #sc[42,]<-220
  rowDataSc<-rowData(sce)
  itr <- round(runif(1) * 100)
  for(n in 1:itr)
  {
    i <- sample(1:200, 1, replace=T) # row number
    j <- sample(1:10, 1, replace=T) #col number
    sc[i,j] <- 0
  }
  itr <- round(runif(1) * 100)
  for(n in 1:itr)
  {
    i <- sample(1:200, 1, replace=T) # row number
    j <- sample(1:10, 1, replace=T) #col number
    sc[i,j] <- 200
  }
  rowDataSc$BaseGeneMean <- apply(sc, 1, mean)
  rowDataSc$GeneMean <- apply(sc, 1, mean)

  scVar1 <- SingleCellExperiment(assays = list(counts = sc))
  rowData(scVar1) <- rowDataSc
  expect_warning(sce1 <- qc_metrics(scVar1, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE))
  expect_true(dim(sce1)[2] < dim(scVar1)[2])

  #not_by_nmads
  expect_warning(sce1 <- qc_metrics(scVar1, sym_col="Gene", by_nmads=FALSE, thresholds=c(198700,1980,800), plot=FALSE, write=FALSE,
		ncores=1, verbose=FALSE))
  expect_true(dim(sce1)[2] < dim(scVar1)[2])

  #non mito genes
  sceNew <- scVar1
  rowData(sceNew)[, "Gene"] <- row.names(scVar1)
  sce1 <- qc_metrics(sceNew, sym_col="Gene", by_nmads=TRUE, thresholds=c(3,3,3), ncores=1, plot=FALSE, write=FALSE,
		verbose=FALSE)
  expect_true(dim(sce1)[2] < dim(scVar1)[2])
})

