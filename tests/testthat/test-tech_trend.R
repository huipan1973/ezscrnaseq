context("tech_trend")

test_that("default", {
  sce1 <- filter_genes(sce, cutoff=0, plot=FALSE, write=FALSE, verbose=FALSE)
  sce1 <- size_factors(sce1, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce1 <- scater::logNormCounts(sce1)
  metadata(sce1)$log.exprs.offset <- 1
  trend1 <- tech_trend(sce1, plot=FALSE)

  sce2 <- filter_genes(sce, cutoff=0, ncores=2, plot=FALSE, write=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce2, min.size=5, min.mean=0, ncores=2,plot=FALSE, verbose=FALSE)
  sce2 <- scater::logNormCounts(sce2)
  metadata(sce2)$log.exprs.offset <- 1
  trend2 <- tech_trend(sce2, ncores=2, plot=FALSE)
  expect_equal(trend1, trend2)

  #test for ncores
  trend1 <- tech_trend(sce2, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, ncores=2, plot=FALSE)
  expect_equal(trend1, trend2)
  # test for default dispersion
  trend1 <- tech_trend(sce2, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, dispersion=0, ncores=1, plot=FALSE)
  expect_equal(trend1, trend2)
  # test for assay_type
  trend1 <- tech_trend(sce2, dispersion=0, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, dispersion=0, assay_type="logcounts", ncores=1, plot=FALSE)
  expect_equal(trend1, trend2)
  # test for size.factors
  trend1 <- tech_trend(sce2, dispersion=0, assay_type="logcounts", ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, dispersion=0, size.factors=1, assay_type="logcounts", ncores=1, plot=FALSE)
  expect_equal(trend1, trend2)
  #negative test
  expect_error(tech_trend(sce2, dispersion=0, size.factors=1, assay_type="logcounts", ncores=1, plot=1))
  expect_error(tech_trend(sce2, dispersion=0, size.factors=1, assay_type="logcounts", ncores=0))
  expect_error(tech_trend(sce2, dispersion=0, size.factors=-1))
  expect_error(tech_trend(sce2, dispersion=-1))
})


test_that("truth table", {
 sc <- matrix(sample(0:200, 2000, replace=TRUE), nrow = 200, ncol = 10)
 colnames(sc) <- c(paste0("Cell_",1:10))
 row.names(sc) <- row.names(sce)[1:200]
 sc[1:100,1:5] <- 7
 sc[101:200,6:10] <- 3
 sc[1:100,6:10] <- 17
 sc[101:200,1:5] <- 13

  scVar <- sc
  itr <- round(runif(1) * 100)
  for(n in 1:itr){
    i <- sample(1:200, 1, replace=TRUE) # row number
    j <- sample(1:10, 1, replace=TRUE) #col number
    v <- sample(1:10, 1, replace=TRUE) # value to insert
    scVar[i,j] <- v
  }
  scVar1 <- SingleCellExperiment(assays = list(counts = scVar))
  sce1 <- scater::logNormCounts(scVar1)
  metadata(sce1)$log.exprs.offset <- 1
  expect_warning(trend2 <- tech_trend(sce1, plot=FALSE))
})
