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
  expect_equal(trend1, trend2 )

  #test for ncores
  trend1 <- tech_trend(sce2, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, ncores=2, plot=FALSE)
  expect_equal(trend1$trend, trend2$trend)
  # test for default dispersion
  trend1 <- tech_trend(sce2, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, dispersion=0, ncores=1, plot=FALSE)
  expect_equal(trend1$trend, trend2$trend)
  # test for assay_type
  trend1 <- tech_trend(sce2, dispersion=0, ncores=1, plot=FALSE)
  trend2 <- tech_trend(sce2, dispersion=0, assay_type="logcounts", ncores=1, plot=FALSE)
  expect_equal(trend1$trend, trend2$trend)
  #negative test
  expect_error(tech_trend(sce2, dispersion=0, assay_type="logcounts", ncores=1, plot=1))
  expect_error(tech_trend(sce2, dispersion=0, assay_type="logcounts", ncores=0, plot=FALSE))
  expect_error(tech_trend(sce2, dispersion=-1, plot=FALSE))
})








