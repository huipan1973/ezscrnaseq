context("size_factors")

test_that("default vs user provided args", {
  #test for ncore
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for min.size = 10
  sce1 <- size_factors(sce, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=10, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   max.size = 3000
  sce1 <- size_factors(sce, min.size=5, min.mean=0, max.size = 3000, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   min.mean = 0.1
  sce1 <- size_factors(sce, min.size=5, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0.1 , ncores=1, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   seed = 100
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, seed = 100, plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #test for   method = igraph
  sce1 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, plot=FALSE, verbose=FALSE)
  sce2 <- size_factors(sce, min.size=5, min.mean=0, ncores=1, method="igraph", plot=FALSE, verbose=FALSE)
  expect_equal(sce1, sce2)

  #negative test for min.size > max.size
  expect_error(size_factors(sce, min.size=500, max.size=5, min.mean=0, ncores=1, method="mycluster", plot=FALSE))
  # min.mean < 0
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=1 , seed=1 , min.mean= -1, plot=FALSE))
  # ncores = 0
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=0 , seed=1 , min.mean= 0, plot=FALSE))
  # seed =abc
  expect_error(size_factors(sce, min.size=500, max.size=500, ncores=1 , seed="abc" , min.mean= 1, plot=FALSE))
  # logical tests for argument
  expect_error(size_factors(sce, plot=1))
  expect_error(size_factors(sce, verbose=1))
})

test_that("truth table", {
  sc <- matrix(sample(0:200, size=2000, replace=TRUE), nrow = 200, ncol = 10)
  colnames(sc) <- c(paste0("Cell_", 1:10))
  row.names(sc) <- row.names(sce)[1:200]
  sc[1:100, 1:5] <- 7
  sc[101:200, 6:10] <- 3
  sc[1:100, 6:10] <- 17
  sc[101:200, 1:5] <- 13

  scVar <- sc
  itr <- round(runif(1) * 100)
  for(n in 1:itr){
    i <- sample(1:200, size=1, replace=TRUE) # row number
    j <- sample(1:10, size=1, replace=TRUE) #col number
    v <- sample(1:10, size=1, replace=TRUE) # value to insert
    scVar[i, j] <- v
  }
  scVar1 <- SingleCellExperiment(assays = list(counts = scVar))
  sce1 <- size_factors(scVar1, min.size=5, min.mean=0, ncores=2, plot=FALSE, verbose=FALSE)
  size.factor <- sizeFactors(sce1)
  expect_true(size.factor[1] < 1)
  expect_true(size.factor[2] < 1)
  expect_true(size.factor[3] < 1)
  expect_true(size.factor[4] < 1)
  expect_true(size.factor[5] < 1)

  expect_true(size.factor[6] > 1)
  expect_true(size.factor[7] > 1)
  expect_true(size.factor[8] > 1)
  expect_true(size.factor[9] > 1)
  expect_true(size.factor[10] > 1)
})

