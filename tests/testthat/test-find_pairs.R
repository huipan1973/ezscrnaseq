context("find_pairs")

test_that("tests", {
  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M)
  genes2 <- row.names(sce)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes2)
  expect_equal(mypairs, mypairs2)

  mypairs <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes)
  genes2 <- rbind(genes,genes)
  mypairs2 <- find_pairs(sce, G1=G1, S=S, G2M=G2M, genes=genes)
  expect_equal(mypairs, mypairs2)

  #test for nrow(sce) <1000
  sce2 <- head(sce, 999)
  expect_error(find_pairs(sce2, G1=G1, S=S, G2M=G2M))

  #test for empty pairs
  sce3 <- rbind(rep(head(sce),times=200))
  expect_error(find_pairs(sce3, G1=G1, S=S, G2M=G2M))
})

test_that("truth table", {
  sc <- matrix(sample(0:200, 12000, replace=T), nrow = 1200, ncol = 10)
  colnames(sc) <- c(paste0("Cell_", 1:10))
  row.names(sc) <- row.names(sce)[1:1200]
  scVar1 <- SingleCellExperiment(assays = list(counts = sc))
  G1 <- 1:3
  S <- 4:8
  G2M <- 9:10
  pairs <- scran::sandbag(scVar1, list(G1=G1, S=S, G2M=G2M ))
  mypairs <- find_pairs(scVar1, G1=G1, S=S, G2M=G2M)
  expect_equal(pairs, mypairs)
})
