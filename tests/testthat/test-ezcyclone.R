context("ezcyclone")

test_that("ncore=1 vs ncore=2", {
  assignments <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  assignments2 <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)

  mypairs <- find_pairs(sce, genes=genes, G1=G1, S=S, G2M=G2M)
  assignments <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10, pairs=mypairs)
  assignments2 <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=100, min.iter=10, pairs=mypairs)
  expect_equal(assignments$phases, assignments2$phases)
})

test_that("default vs user input", {
  #ncore1 vs user  provided ncore1
  assignments <- ezcyclone(sce, organism="hsa", min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  assignments2 <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)

  #orgamism as hsa vs user  provided  organism="hsa"
  assignments <- ezcyclone(sce, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  assignments2 <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)

  #default seed vs user  provided  seed=100
  assignments <- ezcyclone(sce, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  assignments2 <- ezcyclone(sce, organism="hsa", seed=100, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)

  #default iter vs user  provided iter = 1000,
  assignments <- ezcyclone(sce, min.pairs=5, verbose=FALSE, min.iter=10)
  assignments2 <- ezcyclone(sce, organism="hsa", seed=100, min.pairs=5, verbose=FALSE, iter=1000, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)

  #default min.iter vs user  provided min.iter = 100,
  assignments <- ezcyclone(sce, min.pairs=5, verbose=FALSE)
  assignments2 <- ezcyclone(sce, organism="hsa", seed=100, min.pairs=5, verbose=FALSE, iter=1000, min.iter=100)
  expect_equal(assignments$phases, assignments2$phases)

  #default  min.pairs =5  vs user  provided min.pairs = 5 with diff ncore ,
  assignments <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5 , verbose=FALSE)
  assignments2 <- ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE)
  expect_equal(assignments$phases, assignments2$phases)

  #default  gene.names vs user  provided gene.names
  assignments <- ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5 , verbose=FALSE)
  assignments2 <- ezcyclone(sce, organism="hsa", gene.names=row.names(sce),  ncores=2, min.pairs=5, verbose=FALSE)
  expect_equal(assignments$phases, assignments2$phases)
})

test_that("mmu test", {
  #genes in pairs not present in sce
  expect_error(assignments <- ezcyclone(sce, organism="mmu", min.pairs=5, verbose=FALSE, iter=10, min.iter=3))

  #converting hsa sce to mmu sce by replacing with mmu genes
  mmu.sce <- sce
  pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
  row.names(mmu.sce) <- head(pairs$S[,1], n=2000)
  assignments <- ezcyclone(mmu.sce, organism="mmu", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  assignments2 <- ezcyclone(mmu.sce, organism="mmu", ncores=2, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments$phases, assignments2$phases)
})

test_that("negative tests",{
  emptypairs <- scran::sandbag(head(sce), list(G1=G1, S=S, G2M=G2M))
  expect_error(ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=100, min.iter=10, pairs=emptypairs))

  #pairs which return NA
  p2 <- list()
  p2$S <- as.data.frame(t(genes[3:4]))
  colnames(p2$S) <- c("first","second")
  p2$G1 <- as.data.frame(t(genes[6:7]))
  colnames(p2$G1) <- c("first","second")
  p2$G2M <- as.data.frame(t(genes[1:2]))
  colnames(p2$G2M) <-c ("first","second")
  expect_error(ezcyclone(sce, organism="hsa", ncores=2, min.pairs=5, verbose=FALSE, iter=10, min.iter=3, pairs=p2))

  # nrow(sce) != length(genes)
  expect_error(ezcyclone(sce, organism="hsa", gene.names=genes, ncores=2, min.pairs=5, verbose=FALSE, iter=10))

  #gene.names not in rownames (sce)
  genes2 <- row.names(sce)
  genes2[1] <- "MyFavGene"
  expect_error(ezcyclone(sce, organism="hsa", gene.names=genes2, ncores=2, min.pairs=5, verbose=FALSE, iter=10))

  # organism other then hsa and mmu
  expect_error(assignments <- ezcyclone(sce, organism="myOrganism", min.pairs=5, verbose=FALSE, iter=100, min.iter=10))
  #min.pairs is negative
  expect_error(ezcyclone(sce, organism="hsa", ncores=1, min.pairs= -1, verbose=FALSE, iter=102, min.iter=100 ))
  expect_error(ezcyclone(sce, organism="hsa", ncores=1, min.pairs= 0, verbose=FALSE, iter=102, min.iter=100 ))
  #iter < min.iter
  expect_error(ezcyclone(sce, organism="hsa", ncores=1, min.pairs= -1, verbose=FALSE, iter=10, min.iter=100 ))
  #ncore=0
  expect_error(ezcyclone(sce, organism="hsa", ncores=0, min.pairs= -1, verbose=FALSE, iter=10, min.iter=100 ))
  #verbose not logical
  expect_error(ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=1, iter=100, min.iter=10 ))
  #seed not numeric
  expect_error(ezcyclone(sce, organism="hsa", ncores=1, min.pairs=5, verbose=1, iter=100, min.iter=10, seed="abc"))
})

test_that("truth table", {
  sc <- matrix(sample(0:200, 12000, replace=TRUE), nrow = 1200, ncol = 10)
  colnames(sc) <- c(paste0("Cell_", 1:10))
  row.names(sc) <- row.names(sce)[1:1200]
  scVar1 <- SingleCellExperiment(assays = list(counts = sc))
  G1 <- 1:3
  S <- 4:8
  G2M <- 9:10
  pairs <- scran::sandbag(scVar1, list(G1=G1, S=S, G2M=G2M ))

  assignments <- ezcyclone(scVar1, organism="hsa", ncores=1, min.pairs=5, verbose=FALSE, iter=100, min.iter=10, pairs=pairs)
  assignments.2 <- cyclone(scVar1, pairs=pairs, min.pairs=5, verbose=FALSE, iter=100, min.iter=10)
  expect_equal(assignments, assignments.2)
})
