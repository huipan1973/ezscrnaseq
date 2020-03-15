library(BiocParallel)
library(knitr)
library(pheatmap)
library(rprojroot)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(splatter)
library(testthat)
library(reshape2)

options(stringsAsFactors = FALSE)

suppressWarnings(set.seed(seed = 100, sample.kind = "Rounding"))
#Define phase and pairs in data

sce <- readRDS(find_testthat_root_file("testing_data", "sce.rds", path=".."))
genes <- head(rownames(sce),1200)
G1 <- c(3, 5, 6, 7, 8, 10, 11, 12, 13, 15, 17, 18, 19, 21, 23, 24, 26, 27, 28, 30, 31, 32, 34, 35, 36, 37, 39, 40, 42,
43, 44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 60, 61, 63, 64, 65, 66, 69, 73, 74, 75, 77, 79, 80)
S <- c(1, 4, 9, 14, 16, 20, 25, 29, 33, 38, 41, 47, 59, 62, 67, 68, 70, 71, 72, 76, 78)
G2M <- c(2, 22)
