## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
# BiocManager::install("SingleCellExperiment")

## ----fixture------------------------------------------------------------------
library(cudacellr)
library(Matrix)
library(SingleCellExperiment)

set.seed(42)
counts <- matrix(
  rpois(80 * 30, lambda = 2),
  nrow = 80,
  ncol = 30,
  dimnames = list(
    paste0("gene_", seq_len(80)),
    paste0("cell_", seq_len(30))
  )
)
counts[1, ] <- counts[1, ] + 1
counts <- Matrix(counts, sparse = TRUE)

sce <- SingleCellExperiment(
  assays = list(
    counts = counts,
    untouched = counts * 2
  ),
  rowData = S4Vectors::DataFrame(
    symbol = paste0("symbol_", seq_len(nrow(counts)))
  ),
  colData = S4Vectors::DataFrame(
    batch = rep(c("one", "two"), each = 15)
  ),
  metadata = list(owner = "unchanged"),
  reducedDims = list(
    EXISTING = matrix(rnorm(60), nrow = 30)
  )
)
sizeFactors(sce) <- rep(1, ncol(sce))
altExp(sce, "spike") <- SingleCellExperiment(
  assays = list(
    counts = matrix(
      rpois(5 * ncol(sce), 1),
      nrow = 5,
      dimnames = list(paste0("spike_", 1:5), colnames(sce))
    )
  )
)

## ----workflow-----------------------------------------------------------------
result <- cudacell_sce(
  sce,
  assay = "counts",
  n_hvg = 30,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

assayNames(result)
reducedDimNames(result)
colPairNames(result)

## ----outputs------------------------------------------------------------------
assay(result, "cudacell_logcounts")[1:3, 1:3]
reducedDim(result, "CUDACELL_PCA")[1:3, , drop = FALSE]
rowData(result)[
  1:3,
  c(
    "cudacell_mean",
    "cudacell_variance",
    "cudacell_dispersion",
    "cudacell_hvg",
    "cudacell_hvg_rank"
  )
]
colPair(result, "CUDACELL_KNN")

## ----preservation-------------------------------------------------------------
identical(assay(sce, "counts"), assay(result, "counts"))
identical(assay(sce, "untouched"), assay(result, "untouched"))
identical(rowData(sce)$symbol, rowData(result)$symbol)
identical(colData(sce)$batch, colData(result)$batch)
identical(reducedDim(sce, "EXISTING"), reducedDim(result, "EXISTING"))
identical(altExp(sce, "spike"), altExp(result, "spike"))
identical(sizeFactors(sce), sizeFactors(result))
identical(metadata(sce)$owner, metadata(result)$owner)

## ----normalization, eval=FALSE------------------------------------------------
# result <- cudacell_sce(
#   sce,
#   log1p = FALSE,
#   normalized_assay = "normalized_counts",
#   scale_factor = 10000,
#   n_hvg = 30,
#   n_components = 5,
#   k = 5,
#   device = "cpu"
# )

## ----collision, error=TRUE----------------------------------------------------
try({
cudacell_sce(
  result,
  n_hvg = 30,
  n_components = 5,
  k = 5,
  device = "cpu"
)
})

## ----provenance---------------------------------------------------------------
cuda_provenance(result)
metadata(result)$cudacellr$parameters

