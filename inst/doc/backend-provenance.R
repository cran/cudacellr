## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----install, eval=FALSE------------------------------------------------------
# # install.packages("pak")
# pak::pak("cudaverse/cudacellr")
# library(cudacellr)
# 
# health <- cudaverse::cuda_diagnostics()
# health$summary
# health$next_steps
# health$selected_backend  # "native" for the lightweight CUDA path
# cudaverse::cuda_select_device("cuda")

## ----counts-------------------------------------------------------------------
library(cudacellr)
library(Matrix)

set.seed(42)
counts <- matrix(
  rpois(60 * 30, lambda = 2),
  nrow = 60,
  ncol = 30,
  dimnames = list(
    paste0("gene_", seq_len(60)),
    paste0("cell_", seq_len(30))
  )
)
counts[1, ] <- counts[1, ] + 1
counts <- Matrix(counts, sparse = TRUE)

## ----cuda-workflow, eval=FALSE------------------------------------------------
# fit <- cudacell_workflow(
#   counts,
#   n_hvg = 20,
#   n_components = 5,
#   k = 5,
#   batch_size = 8,
#   device = "cuda"
# )
# 
# fit
# head(fit$pca$x)
# head(fit$neighbors$index)
# head(fit$neighbors$distance)
# 
# identical(rownames(fit$pca$x), colnames(counts))
# identical(rownames(fit$neighbors$index), colnames(counts))

## ----cuda-stages, eval=FALSE--------------------------------------------------
# pca <- cuda_cell_pca(
#   counts,
#   n_hvg = 20,
#   n_components = 5,
#   device = "cuda"
# )
# neighbors <- cuda_cell_neighbors(
#   pca,
#   k = 5,
#   batch_size = 8,
#   device = "cuda"
# )

## ----cuda-provenance, eval=FALSE----------------------------------------------
# workflow_provenance <- cuda_provenance(fit)
# workflow_provenance
# attr(workflow_provenance, "compute_device")

## ----cpu-reference------------------------------------------------------------
reference <- cudacell_workflow(
  counts,
  n_hvg = 20,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

reference
identical(rownames(reference$pca$x), colnames(counts))
identical(rownames(reference$neighbors$index), colnames(counts))
cuda_provenance(reference)

