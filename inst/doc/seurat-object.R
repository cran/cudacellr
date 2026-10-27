## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("SeuratObject", quietly = TRUE)
)

## ----install, eval=FALSE------------------------------------------------------
# install.packages("SeuratObject")

## ----fixture------------------------------------------------------------------
library(cudacellr)
library(Matrix)

set.seed(42)
counts <- Matrix(
  matrix(
    rpois(80 * 30, lambda = 2),
    nrow = 80,
    dimnames = list(
      paste0("gene_", seq_len(80)),
      paste0("cell_", seq_len(30))
    )
  ),
  sparse = TRUE
)
counts[1, ] <- counts[1, ] + 1

object <- SeuratObject::CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  meta.data = data.frame(
    batch = rep(c("one", "two"), each = 15),
    row.names = colnames(counts)
  ),
  project = "preserved-project"
)
object[["ADT"]] <- SeuratObject::CreateAssay5Object(
  counts = counts[seq_len(5), , drop = FALSE]
)
object[["existing_pca"]] <- SeuratObject::CreateDimReducObject(
  embeddings = matrix(
    rnorm(ncol(counts) * 2),
    nrow = ncol(counts),
    dimnames = list(colnames(counts), c("OLDPC_1", "OLDPC_2"))
  ),
  assay = "RNA",
  key = "OLDPC_"
)
SeuratObject::Idents(object) <- factor(
  rep(c("A", "B"), each = 15)
)
SeuratObject::Misc(object, slot = "owner") <- "unchanged"
add_example_tool <- function(x) {
  SeuratObject::Tool(x) <- list(owner = "unchanged")
  x
}
object <- add_example_tool(object)
existing_tool_names <- SeuratObject::Tool(object)

## ----workflow-----------------------------------------------------------------
result <- cudacell_seurat(
  object,
  assay = "RNA",
  layer = "counts",
  n_hvg = 30,
  n_components = 5,
  k = 5,
  batch_size = 8,
  device = "cpu"
)

SeuratObject::Assays(result)
SeuratObject::Reductions(result)
SeuratObject::Neighbors(result)

## ----outputs------------------------------------------------------------------
normalized <- SeuratObject::LayerData(
  result[["CUDACELL"]],
  layer = "data",
  fast = FALSE
)
normalized[1:3, 1:3]

SeuratObject::Embeddings(result[["cudacell_pca"]])[1:3, , drop = FALSE]
SeuratObject::Indices(result[["cudacell_knn"]])[1:3, , drop = FALSE]
result[["CUDACELL"]][[]][
  1:3,
  c(
    "cudacell_mean",
    "cudacell_variance",
    "cudacell_dispersion",
    "cudacell_hvg",
    "cudacell_hvg_rank"
  )
]

## ----preservation-------------------------------------------------------------
identical(
  SeuratObject::LayerData(object[["RNA"]], "counts"),
  SeuratObject::LayerData(result[["RNA"]], "counts")
)
identical(object[["ADT"]], result[["ADT"]])
identical(object[["existing_pca"]], result[["existing_pca"]])
identical(SeuratObject::Idents(object), SeuratObject::Idents(result))
identical(
  SeuratObject::Misc(object, slot = "owner"),
  SeuratObject::Misc(result, slot = "owner")
)
all(vapply(
  existing_tool_names,
  function(name) {
    identical(
      SeuratObject::Tool(object, slot = name),
      SeuratObject::Tool(result, slot = name)
    )
  },
  logical(1)
))

## ----collision, error=TRUE----------------------------------------------------
try({
cudacell_seurat(
  result,
  n_hvg = 30,
  n_components = 5,
  k = 5,
  device = "cpu"
)
})

## ----provenance---------------------------------------------------------------
cuda_provenance(result)
SeuratObject::Tool(result, slot = "cudacell_seurat")$outputs
SeuratObject::Tool(result, slot = "cudacell_seurat")$parameters

