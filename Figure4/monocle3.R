# Monocle3 trajectory inference of granulosa cell subtypes
# Ref: https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-case_monocle3-rstudio/tutorial.html

rm(list = ls())
options(stringsAsFactors = FALSE)

# Python environment for reading the h5ad object
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "sc_omicverse")
reticulate::use_condaenv("sc_omicverse", required = TRUE)

library(reticulate)
library(monocle3)
library(Matrix)
library(tidyverse)
library(data.table)
library(qs)

dir_pt <- './10_cytotrace2'

sc <- import('scanpy')

adata_gr <- sc$read_h5ad(file.path(dir_pt, '10_after_hb_Gr_subanndata.h5ad'))

# Raw counts (transpose: python stores genes as columns)
counts <- adata_gr$layers['counts'] %>% t()

# Metadata
meta <- adata_gr$obs
var <- adata_gr$var

# Add row/column names
colnames(counts) <- rownames(meta)
rownames(counts) <- rownames(var)

pheno.data <- meta
feature.data <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
identical(colnames(counts), rownames(pheno.data))

# Monocle3 requires a sparse matrix
data.input <- counts %>% as(., 'CsparseMatrix')

cds <- new_cell_data_set(expression_data = data.input,
                         cell_metadata = pheno.data,
                         gene_metadata = feature.data)
gc()

# Preprocess (PCA) and reduce dimension (UMAP)
cds <- cds %>%
  monocle3::preprocess_cds(cds = .,
                           method = 'PCA',
                           num_dim = 100,
                           norm_method = 'log',
                           scaling = TRUE,
                           verbose = TRUE) %>%
  monocle3::reduce_dimension(., reduction_method = 'UMAP')

colnames(colData(cds))

# Replace the UMAP with the scanpy embedding
X_umap <- adata_gr$obsm['X_umap']
rownames(X_umap) <- rownames(meta)
colnames(X_umap) <- c('UMAP-1', 'UMAP-2')
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- X_umap %>% as.matrix()
int.embed <- int.embed[rownames(cds.embed), ]
cds@int_colData$reducedDims$UMAP <- int.embed

# Cluster and partition cells (leiden)
cds <- cluster_cells(cds,
                     reduction_method = 'UMAP',
                     partition_qval = .05,
                     resolution = 0.0005)

# Learn the trajectory graph (across partitions, no loops)
cds_trajectory <- learn_graph(cds,
                              verbose = TRUE,
                              use_partition = FALSE,
                              close_loop = FALSE)

# Helper: identify the root principal point for a given cell state
get_correct_root_state <- function(cds, cell_phenotype, root_type) {
  cell_ids <- which(pData(cds)[, cell_phenotype] == root_type)

  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                             (which.max(table(closest_vertex[cell_ids, ]))))]
  root_pr_nodes
}

# FIXME: `cds_order_1` is used below but never created in this script
# (the trajectory object above is `cds_trajectory`); the root node 'Y_504' was
# also chosen manually. Confirm the intended object before rerunning.
cds_order_1 <- order_cells(cds_order_1, root_pr_nodes = 'Y_504')

# Extract pseudotime
pseudotime <- monocle3::pseudotime(cds_order_1) %>% as.data.frame()
colnames(pseudotime) <- 'monocle3_pseudotime'

# Export for the python/scanpy side
data.table::fwrite(pseudotime,
                   file = file.path(dir_pt, '10_gr_monocle3_pseudotime.csv'),
                   quote = FALSE, row.names = TRUE)

qsave(cds_order_1, file = file.path(dir_pt, '10_monocle3_cds_order.qs'))
