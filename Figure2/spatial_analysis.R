# Spatial transcriptomics (Slide-seq, GSE240271 ovary): RCTD deconvolution of
# the 0h and Immature (untreated) samples with the single-cell reference.

rm(list = ls())
options(stringsAsFactors = FALSE)

# Python environment for reading h5ad files
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "sc_ss_ov")
reticulate::use_condaenv("sc_ss_ov", required = TRUE)

library(qs)
library(tidyverse)
library(data.table)
library(Seurat)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(future)
library(spacexr)
require(Matrix)

dir_st <- './16_spatial/GSE240271'
dir_figures <- './figures'

sc <- reticulate::import('scanpy')

# ---- Sample 1 (0h): build the spatial object and run RCTD ---------------------

all_meta <- fread(file.path(dir_st, 'GSE240271_meta_all.csv'), data.table = FALSE)

meta_0h <- all_meta %>% filter(Sample == '0h')
rownames(meta_0h) <- meta_0h$barcode

adata_dir <- file.path(dir_st, 'GSE240271_untar/GSM7689280_adata_ovary_0hr_spatial_raw_counts.h5ad')

h5ad2seurat <- function(adata_dir) {
  adata <- sc$read_h5ad(adata_dir)  # annotated data

  # metadata for each spot
  meta <- adata$obs
  # metadata for each gene
  var <- adata$var
  # expression matrix; transpose (python stores genes as columns)
  counts_data <- adata$X %>% Matrix::t()

  colnames(counts_data) <- rownames(meta)
  rownames(counts_data) <- rownames(var)

  counts_data <- counts_data[, rownames(meta_0h)]
  seurat_data <- CreateSeuratObject(counts = counts_data,
                                    meta.data = meta_0h,
                                    project = 'Ovary_ST',
                                    assay = 'Spatial')
}

seurat_data <- h5ad2seurat(adata_dir)
adata <- sc$read_h5ad(adata_dir)
total_count <- as.matrix(seurat_data@assays[["Spatial"]]@counts)

coords <- adata$obsm['spatial'] %>% as.data.frame()
names(coords) <- c("xcoord", "ycoord")
rownames(coords) <- adata$obs %>% rownames(.)
coords <- coords[seurat_data@meta.data %>% rownames(.), ]

nUMI <- colSums(total_count)
puck <- SpatialRNA(coords, total_count, nUMI)

# ---- Single-cell reference (shared by both samples) ---------------------------

seurat_dat <- qread('./17_cell_communication/seurat_dat_for_cellchat.qs')

# Granulosa subtype annotation
Gr_subannotation <- data.table::fread('./10_cytotrace2/pseudotime_res_all.csv', data.table = FALSE)
rownames(Gr_subannotation) <- Gr_subannotation[, 1]
Gr_subannotation <- Gr_subannotation[, -1]
Gr_subannotation <- Gr_subannotation %>% dplyr::select(Gr_subannotation)

meta_all <- seurat_dat@meta.data[, 'Major_annotation_d', drop = FALSE]

meta_use <- meta_all %>% merge(., Gr_subannotation, by = 0, all.x = TRUE, sort = FALSE)
# New column Cluster_all: start from major annotation
meta_use[['Cluster_all']] <- meta_use[['Major_annotation_d']] %>% as.character()
# Inject granulosa subtype labels into the major annotation
meta_use[['Cluster_all']][meta_use[['Cluster_all']] == 'Granulosa'] <-
  meta_use$Gr_subannotation[meta_use[['Cluster_all']] == 'Granulosa']
# Cells without subtype info (e.g. red blood cells) are named GGGG
meta_use[['Cluster_all']][is.na(meta_use[['Cluster_all']])] <- 'GGGG'

table(meta_use$Cluster_all)

rownames(meta_use) <- meta_use$Row.names
identical(meta_use$Row.names, seurat_dat@meta.data %>% rownames(.))
# merge() reorders rows; restore the original cell order
meta_use <- meta_use[seurat_dat@meta.data %>% rownames(.), ]
identical(meta_use$Row.names, seurat_dat@meta.data %>% rownames(.))

seurat_dat@meta.data$Cluster_all <- meta_use$Cluster_all
# Remove GGGG cells
seurat_dat <- seurat_dat %>% subset(Cluster_all != 'GGGG')
seurat_dat$Cluster_all <- factor(seurat_dat$Cluster_all)

sc_counts <- seurat_dat@assays$RNA@counts

cell_types <- seurat_dat$Cluster_all
names(cell_types) <- rownames(seurat_dat@meta.data)
nUMI <- seurat_dat$nCount_RNA
names(nUMI) <- rownames(seurat_dat@meta.data)

reference <- Reference(sc_counts, cell_types, nUMI, n_max_cells = 1000000)

myRCTD <- create.RCTD(puck, reference, max_cores = 12) %>%
  run.RCTD(., doublet_mode = 'doublet')

qsave(myRCTD, file.path(dir_st, 'GSE240271_0h_ST_myRCTD.qs'))

# ---- Sample 2 (Immature/untreated): RCTD --------------------------------------

untreated_meta <- fread('./16_spatial/GSE240271_untreated_metadata.csv', data.table = FALSE)
rownames(untreated_meta) <- untreated_meta$barcode

all_meta <- fread(file.path(dir_st, 'GSE240271_meta_all.csv'), data.table = FALSE)
table(all_meta$Celltypes)

adata_dir <- file.path(dir_st, 'GSE240271_untar/GSM7689279_adata_ovary_Immature_spatial_raw_counts.h5ad')

h5ad2seurat <- function(adata_dir) {
  adata <- sc$read_h5ad(adata_dir)  # annotated data

  # metadata for each spot
  meta <- adata$obs
  # metadata for each gene
  var <- adata$var
  # expression matrix; transpose (python stores genes as columns)
  counts_data <- adata$X %>% Matrix::t()

  colnames(counts_data) <- rownames(meta)
  rownames(counts_data) <- rownames(var)

  counts_data <- counts_data[, rownames(untreated_meta)]
  seurat_data <- CreateSeuratObject(counts = counts_data,
                                    meta.data = untreated_meta,
                                    project = 'Ovary_ST',
                                    assay = 'Spatial')
}

seurat_data <- h5ad2seurat(adata_dir)
adata <- sc$read_h5ad(adata_dir)
total_count <- as.matrix(seurat_data@assays[["Spatial"]]@counts)

coords <- adata$obsm['spatial'] %>% as.data.frame()
names(coords) <- c("xcoord", "ycoord")
rownames(coords) <- adata$obs %>% rownames(.)
coords <- coords[seurat_data@meta.data %>% rownames(.), ]

nUMI <- colSums(total_count)
puck <- SpatialRNA(coords, total_count, nUMI)

# Single-cell reference (same as above)
seurat_dat <- qread('./17_cell_communication/seurat_dat_for_cellchat.qs')
Gr_subannotation <- data.table::fread('./10_cytotrace2/pseudotime_res_all.csv', data.table = FALSE)
rownames(Gr_subannotation) <- Gr_subannotation[, 1]
Gr_subannotation <- Gr_subannotation[, -1]
Gr_subannotation <- Gr_subannotation %>% dplyr::select(Gr_subannotation)

meta_all <- seurat_dat@meta.data[, 'Major_annotation_d', drop = FALSE]
meta_use <- meta_all %>% merge(., Gr_subannotation, by = 0, all.x = TRUE, sort = FALSE)
meta_use[['Cluster_all']] <- meta_use[['Major_annotation_d']] %>% as.character()
meta_use[['Cluster_all']][meta_use[['Cluster_all']] == 'Granulosa'] <-
  meta_use$Gr_subannotation[meta_use[['Cluster_all']] == 'Granulosa']
meta_use[['Cluster_all']][is.na(meta_use[['Cluster_all']])] <- 'GGGG'

rownames(meta_use) <- meta_use$Row.names
# merge() reorders rows; restore the original cell order
meta_use <- meta_use[seurat_dat@meta.data %>% rownames(.), ]

seurat_dat@meta.data$Cluster_all <- meta_use$Cluster_all
seurat_dat <- seurat_dat %>% subset(Cluster_all != 'GGGG')
seurat_dat$Cluster_all <- factor(seurat_dat$Cluster_all)

sc_counts <- seurat_dat@assays$RNA@counts

cell_types <- seurat_dat$Cluster_all
names(cell_types) <- rownames(seurat_dat@meta.data)
nUMI <- seurat_dat$nCount_RNA
names(nUMI) <- rownames(seurat_dat@meta.data)

reference <- Reference(sc_counts, cell_types, nUMI, n_max_cells = 1000000)

myRCTD <- create.RCTD(puck, reference, max_cores = 8) %>%
  run.RCTD(., doublet_mode = 'doublet')

qsave(myRCTD, file.path(dir_st, 'GSE240271_untreat_ST_myRCTD.qs'))

# Attach normalized RCTD weights to the untreated spatial object
myRCTD <- qread(file.path(dir_st, 'GSE240271_untreat_ST_myRCTD.qs'))

norm_weights <- normalize_weights(myRCTD@results[["weights"]]) %>%
  as.matrix() %>% as.data.frame()
identical(rownames(norm_weights), seurat_data@meta.data %>% rownames(.))

seurat_data@meta.data <- seurat_data@meta.data %>% cbind(., norm_weights)
qsave(seurat_data, file = file.path(dir_st, 'GSE240271_untreat_ST_sce.qs'))

seurat_data <- qread(file.path(dir_st, 'GSE240271_untreat_ST_sce.qs'))

# Spot annotation table (untreated)
df <- table(seurat_data$Celltypes, seurat_data$BroadCelltype) %>%
  as.data.frame() %>% filter(Freq != 0) %>% t() %>% as.data.frame()
data.table::fwrite(x = df, './16_spatial/untreat_Table_S1.csv', quote = FALSE, row.names = TRUE)

seurat_data@images$image <- new(Class = 'SlideSeq',
                                assay = 'Spatial',
                                key = 'image_',
                                coordinates = seurat_data@meta.data[, c('x', 'y')])

# ---- Sample 1 (0h): attach RCTD weights ---------------------------------------

all_meta <- fread(file.path(dir_st, 'GSE240271_meta_all.csv'), data.table = FALSE)
meta_0h <- all_meta %>% filter(Sample == '0h')
rownames(meta_0h) <- meta_0h$barcode

adata_dir <- file.path(dir_st, 'GSE240271_untar/GSM7689280_adata_ovary_0hr_spatial_raw_counts.h5ad')

h5ad2seurat <- function(adata_dir) {
  adata <- sc$read_h5ad(adata_dir)  # annotated data

  # metadata for each spot
  meta <- adata$obs
  # metadata for each gene
  var <- adata$var
  # expression matrix; transpose (python stores genes as columns)
  counts_data <- adata$X %>% Matrix::t()

  colnames(counts_data) <- rownames(meta)
  rownames(counts_data) <- rownames(var)

  counts_data <- counts_data[, rownames(meta_0h)]
  seurat_data <- CreateSeuratObject(counts = counts_data,
                                    meta.data = meta_0h,
                                    project = 'Ovary_ST',
                                    assay = 'Spatial')
}

seurat_data <- h5ad2seurat(adata_dir)

myRCTD <- qread(file.path(dir_st, 'GSE240271_0h_ST_myRCTD.qs'))

norm_weights <- normalize_weights(myRCTD@results[["weights"]]) %>%
  as.matrix() %>% as.data.frame()
identical(rownames(norm_weights), seurat_data@meta.data %>% rownames(.))

seurat_data@meta.data <- seurat_data@meta.data %>% cbind(., norm_weights)
qsave(seurat_data, file = file.path(dir_st, 'GSE240271_0h_ST_sce.qs'))

seurat_st1 <- qread(file.path(dir_st, 'GSE240271_0h_ST_sce.qs'))

# Spot annotation table (0h)
df <- table(seurat_st1$Celltypes, seurat_st1$BroadCelltype) %>%
  as.data.frame() %>% filter(Freq != 0) %>% t() %>% as.data.frame()
data.table::fwrite(x = df, './16_spatial/0h_Table_S1.csv', quote = FALSE, row.names = TRUE)

seurat_st1@images$image <- new(Class = 'SlideSeq',
                               assay = 'Spatial',
                               key = 'image_',
                               coordinates = seurat_st1@meta.data[, c('x', 'y')])

Idents(seurat_st1) <- seurat_st1$Celltypes
seurat_st2 <- qread(file.path(dir_st, 'GSE240271_untreat_ST_sce.qs'))
seurat_st2@images$image <- new(Class = 'SlideSeq',
                               assay = 'Spatial',
                               key = 'image_',
                               coordinates = seurat_st2@meta.data[, c('x', 'y')])
Idents(seurat_st2) <- seurat_st2$Celltypes

seurat_st_ls <- list(seurat_st1 = seurat_st1,
                     seurat_st2 = seurat_st2)

# ---- Broad cell type spatial plots --------------------------------------------

col_use <- c('#ff7f0c', '#279e68', '#d52425', '#ffff33', '#aec7e8', '#b5bd61')

color_Broad_cells <- c('#1e77b4', '#FFA500', '#008000', '#FF0000', '#800080',
                       '#8B4513', '#FFC0CB', '#b5bd60', '#00CED1', '#B0C4DE',
                       '#FFD700')

for (i in seq_along(seurat_st_ls)) {
  st_seurat <- seurat_st_ls[[i]]
  Idents(st_seurat) <- st_seurat$BroadCelltype

  p_spatial_dim <- SpatialDimPlot(st_seurat, group.by = 'BroadCelltype',
                                  crop = TRUE, combine = TRUE,
                                  stroke = 0.0000001) +
    ggtitle(label = st_seurat$Sample[1]) +
    scale_color_manual(values = col_use) +
    scale_fill_manual(values = col_use) +
    guides(fill = guide_legend(override.aes = list(size = 2, color = 'black'))) +
    theme(legend.position = c(0.1, 0.15),
          legend.background = element_blank(),
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          legend.title = element_blank(),
          plot.title = element_text(size = 16, hjust = .5))
  p_spatial_dim <- p_spatial_dim[[1]]

  # Sum of the six GC subtype weights as overall granulosa proportion
  GC_indexs <- paste0('GC_', 0:5)
  st_seurat$Granulosa <- rowSums(st_seurat@meta.data[, GC_indexs])

  index_features <- c("B", "Endothelial", 'Epithelial', 'Granulosa', 'Macrophage',
                      'Neutrophil', 'NK T cell', 'Oocyte', 'Perivascular',
                      'Stromal', 'Theca')
  p_feature <- SpatialFeaturePlot(st_seurat, features = index_features,
                                  ncol = 6, pt.size.factor = 1.8,
                                  combine = FALSE, stroke = 0)
  p_feature <- map(seq_along(p_feature), function(x) {
    a <- p_feature[[x]] +
      scale_fill_gradientn(colors = colorRampPalette(c('black', color_Broad_cells[x]))(100)) +
      ggtitle(label = paste0('RCTD: ', index_features[x])) +
      theme(panel.background = element_rect(fill = 'black'),
            panel.grid = element_blank(),
            legend.key.size = unit(.25, 'cm'),
            plot.title = element_text(size = 16, hjust = .5),
            legend.position = c(0.1, 0.15),
            legend.title = element_blank())
  })

  p_ls <- list()
  p_ls[[1]] <- p_spatial_dim
  p_ls[2:12] <- p_feature

  p_ls_comb <- p_ls %>% patchwork::wrap_plots(., ncol = 6)

  save_dir_pdf <- file.path(dir_figures, paste0('18_spatial_', st_seurat$Sample[1], '.pdf'))
  save_dir_png <- file.path(dir_figures, paste0('18_spatial_', st_seurat$Sample[1], '.png'))

  ggsave(save_dir_pdf, p_ls_comb, width = 3 * 6, height = 3 * 2, limitsize = FALSE)
  ggsave(save_dir_png, p_ls_comb, width = 3 * 6, height = 3 * 2, dpi = 600, limitsize = FALSE)
}

# ---- Granulosa subtype spatial plots ------------------------------------------

color_GC_cells <- c('#1e77b4', '#FFA500', '#008000', '#FF0000', '#800080',
                    '#8B4513')

p_all <- list()

for (i in seq_along(seurat_st_ls)) {
  st_seurat <- seurat_st_ls[[i]]
  Idents(st_seurat) <- st_seurat$Celltypes

  # Fine-grained spot annotation
  p_spatial_dim <- SpatialDimPlot(st_seurat, group.by = 'Celltypes',
                                  crop = TRUE, combine = TRUE,
                                  stroke = 0.0000001) +
    ggtitle(label = st_seurat$Sample[1]) +
    ggsci::scale_fill_d3('category20') +
    guides(fill = guide_legend(override.aes = list(size = 2, color = 'black'))) +
    theme(legend.position = c(0.1, 0.15),
          legend.background = element_blank(),
          legend.key.height = unit(0.2, 'cm'),
          legend.key.width = unit(.2, 'cm'),
          legend.title = element_blank(),
          plot.title = element_text(size = 16, hjust = .5))
  p_spatial_dim <- p_spatial_dim[[1]]

  GC_st_indexs <- grep('GC', levels(Idents(st_seurat)), value = TRUE)

  p_spatial_GC <- SpatialDimPlot(st_seurat,
                                 cells.highlight = CellsByIdentities(object = st_seurat,
                                                                     idents = GC_st_indexs),
                                 facet.highlight = TRUE,
                                 cols.highlight = c('red', '#d4d4d4'),
                                 stroke = 0.000001,
                                 crop = TRUE)
  p_spatial_GC <- map(seq_along(p_spatial_GC), function(x) {
    a <- p_spatial_GC[[x]] +
      theme(panel.background = element_rect(fill = 'black'),
            panel.grid = element_blank(),
            plot.title = element_text(size = 16, hjust = .5),
            legend.title = element_blank())
  })

  GC_indexs <- paste0('GC_', 0:5)
  p_feature <- SpatialFeaturePlot(st_seurat, features = GC_indexs,
                                  ncol = 6, pt.size.factor = 1.4,
                                  combine = FALSE, stroke = 0)
  p_feature <- map(seq_along(p_feature), function(x) {
    a <- p_feature[[x]] +
      scale_fill_gradientn(colors = colorRampPalette(c('white', color_GC_cells[x]))(100)) +
      ggtitle(label = paste0('RCTD: ', GC_indexs[x])) +
      theme(panel.background = element_rect(fill = 'black'),
            panel.grid = element_blank(),
            legend.key.size = unit(.25, 'cm'),
            plot.title = element_text(size = 16, hjust = .5),
            legend.position = c(0.1, 0.15),
            legend.title = element_blank())
  })

  p_all[[st_seurat$Sample[1]]] <- list(p_spatial_dim = p_spatial_dim,
                                       p_spatial_GC = p_spatial_GC,
                                       p_feature = p_feature)
}

p_save <- p_all[[1]][[1]]
ggview:::ggview(p_save, width = 7, height = 7)

p_save <- p_all[[1]][[1]]
ggsave(file.path(dir_figures, '18_spatial_0h_celltype_all.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

p_save <- p_all[[2]][[1]]
ggsave(file.path(dir_figures, '18_spatial_untreat_celltype_all.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

p_save <- p_all[[1]][[2]] %>% patchwork::wrap_plots(., ncol = 3)
ggsave(file.path(dir_figures, '18_spatial_0h_GC_each.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

p_save <- p_all[[2]][[2]] %>% patchwork::wrap_plots(., ncol = 3)
ggsave(file.path(dir_figures, '18_spatial_untreat_GC_each.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

p_save <- p_all[[1]][[3]] %>% patchwork::wrap_plots(., ncol = 3)
ggsave(file.path(dir_figures, '18_spatial_0h_GC_RCTD.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

p_save <- p_all[[2]][[3]] %>% patchwork::wrap_plots(., ncol = 3)
ggsave(file.path(dir_figures, '18_spatial_untreat_GC_RCTD.png'), p_save,
       width = 7, height = 7, limitsize = FALSE, dpi = 600)

# GC subtype proportion per spot (single sample preview)
Idents(seurat_st1) <- seurat_st1$Celltypes
GC_st_indexs <- grep('GC', levels(Idents(seurat_st1)), value = TRUE)
VlnPlot(seurat_st1, features = paste0('GC_', 0:5),
        idents = GC_st_indexs, pt.size = 0, stack = TRUE)

# ---- Merge both samples: stacked GC proportion violin -------------------------

dat_st1 <- seurat_st1@meta.data
dat_st2 <- seurat_st2@meta.data
dat_st2 <- dat_st2[, -c(70:75)]

seurat_st_for_merge <- list()
for (i in seq_along(seurat_st_ls)) {
  count_data <- GetAssayData(seurat_st_ls[[i]], slot = 'counts')
  metadata <- seurat_st_ls[[i]]@meta.data
  rownames(metadata) <- metadata[['V1']]
  colnames(count_data) <- rownames(metadata)

  sample_name_index <- metadata[['Sample']][1]
  seurat_data <- CreateSeuratObject(counts = count_data,
                                    meta.data = metadata,
                                    project = sample_name_index,
                                    assay = 'Spatial')
  seurat_data@images$image <- new(Class = 'SlideSeq',
                                  assay = 'Spatial',
                                  key = 'image_',
                                  coordinates = seurat_data@meta.data[, c('x', 'y')])
  seurat_st_for_merge[[sample_name_index]] <- seurat_data
}

merge_st_dat <- seurat_st_for_merge[[1]] %>% merge(., seurat_st_for_merge[-1])
Idents(merge_st_dat) <- merge_st_dat$Celltypes

p_vlnplot <- VlnPlot(merge_st_dat, features = paste0('GC_', 0:5),
                     idents = GC_st_indexs, pt.size = 0, stack = TRUE,
                     cols = color_GC_cells) +
  labs(x = 'RCTD_cell_proportion', y = 'ST_cell_type') +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none',
        strip.text.x = element_text(angle = 0))

library(aplot)
p.rect <- ggplot() +
  geom_rect(xmin = 0, xmax = 25.5, ymin = 0, ymax = 13.5, fill = 'transparent') +
  theme_bw(base_rect_size = 1.5) +
  theme(plot.background = element_blank(),
        panel.background = element_blank())
p.rect

p.use <- p_vlnplot + patchwork::inset_element(p.rect, left = -0.02, bottom = -0.05,
                                              right = 1.02, top = 1.05, align_to = 'plot')
ggview:::ggview(p.use, width = 5 * 2.2, height = 3.5, bg = 'white')

ggsave(file.path(dir_figures, '18_spatial_both_GC_vlnplot.png'), p.use,
       width = 5 * 2.2, height = 3.5, bg = 'white')
ggsave(file.path(dir_figures, '18_spatial_both_GC_vlnplot.pdf'), p.use,
       width = 5 * 2.2, height = 3.5, bg = 'white')
