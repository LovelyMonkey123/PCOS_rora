# External cohort validation (GSE240688): RORA(+) cell signature AUCell score

rm(list = ls())
options(stringsAsFactors = FALSE)

# Python environment for reading the h5ad object
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "celloracle_env")
reticulate::use_condaenv("celloracle_env", required = TRUE)

pkgs <- c("ggplot2", "dplyr", "tidyr", "tibble", "reshape2", "Seurat", "reticulate")
sapply(pkgs, require, character.only = TRUE)
sc <- import('scanpy')

library(qs)
library(patchwork)
library(ggsci)

dir_val <- './19_validation'
dir_figures <- './figures'

source('./source/analysis_cellchat_ssh.R')

# Convert the external cohort h5ad to a Seurat object
seurat_dat <- h5ad2seurat_for_cellchat(adata_dir = file.path(dir_val, 'GSE240688_raw/03_GSE240688_remove_double.h5ad'),
                                       add_data = TRUE, add_reduction = TRUE)

# Cluster labels (leiden 1.0)
meta_data <- data.table::fread(file.path(dir_val, 'GSE240688_raw/GSE240688_leiden_1.0.csv'),
                               data.table = FALSE)
rownames(meta_data) <- meta_data[, 1]
meta_data <- meta_data[, -1]

identical(colnames(seurat_dat), rownames(meta_data))

seurat_dat$leiden_1.0 <- meta_data$leiden_1.0
Idents(seurat_dat) <- seurat_dat$leiden_1.0
print(levels(Idents(seurat_dat)))

# ---- AUCell signature score of RORA(+) cell markers ---------------------------

seurat_dat <- qread(file.path(dir_val, '19_seurat_dat.qs'))
Idents(seurat_dat) <- seurat_dat$leiden_1.0

# Top-50 markers of RORA(+) cells from the discovery cohort
top50_cell_markers <- qread('./12_SCENIC/12_top_50_RORA_cells_markers.qs')
rora_markers <- top50_cell_markers$rora_cells$feature

source("./source/compute_module_score.R")

marker_ls <- list(
  RORA_pos_cells = rora_markers)

# Mouse-to-human ortholog conversion
human_marker_ls <- marker_ls %>%
  map(., function(x) {
    a <- convertMouseGeneList(x)
    a <- a[, 2]
  })
names(human_marker_ls) <- names(marker_ls)

qsave(human_marker_ls, file = file.path(dir_val, 'RORA_pos_cells_signature_human.qs'))

human_marker_ls <- qread(file.path(dir_val, 'RORA_pos_cells_signature_human.qs'))

seurat_dat <- ComputeModuleScore(seurat_dat, gene.sets = human_marker_ls,
                                 min.size = 10, cores = 6)

# Violin plot of the signature score by group
p_vlnplot <- VlnPlot(seurat_dat, features = "RORA-pos-cells", group.by = 'group',
                     pt.size = 0.0001,
                     cols = c('#ff7f0e', '#1e76b4') %>% rev()) +
  labs(x = NULL, title = 'RORA(+)_pos_cells_AUCell') +
  theme(axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        legend.position.inside = c(0.8, 0.9),
        legend.position = 'none',
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = .5, size = 18, face = 'bold'))

# NOTE: second VlnPlot call overwrites p_vlnplot with feature "RORA"
# (TF expression instead of the signature score); kept as in the original.
p_vlnplot <- VlnPlot(seurat_dat, features = "RORA", group.by = 'group',
                     pt.size = 0.0001,
                     cols = c('#ff7f0e', '#1e76b4') %>% rev()) +
  labs(x = NULL, title = 'RORA(+)_pos_cells_AUCell') +
  theme(axis.title = element_text(size = 18),
        panel.grid = element_blank(),
        legend.position.inside = c(0.8, 0.9),
        legend.position = 'none',
        legend.text = element_text(size = 16),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = .5, size = 18, face = 'bold'))

# UMAP of external cohort clusters
p_umap <- DimPlot(seurat_dat, cols = ggsci::pal_d3(palette = 'category20')(20),
                  label = TRUE, label.size = 5) +
  ggtitle('Cluster') +
  NoLegend() +
  theme_bw(base_rect_size = 1.5) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 18),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 16),
        legend.box.background = element_rect(color = 'black', linetype = 2, linewidth = 0.8),
        legend.key.width = unit(0.6, 'cm'),
        legend.key.height = unit(0.5, 'cm'),
        legend.title = element_text(size = 16),
        plot.title = element_text(hjust = .5, size = 18, face = 'bold'))

p_combine <- p_umap | p_vlnplot + plot_layout(widths = c(1.2, 1))

ggview:::ggview(p_combine, width = 2.2 * 4, height = 1.2 * 4)

ggsave(file.path(dir_figures, '19_GSE240688_RORA_validation.pdf'), p_combine,
       width = 2.2 * 4, height = 1.2 * 4)

capture.output(sessionInfo(), file = file.path(dir_val, '19_sessionInfo.txt'))
