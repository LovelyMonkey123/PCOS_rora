# CellChat analysis: Control vs PCOS, and RORA(+) GC communication in PCOS

rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(CellChat)
library(reticulate)
library(qs)
library(data.table)
library(future)
library(ggsci)

# Python environment (used upstream for the h5ad conversion helper)
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "celloracle_env")
reticulate::use_condaenv("celloracle_env", required = TRUE)

dir_comm <- './17_cell_communication'
dir_figures <- './figures'

source('./source/analysis_cellchat_ssh.R')

sc <- import('scanpy')

# ---- Prepare the Seurat object for CellChat ------------------------------------

seurat_dat <- h5ad2seurat_for_cellchat('./04_doublet_removal/adata_after_double.h5ad')
qsave(seurat_dat, file = file.path(dir_comm, 'seurat_dat_for_cellchat.qs'))

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

# ---- CellChat per group ---------------------------------------------------------

source('./source/analysis_cellchat_ssh.R')

group_index <- seurat_dat@meta.data$group %>% unique()

seurat_ls <- list()
for (i in group_index) {
  seu <- seurat_dat %>% subset(group == i) %>%
    createCellChat(object = ., group.by = "Cluster_all", assay = "RNA")
  seurat_ls[[i]] <- seu
}

qsave(seurat_ls, file = file.path(dir_comm, 'cellchat_object.list.qs'))

options(future.globals.maxSize = 10000 * 1024^2)
future::plan("multisession", workers = 4)
min.cells <- 10

for (i in seq_along(seurat_ls)) {
  object <- seurat_ls[[i]]

  name_use <- names(seurat_ls)[i]
  dir_output <- file.path(dir_comm, paste0(name_use, '_cellchat_res.qs'))
  cellchat.res <- cellchat_pipeline.ssh(cellchat.input = object,
                                        type = 'triMean',
                                        cellchat.output = dir_output,
                                        category = 'mouse')
  print(paste0(i, ' over!'))
}

# Merge Control and PCOS results
name_use <- c('Control', 'PCOS')
index_c <- file.path(dir_comm, paste0(name_use, '_cellchat_res.qs'))

cellchat_ls <- map(index_c, function(x) { a <- qread(x) })
names(cellchat_ls) <- name_use

cellchat <- mergeCellChat(cellchat_ls,
                          add.names = names(cellchat_ls),
                          cell.prefix = TRUE)

# ---- Interaction number/strength comparison between groups ----------------------

source('./source/analysis_cellchat_ssh.R')

all_counts <- compareInteractions_for_cell_type_ssh(cellchat, measure = 'weight',
                                                    cell_type_index = 'Cluster_all')

plot_a <- compareInteractions_plot_for_cell_type_ssh(input = all_counts,
                                                     measure = 'weight',
                                                     cols = rev(c('#ff7f0e', '#1f77b4')))

ggsave(file.path(dir_figures, '17_cellchat_interaction_strength.pdf'),
       plot_a, width = 2 * 4, height = 2.3 * 4)

# ---- CellChat on RORA(+) labels (PCOS) ------------------------------------------

sce_dat <- qread(file.path(dir_comm, 'RORA_sce_dat_PCOS_for_cellchat.qs'))

source('./source/analysis_cellchat_ssh.R')

options(future.globals.maxSize = 10000 * 1024^2)
future::plan("multisession", workers = 4)
min.cells <- 10  # filter out cell types with too few cells

cc_dat <- sce_dat %>%
  createCellChat(object = ., group.by = "Cluster_PCOS", assay = "RNA")

name_use <- 'RORA_PCOS'
dir_output <- file.path(dir_comm, paste0(name_use, '_cellchat_res.qs'))
cellchat.res <- cellchat_pipeline.ssh(cellchat.input = cc_dat,
                                      type = 'triMean',
                                      cellchat.output = dir_output,
                                      category = 'mouse')

groupSize <- as.numeric(table(cellchat.res@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat.res@net$count, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat.res@net$weight, vertex.weight = groupSize, weight.scale = TRUE,
                 label.edge = FALSE, title.name = "Interaction weights/strength")

mat <- cellchat.res@net$weight

# Circle plots: interactions to and from the RORA(+) population (index 15)
pdf(file.path(dir_figures, '17_cellchat_RORA_Celltypes_interactions_circos.pdf'),
    width = 4.6 * 2, height = 4.6, onefile = TRUE)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(mat,
                 sources.use = setdiff(1:17, 15),
                 targets.use = 15,
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 edge.weight.max = max(mat),
                 title.name = paste0('Target-', rownames(mat)[15]))
netVisual_circle(mat,
                 sources.use = 15,
                 targets.use = setdiff(1:17, 15),
                 vertex.weight = groupSize,
                 weight.scale = TRUE,
                 edge.weight.max = max(mat),
                 title.name = paste0('Source-', rownames(mat)[15]))
dev.off()

# ---- RORA(+) interaction lollipop plots -----------------------------------------

cells.level <- rownames(mat)
df.net <- reshape2::melt(mat, value.name = "value")
colnames(df.net)[1:2] <- c("source", "target")

# Outgoing: RORA(+) -> other cell types
RORA_source.net <- df.net %>% filter(source == 'RORApos_GCs', target != 'RORApos_GCs') %>%
  arrange(desc(value)) %>%
  mutate(target = factor(.$target, levels = .$target))

# Incoming: other cell types -> RORA(+)
RORA_target.net <- df.net %>% filter(target == 'RORApos_GCs', source != 'RORApos_GCs') %>%
  arrange(desc(value)) %>%
  mutate(source = factor(.$source, levels = .$source))

lollipop_theme <- theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

p_RORA_source <- ggplot(RORA_source.net, aes(x = target, y = value)) +
  geom_segment(aes(x = target, xend = target, y = 0, yend = value),
               color = 'black', linetype = 2) +
  geom_point(size = 5, pch = 21, fill = '#f6bea8', color = 'black') +
  labs(title = 'Rora(+)_pos_cells --> Other cell types', x = NULL,
       y = "Interaction Weights") +
  scale_y_continuous(expand = c(0.05, 0)) +
  lollipop_theme

p_RORA_target <- ggplot(RORA_target.net, aes(x = source, y = value)) +
  geom_segment(aes(x = source, xend = source, y = 0, yend = value),
               color = 'black', linetype = 2) +
  geom_point(size = 5, pch = 21, fill = '#b3b6cf', color = 'black') +
  labs(title = 'Other cell types-->Rora(+)_pos_cells', x = NULL,
       y = "Interaction Weights") +
  scale_y_continuous(expand = c(0.05, 0)) +
  lollipop_theme

p.target.source_RORA <- p_RORA_source / p_RORA_target

ggsave(file.path(dir_figures, '17_cellchat_RORA_Celltypes_interactions_lolliplot.pdf'),
       p.target.source_RORA, width = 5.6, height = 10.8)

# ---- Outgoing signaling role heatmap --------------------------------------------

pdf(file.path(dir_figures, '17_cellchat_RORA_Celltype_interactions_ht.pdf'),
    width = 6, height = 8)
ht1 <- netAnalysis_signalingRole_heatmap2(cellchat.res, pattern = "outgoing",
                                          width = 8, height = 14,
                                          color.heatmap.use = colorRampPalette(c('white', '#2166ac'))(100))
ht1
dev.off()

# ---- MK signaling pathway ---------------------------------------------------------

pathways.show <- c("MK")

p_bubble_MK <- netVisual_bubble(cellchat.res,
                                sources.use = c(1:14, 16:17),
                                targets.use = 15,
                                remove.isolate = TRUE,
                                signaling = pathways.show) +
  theme(axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

ggsave(file.path(dir_figures, '17_cellchat_RORA_MK_path_bubble.pdf'),
       p_bubble_MK, width = 5.2, height = 4.2)

# L-R pair contribution to the MK signaling
data_contributation <- netAnalysis_contribution(cellchat.res, signaling = pathways.show,
                                                sources.use = 16, targets.use = 15,
                                                return.data = TRUE)[[1]]

data_contributation$name <- factor(data_contributation$name,
                                   levels = rev(data_contributation$name %>% as.character()))

p_contribution <- ggplot(data_contributation, aes(x = contribution, y = name)) +
  geom_bar(stat = "identity", width = .5, linetype = 2, linewidth = .6,
           color = 'black', fill = '#b2a0ca') +
  labs(title = 'Contribution of each L-R pair', x = 'Relative contribution', y = NULL) +
  scale_x_continuous(expand = c(0, 0)) +
  lollipop_theme

p_contribution

ggsave(file.path(dir_figures, '17_cellchat_RORA_MK_path_contribution.pdf'),
       p_contribution, width = 5.2, height = 3)

# Expression of the MK ligand-receptor genes
plot_a <- dotPlot(sce_dat,
                  features = c('Mdk', 'Ncl', 'Sdc1', 'Sdc4', 'Itga6', 'Itga1', 'Sdc2'),
                  group.by = 'Cluster_PCOS', rotation = FALSE) +
  theme(axis.title = element_text(size = 18),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(color = c('#d62728', rep('#2323f7', 6)), size = 14),
        panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

ggsave(file.path(dir_figures, '17_cellchat_RORA_MK_path_dotplot.pdf'),
       plot_a, width = 5.4, height = 4.2, bg = 'white')

# Individual L-R pair circle plot (top MK pair)
pairLR.MK <- extractEnrichedLR(cellchat.res, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.MK[1, ]

vertex.receiver <- seq(1, 17)

pdf(file.path(dir_figures, '17_cellchat_Mdk_Sdc1_circleplot.pdf'), width = 4, height = 4)
netVisual_individual(cellchat.res, signaling = pathways.show,
                     pairLR.use = LR.show, vertex.receiver = vertex.receiver)
dev.off()
