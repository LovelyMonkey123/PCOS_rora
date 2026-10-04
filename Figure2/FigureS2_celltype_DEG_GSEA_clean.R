# Major cell type DEG analysis (PCOS vs Control) and GO-BP GSEA visualization

rm(list = ls())
options(stringsAsFactors = FALSE)

library(qs)
library(tidyverse)
library(data.table)
library(presto)
library(clusterProfiler)
library(msigdbr)
library(ggrepel)
library(yulab.utils)

dir_comm <- './17_cell_communication'
dir_figures <- './figures'

# ---- DEG analysis per major cell type (wilcoxauc) -----------------------------

seurat_dat_for_deg <- qread(file.path(dir_comm, 'seurat_dat_for_cellchat.qs'))

fast_FindMarkers <- function(object, group_by = 'group', cluster, cluster_name) {
  object$Cluster_use <- object@meta.data[[cluster_name]]
  seu <- object %>% subset(., Cluster_use == cluster)
  deg_res <- presto::wilcoxauc(seu, group_by = 'group', assay = 'data')
  colnames(deg_res) <- c('gene', 'cluster', 'avg_Expr', 'avg_log2FC', 'U_statistic',
                         'AUC', 'p_val', 'p_val_adj', 'pct.1', 'pct.2')
  deg_res$celltype_label <- cluster
  return(deg_res)
}

source('./source/sc_functions.R')

deg_res_ls <- map(seurat_dat_for_deg$Major_annotation_d %>% levels(), function(x) {
  fast_FindMarkers(seurat_dat_for_deg, group_by = 'group',
                   cluster = x, cluster_name = 'Major_annotation_d') %>%
    filter(cluster == 'PCOS')
})

qsave(deg_res_ls, file = file.path(dir_comm, '17_2_celltype_group_DEG_results.qs'))

# ---- DEG scatter plot: logFC vs diff_pct, top 5 labeled per cell type ---------

plot_dat <- deg_res_ls %>% map(., function(x) {
  a <- x %>% mutate(diff_pct = (pct.1 - pct.2) / 100) %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(labeled_color = 'others')
  a[['labeled_color']][1:5] <- a[['celltype_label']][1]
  a[['labeled_color']][seq(to = length(a[['labeled_color']]), length.out = 5)] <- a[['celltype_label']][1]
  return(a)
})

plot_dat <- plot_dat %>%
  Reduce(rbind, .)

p <- ggplot(plot_dat) +
  ggrastr::rasterise(geom_point(aes(x = diff_pct, y = avg_log2FC, color = labeled_color)),
                     dpi = 500, scale = 0.75) +
  ggrepel::geom_text_repel(aes(x = diff_pct, y = avg_log2FC, label = gene, color = labeled_color),
                           data = plot_dat %>% filter(labeled_color != 'others')) +
  geom_hline(yintercept = c(0.25, -0.25), linetype = 2, linewidth = 0.8, color = 'black') +
  scale_color_manual(values = c('B' = '#1f77b4',
                                'Endothelial' = '#ff7f0e',
                                'Epithelial' = '#279e68',
                                'Granulosa' = '#d62728',
                                'Macrophage' = '#aa40fc',
                                'NK T cell' = '#8c564b',
                                'Neutrophil' = '#e377c2',
                                'Oocyte' = '#b5bd61',
                                'Perivascular' = '#17becf',
                                'Stromal' = '#aec7e8',
                                'Theca' = '#ffbb78',
                                'others' = '#cccccc')) +
  labs(y = 'log2(Fold change)', x = NULL) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-2.5, 2.5)) +
  facet_wrap(~ celltype_label, nrow = 2) +
  theme(panel.background = element_rect(fill = 'transparent', color = 'black', linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = 'none',
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16),
        strip.background = element_blank(),
        strip.text = element_text(size = 16))

ggsave(file.path(dir_figures, '17_celltype_DEG_point_plot.pdf'), p,
       width = 6 * 2.8, height = 12)

# ---- GO-BP GSEA on the ranked gene lists --------------------------------------

# Exclude mitochondrial, ribosomal, and hemoglobin genes
genes_all <- deg_res_ls[[1]][['gene']]
MT_genes <- genes_all[str_starts(genes_all, 'mt-')]
ribo_genes <- genes_all[grep("^(RP[SL]|Rps|Rpl)", genes_all)]
hb_genes <- genes_all[grep("^Hb[^(p)]", genes_all)]

dat_for_function <- deg_res_ls %>% map(., function(x) {
  a <- x %>% mutate(diff_pct = (pct.1 - pct.2) / 100) %>%
    filter(!gene %in% c(MT_genes, ribo_genes, hb_genes)) %>%
    arrange(desc(avg_log2FC))
  return(a)
})

# 1. Ranked gene lists per cell type
gsea.input.ls <- dat_for_function %>% map(., function(x) {
  ranks <- x[['avg_log2FC']]
  names(ranks) <- x[['gene']]
  return(ranks)
})
names_use <- dat_for_function %>% map(., function(x) {
  name_use <- x[['celltype_label']][1]
}) %>% Reduce(c, .)
names(gsea.input.ls) <- names_use

# 2. Gene sets (MSigDB GO-BP, mouse)
pathway_genesets <- msigdbr::msigdbr(species = 'Mus musculus', subcollection = c('GO:BP'))
GOBP_genesets <- pathway_genesets %>% dplyr::select(9, 1)
colnames(GOBP_genesets) <- c('term', 'name')
# Clean up pathway names (drop the GO prefix, title case)
GOBP_genesets[['term']] <- str_split_fixed(GOBP_genesets[['term']], '_', 2)[, 2, drop = TRUE] %>%
  str_to_title() %>% gsub('_', ' ', .)

# 3. Run GSEA per cell type
GOBP.res.ls <- list()
for (i in seq_along(gsea.input.ls)) {
  ranks <- gsea.input.ls[[i]]
  gsea.res <- GSEA(ranks, TERM2GENE = GOBP_genesets,
                   minGSSize = 10, maxGSSize = 1000, seed = 1234)
  name_out <- names(gsea.input.ls)[i]
  GOBP.res.ls[[name_out]] <- gsea.res
}

GOBP.res.ls_a <- GOBP.res.ls %>%
  map(., function(x) {
    a <- x@result %>% arrange(desc(NES))
    return(a)
  })

qsave(GOBP.res.ls_a, file = file.path(dir_comm, '17_2_celltype_group_GSEA_results.qs'))

GOBP.res.ls_a <- qread(file.path(dir_comm, '17_2_celltype_group_GSEA_results.qs'))

# ---- GSEA heatmaps: top up (PCOS) and top down (Control) terms ----------------

top5_up_GOBP.res_a <- seq_along(GOBP.res.ls_a) %>%
  map(., function(x) {
    a <- GOBP.res.ls_a[[x]] %>% head(., 4) %>%
      mutate(cluster = names(GOBP.res.ls_a)[x]) %>%
      dplyr::select(cluster, ID, NES)
    return(a)
  }) %>% Reduce(rbind, .)

plot_data_top5 <- top5_up_GOBP.res_a
ID <- plot_data_top5$ID %>%
  as.character() %>%
  map_chr(., function(char.use) {
    ifelse(nchar(char.use) > 60, yulab.utils::str_wrap(char.use, 60), char.use)
  })
plot_data_top5$ID <- ID
levels_use <- plot_data_top5$ID %>% unique()
plot_data_top5[['ID']] <- factor(plot_data_top5[['ID']], levels = levels_use)

p_gsea_a <- ggplot(data = plot_data_top5, aes_string(x = 'cluster', y = 'ID', fill = 'NES')) +
  geom_tile(color = 'white', width = 1, height = 1) +
  labs(y = NULL, x = NULL, title = 'GSEA for PCOS') +
  scale_fill_gradientn(colors = colorRampPalette(colors = c('#f9f3f7', '#ff9638'))(100),
                       limits = c(0, 2.7)) +
  theme(panel.background = element_rect(fill = 'transparent', color = 'black', linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.1, 0.9),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = 'black', linetype = 2, linewidth = 0.8),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 20))

top5_down_GOBP.res_a <- seq_along(GOBP.res.ls_a) %>%
  map(., function(x) {
    a <- GOBP.res.ls_a[[x]] %>% tail(., 4) %>%
      mutate(cluster = names(GOBP.res.ls_a)[x]) %>%
      dplyr::select(cluster, ID, NES)
    return(a)
  }) %>% Reduce(rbind, .)

plot_data_top5 <- top5_down_GOBP.res_a
ID <- plot_data_top5$ID %>%
  as.character() %>%
  map_chr(., function(char.use) {
    ifelse(nchar(char.use) > 60, yulab.utils::str_wrap(char.use, 60), char.use)
  })
plot_data_top5$ID <- ID
levels_use <- plot_data_top5$ID %>% unique()
plot_data_top5[['ID']] <- factor(plot_data_top5[['ID']], levels = levels_use)

p_gsea_b <- ggplot(data = plot_data_top5, aes_string(x = 'cluster', y = 'ID', fill = 'NES')) +
  geom_tile(color = 'white', width = 1, height = 1) +
  labs(y = NULL, x = NULL, title = 'GSEA for Control') +
  scale_fill_gradientn(colors = colorRampPalette(colors = c('#468fc1', '#f9f3f7'))(100),
                       limits = c(-2.5, 0)) +
  scale_y_discrete(position = 'left') +
  theme(panel.background = element_rect(fill = 'transparent', color = 'black', linewidth = 0.8),
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 18),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.1, 0.9),
        legend.text = element_text(size = 14),
        legend.box.background = element_rect(color = 'black', linetype = 2, linewidth = 0.8),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 20))

ggsave(file.path(dir_figures, '17_celltype_GSEA_heatmap_plot.pdf'),
       p_gsea_b + p_gsea_a, width = 24, height = 13)
