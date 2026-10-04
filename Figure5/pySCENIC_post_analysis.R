# ==============================================================================
# SCENIC downstream analysis: granulosa cell regulons (PCOS vs Control)
# ==============================================================================


# ==============================================================================
# Part 1. Package loading
# ==============================================================================

suppressPackageStartupMessages({
  # SCENIC / single-cell core
  library(SCENIC)
  library(AUCell)
  library(Seurat)
  library(SeuratObject)
  library(qs)
  library(presto)

  # Data manipulation
  library(tidyverse)
  library(data.table)
  library(reshape2)

  # Plotting
  library(ggplot2)
  library(ComplexHeatmap)
  library(colorRamp2)
  library(ggsci)
  library(corrplot)
  library(ggraph)
  library(tidygraph)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(scales)
  library(tidydr)
})

# Custom helper functions (project-local sources; keep under version control)
source('./source/SCENIC_function.R') 
source('./source/sc_functions.R')    

# Global parameters
set.seed(123)
n_cores <- 10

dir_scenic <- './12_SCENIC'
dir_scenic_out <- './12_SCENIC/output'
dir_comm <- './17_cell_communication'
dir_figures <- './figures'


# ==============================================================================
# Part 2. Data processing
# ==============================================================================

# ---- 2.1 Binarize regulon AUC activity ---------------------------------------
seu <- qread(file.path(dir_scenic_out, 'gc_annotation_seu.qs'))
regulon_auc <- seu@assays[["AUCell"]]@data %>% as.matrix()

# Explore an AUC threshold per regulon (top 1% of cells as default cut)
threshold_res <- AUCell_exploreThresholds(regulon_auc, thrP = 0.01,
                                          plotHist = TRUE, nCores = n_cores)
qsave(threshold_res, file = file.path(dir_scenic_out, '04_threshold_binary_res.qs'))

thresholds <- AUCell::getThresholdSelected(threshold_res)

# Cells passing the threshold per regulon
regulonsCells <- setNames(lapply(names(thresholds), function(x) {
  trh <- thresholds[[x]]$selected
  names(which(regulon_auc[gsub('.aucThr', '', x), ] > trh))
}), names(thresholds))

# Melt to long format and cast to a regulon x cell binary matrix
regulonActivity <- reshape2::melt(regulonsCells)
binaryRegulonActivity <- t(table(regulonActivity[, 1], regulonActivity[, 2]))
class(binaryRegulonActivity) <- "matrix"
rownames(binaryRegulonActivity) <- gsub('.aucThr', '', rownames(binaryRegulonActivity))

qsave(binaryRegulonActivity, file = file.path(dir_scenic_out, '04_SCENIC_binary_mtx.qs'))

# ---- 2.2 Cluster cells on regulon AUC PCA and compute RSS --------------------
seu <- qread(file.path(dir_scenic_out, 'gc_annotation_seu.qs'))

# Re-cluster granulosa cells using the regulon-AUC PCA (pcaRAS)
seu <- seu %>%
  FindNeighbors(reduction = 'pcaRAS') %>%
  FindClusters(resolution = 1.0)

binaryRegulonActivity <- qread(file.path(dir_scenic_out, '04_SCENIC_binary_mtx.qs'))
regulon_auc <- seu@assays[["AUCell"]]@data %>% as.matrix()

# Regulon specificity score (RSS) per granulosa-cell subtype
rss_cellType <- calcRSS_ssh(regulon_auc,
                            cellAnnotation = seu@meta.data[['Gr_subannotation']]) %>%
  as.data.frame()

# Top-20 regulons per cell type, restricted to binary-active regulons
regulon_celltype <- sapply(names(rss_cellType), function(x) {
  a <- rss_cellType[, x, drop = TRUE]
  names(a) <- rownames(rss_cellType)
  a <- a[intersect(names(a), rownames(binaryRegulonActivity))]
  b <- a %>% sort(., decreasing = TRUE) %>% names(.) %>% head(20)
})
qsave(list(rss_cellType = rss_cellType, regulon_celltype = regulon_celltype),
      file = file.path(dir_scenic, '12_SCENIC_rss_and_top20_regulons.qs'))

# ---- 2.3 Differential regulon activity between PCOS and Control --------------
# Proportion of active cells per regulon within each group; a regulon is
# differential if |PCOS - Control| active proportion exceeds the cutoff.

# FIXME: this function reads the global object `seu`; pass it explicitly if the
# function is reused outside this script.
celltype_topn_regulon_ssh <- function(cell_type, binaryRegulonActivity,
                                      regulon_celltype, cutoff = 0.2) {
  regulons <- regulon_celltype[, cell_type, drop = TRUE]  # subtype-specific regulons
  binary.mtx1 <- binaryRegulonActivity %>% .[regulons, ] %>% t() %>% as.data.frame()

  # Keep cells of the target subtype, split by group
  group_use <- seu[['group']][seu[['Gr_subannotation']] == cell_type]
  binary.mtx.cell.ls <- binary.mtx1 %>%
    .[seu[['Gr_subannotation']] == cell_type, ] %>%
    split(., as.character(group_use))

  # Active proportion per regulon within one group
  one.cell.significant.tfs <- function(x) {
    prop_tf <- apply(x, 2, function(tf) { sum(tf) / length(tf) })
    return(prop_tf)
  }
  res <- map(binary.mtx.cell.ls, one.cell.significant.tfs)

  dea <- abs(res[['PCOS']] - res[['Control']]) > cutoff
  celltype_dea_regulon <- regulons[dea]
  cat(paste0(c(paste0(cell_type, ':', collapse = ' '), celltype_dea_regulon), '\n'))

  return(celltype_dea_regulon)
}


dea_regulon_celltype <- map(paste0('GC_', 0:5), function(x) {
  res <- celltype_topn_regulon_ssh(cell_type = x,
                                   binaryRegulonActivity = binaryRegulonActivity,
                                   regulon_celltype = regulon_celltype,
                                   cutoff = 0.15,
                                   return.dea = TRUE)
})

# Merge per-subtype results into one data frame
dea_df <- dea_regulon_celltype %>% map(., function(x) { a <- x[[2]] }) %>% Reduce(rbind, .)

qsave(dea_regulon_celltype, file = file.path(dir_scenic, '12_SCENIC_celltype_DEA.qs'))
qsave(dea_df, file = file.path(dir_scenic, '12_SCENIC_celltype_DEA_df.qs'))

# ---- 2.4 Heatmap input matrices ----------------------------------------------
dea_df[['index']] <- paste0('Regulon', seq_len(nrow(dea_df)))
tf_select <- dea_df[['Regulon']]

# Long-format binary matrix (cells x regulons), restricted to selected regulons
bianry_df <- binaryRegulonActivity %>% as.data.frame() %>% t() %>% as.data.frame()
mtx_binary <- map(tf_select, function(x) { bianry_df %>% .[, x, drop = FALSE] }) %>%
  Reduce(cbind, .)
colnames(mtx_binary) <- paste0('Regulon', seq_len(ncol(mtx_binary)))
mtx_binary <- mtx_binary %>% t()  # regulons x cells

# ---- 2.5 Correlation of regulon activities -----------------------------------
mtx_ac <- seu@assays$AUCell@data %>% t()

regulon_plot <- dea_regulon_celltype %>% map(., function(x) { x[[2]] }) %>%
  Reduce(rbind, .) %>%
  filter(type != 'none') %>%
  mutate(color = ifelse(type == 'Up', '#ff7f0e', '#1f77b4'))

cor_pcos <- cor(mtx_ac[seu$group == 'PCOS', regulon_plot$Regulon])
cor_con  <- cor(mtx_ac[seu$group != 'PCOS', regulon_plot$Regulon])
qsave(list(cor_pcos = cor_pcos, cor_con = cor_con),
      file = file.path(dir_scenic, '12_SCENIC_regulon_correlations.qs'))

# ---- 2.6 pySCENIC network input (TF -> target importance) --------------------
data <- LoadpySCENICOutput(
  regulon.gmt = file.path(dir_scenic_out, 'GC.regulons.gmt'),
  adj.mat.file = file.path(dir_scenic_out, '01-step1_adj.tsv'))
data <- subset(data, importance > 1)

# ---- 2.7 Proportion of Rora(+) cells --------------------------------
# Binarized Rora/Prdm1 regulon status appended to metadata
seu@meta.data[['Rora_cells']] <- binaryRegulonActivity['Rora(+)', ]

# Per-subtype, per-group fraction of Rora(+) cells
plot_rora <- seu@meta.data %>%
  dplyr::select(Rora_cells, Gr_subannotation, group) %>%
  mutate(group_use = paste0(group, '.', Gr_subannotation)) %>%
  dplyr::select(-c(Gr_subannotation, group)) %>%
  group_by(group_use) %>%
  mutate(percent = 100 * sum(Rora_cells) / n()) %>%
  ungroup() %>%
  dplyr::select(-1) %>%
  distinct(group_use, .keep_all = TRUE) %>%
  separate(., group_use, into = c('group', 'Gr_subannotation'), sep = '\\.') %>%
  mutate(label_percent = round(percent, 2))

qsave(list(plot_rora = plot_rora),
      file = file.path(dir_scenic, '12_SCENIC_regulon_pos_cell_fractions.qs'))

# ---- 2.8 DEGs of Rora(+) cells (wilcoxauc) ----------------------
res_deg_rora <- presto::wilcoxauc(seu, group_by = 'Rora_cells')
qsave(res_deg_rora, file = file.path(dir_scenic, '12_DEA_for_RORA_cells.qs'))

res_deg_rora_filter <- res_deg_rora %>%
  mutate(Ratio = round(pct_in / pct_out, 3), pct.fc = pct_in - pct_out) %>%
  filter(padj < 0.05, pct.fc > 25, logFC > 0.25) %>%
  filter(group == 1)

# ---- 2.9 Rora regulon targets overlapping DEGs (for Metascape) ---------------
target_genes_Rora <- data %>% filter(TF == 'Rora') %>% arrange(desc(importance))
de_target_genes <- intersect(res_deg_rora_filter$feature, target_genes_Rora$target)
qsave(de_target_genes, file = file.path(dir_scenic, '12_DEA_target_genes_for_Rora.qs'))

# Export as a padded matrix for Metascape input (6 columns, filled by row)
n_col <- 6
n_row <- ceiling(length(de_target_genes) / n_col)
mat <- matrix(c(de_target_genes, rep(NA, n_row * n_col - length(de_target_genes))),
              nrow = n_row, ncol = n_col, byrow = TRUE)
df_metascape <- as.data.frame(mat)
df_metascape %>% data.table::fwrite(
  file.path(dir_scenic, '12_rora_regulon_and_DEA_for_metascape_genes.csv'),
  row.names = FALSE, quote = FALSE, na = '')

# ---- 2.10 RORA(+) metadata for cell communication ----------------------------
seu@meta.data$Gr_RORA <- seu@meta.data$Gr_subannotation %>% as.character()
seu@meta.data$Gr_RORA[seu@meta.data$Rora_cells == 1] <- 'RORApos_GCs'
meta_data <- seu@meta.data[, 'Gr_RORA', drop = FALSE]
fwrite(meta_data, file.path(dir_scenic_out, 'RORA_meta_data.csv'), row.names = TRUE)

# ---- 2.11 Prepare Seurat object for CellChat (PCOS only) ---------------------

seurat_data <- qread(file.path(dir_comm, 'seurat_dat_for_cellchat.qs'))

# Subtype annotation of granulosa cells
Gr_subannotation <- data.table::fread(
  file.path(dir_comm, 'pseudotime_res_all.csv'), data.table = FALSE)
rownames(Gr_subannotation) <- Gr_subannotation[, 1]
Gr_subannotation <- Gr_subannotation[, -1]
Gr_subannotation <- Gr_subannotation %>% dplyr::select(Gr_subannotation)

meta_all <- seurat_data@meta.data[, 'Major_annotation_d', drop = FALSE]
meta_use <- meta_all %>% merge(., Gr_subannotation, by = 0, all.x = TRUE, sort = FALSE)
meta_use[['Cluster_all']] <- meta_use[['Major_annotation_d']] %>% as.character()
# Inject granulosa subtype labels into the major annotation
meta_use[['Cluster_all']][meta_use[['Cluster_all']] == 'Granulosa'] <-
  meta_use$Gr_subannotation[meta_use[['Cluster_all']] == 'Granulosa']
meta_use[['Cluster_all']][is.na(meta_use[['Cluster_all']])] <- 'GGGG'  # red blood cells, removed later

rownames(meta_use) <- meta_use$Row.names
stopifnot(identical(meta_use$Row.names, rownames(seurat_data@meta.data)))
# merge() reorders rows; restore the original cell order
meta_use <- meta_use[rownames(seurat_data@meta.data), ]

seurat_data@meta.data$Cluster_all <- meta_use$Cluster_all
seurat_data <- seurat_data %>% subset(Cluster_all != 'GGGG')

sce_dat <- seurat_data %>% subset(group == 'PCOS')

# Merge RORA(+) labels into the PCOS object
meta_data <- fread(file.path(dir_scenic_out, 'RORA_meta_data.csv'), data.table = FALSE)
rownames(meta_data) <- meta_data[, 1]
meta_data <- meta_data[, -1, drop = FALSE]

sce_dat_meta <- sce_dat@meta.data[, 'Cluster_all', drop = FALSE]
meta_filter <- sce_dat_meta %>% merge(., meta_data, by = 0, all.x = TRUE, sort = FALSE)
meta_filter[['Cluster_PCOS']] <- meta_filter[['Gr_RORA']] %>% as.character()
meta_filter[['Cluster_PCOS']][is.na(meta_filter[['Cluster_PCOS']])] <-
  meta_filter$Cluster_all[is.na(meta_filter[['Cluster_PCOS']])]

rownames(meta_filter) <- meta_filter$Row.names
stopifnot(identical(meta_filter$Row.names, rownames(sce_dat@meta.data)))
meta_filter <- meta_filter[rownames(sce_dat@meta.data), ]
stopifnot(identical(meta_filter$Row.names, rownames(sce_dat@meta.data)))

sce_dat@meta.data$Cluster_PCOS <- meta_filter$Cluster_PCOS
qsave(sce_dat, file = file.path(dir_comm, 'RORA_sce_dat_PCOS_for_cellchat.qs'))

# ==============================================================================
# Part 3. Visualization
# Reads back the objects saved in Part 2; no data processing here.
# ==============================================================================

# Shared color scales
colors_cluster <- c('#1f77b4', '#ff7f0e', '#2c9e6b', '#d62728', '#aa42fc', '#8c564b')
names(colors_cluster) <- paste0("GC_", 0:5)
groupcolor <- c('#1f77b4', '#ff7f0e')
names(groupcolor) <- c('Control', 'PCOS')

# ---- 3.1 RSS dot plots per subtype -------------------------------------------
rss_obj <- qread(file.path(dir_scenic, '12_SCENIC_rss_and_top20_regulons.qs'))
rss_cellType <- rss_obj$rss_cellType

plot.list <- lapply(paste0("GC_", 0:5), function(clu) {
  plotRSS_ssh(rss_cellType, clu, regulons.topn = 20,
              col.use = colors_cluster[clu] %>% as.character())
})
ggsave(file.path(dir_figures, '12_SCENIC_RSSplot.pdf'),
       cowplot::plot_grid(plotlist = plot.list, ncol = 3),
       width = 3 * 5, height = 2 * 5)

# ---- 3.2 Binarized regulon activity heatmap ----------------------------------
seu <- qread(file.path(dir_scenic_out, 'gc_annotation_seu.qs'))
binaryRegulonActivity <- qread(file.path(dir_scenic_out, '04_SCENIC_binary_mtx.qs'))
dea_regulon_celltype <- qread(file.path(dir_scenic, '12_SCENIC_celltype_DEA.qs'))
dea_df <- qread(file.path(dir_scenic, '12_SCENIC_celltype_DEA_df.qs'))

tf_select <- dea_df[['Regulon']]
bianry_df <- binaryRegulonActivity %>% as.data.frame() %>% t() %>% as.data.frame()
mtx_binary <- map(tf_select, function(x) { bianry_df %>% .[, x, drop = FALSE] }) %>%
  Reduce(cbind, .)
colnames(mtx_binary) <- paste0('Regulon', seq_len(ncol(mtx_binary)))
mtx_binary <- mtx_binary %>% t()

# Order columns: Control subtypes left-to-right, PCOS mirrored, each sorted by
# subtype then by regulon-AUC clustering (AUCell_snn_res.1)
meta <- seu@meta.data
meta$ID <- rownames(meta)
meta_PCOS <- subset(meta, group == 'PCOS')
meta_Control <- subset(meta, group == 'Control')

meta_Control <- meta_Control %>%
  mutate(celltype = fct_relevel(Gr_subannotation, c(paste0('GC_', 0:5))),
         Regulon_cluster = fct_relevel(AUCell_snn_res.1, as.character(0:17))) %>%
  arrange(celltype, Regulon_cluster)
meta_PCOS <- meta_PCOS %>%
  mutate(celltype = fct_relevel(Gr_subannotation, rev(c(paste0('GC_', 0:5)))),
         Regulon_cluster = fct_relevel(AUCell_snn_res.1, rev(as.character(0:17)))) %>%
  arrange(celltype, Regulon_cluster)

data_plot <- mtx_binary[, c(meta_Control$ID, meta_PCOS$ID)]

annotation_col <- data.frame(
  celltype = c(meta_Control$celltype, meta_PCOS$celltype),
  Regulon_cluster = c(meta_Control$Regulon_cluster, meta_PCOS$Regulon_cluster),
  group = c(meta_Control$group, meta_PCOS$group))
row.names(annotation_col) <- colnames(data_plot)

annotation_row <- data.frame(
  celltype = factor(dea_df$cluster, paste0('GC_', 0:5)),
  DE_score = dea_df$de_score,
  Type = factor(dea_df$type, c('Down', 'none', 'Up')))
row.names(annotation_row) <- rownames(data_plot)

Regulon_cluster_color <- ggsci::pal_d3(palette = 'category20')(20)[1:18]
names(Regulon_cluster_color) <- 0:17

top_anno <- HeatmapAnnotation(df = annotation_col, border = FALSE,
                              show_annotation_name = FALSE,
                              col = list(group = groupcolor,
                                         Regulon_cluster = Regulon_cluster_color,
                                         celltype = colors_cluster))
col_DE_score <- colorRamp2(c(-0.5, 0, 0.5), c('#0072B2', "white", "#C08372"))
left_anno <- rowAnnotation(df = annotation_row, border = FALSE,
                           show_annotation_name = FALSE,
                           col = list(celltype = colors_cluster,
                                      Type = c('Down' = '#3483ba', 'none' = 'grey',
                                               'Up' = '#e46b4f'),
                                      DE_score = col_DE_score))

p <- Heatmap(as.matrix(data_plot),
             cluster_rows = FALSE, cluster_columns = FALSE,
             show_column_names = FALSE, show_row_names = FALSE,
             column_title = NULL,
             heatmap_legend_param = list(title = 'Binarized AUC',
                                         labels = c('1', '0'),
                                         labels_gp = gpar(fontsize = 10),
                                         border = 'black'),
             col = colorRampPalette(colors = c('#D1D1D1', 'blue'))(2),
             top_annotation = top_anno, use_raster = FALSE,
             column_split = annotation_col$group,
             left_annotation = left_anno)

pdf(file.path(dir_figures, '12_heatmap_binary_SCENIC.pdf'), width = 8, height = 8)
p
dev.off()

# ---- 3.3 Regulon labels flanking the heatmap ---------------------------------
dea_regulon_celltype_df <- dea_regulon_celltype %>% map(., function(x) {
  x[[2]] %>% filter(type != 'none')
}) %>% Reduce(rbind, .)

p1 <- ggplot(data.frame(x = c(rep(1.1, nrow(dea_regulon_celltype_df)),
                              rep(1.2, nrow(dea_regulon_celltype_df)),
                              rep(1.3, nrow(dea_regulon_celltype_df))),
                        y = rep(1:nrow(dea_regulon_celltype_df), 3))) +
  geom_tile(aes(x = x, y = y), fill = "white") +
  geom_text(data = data.frame(x = c(1.1, 1.2), y = c(2.5, 2.5),
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_5') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_5') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_text(data = data.frame(x = c(rep(c(1.1, 1.2, 1.3), 3)[1:4]),
                              y = c(6.25, 8.75, 8.75, 8.75),
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_4') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_4') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_text(data = data.frame(x = c(rep(c(1.1, 1.2, 1.3), 4)[1:10]),
                              y = rep(rev(seq(10.625, 14.375, 1.25)), each = 3)[1:10],
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_3') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_3') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_text(data = data.frame(x = c(rep(c(1.1, 1.2, 1.3), 4)[1:4]),
                              y = c(17.5, 17.5, 17.5, 16.25),
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_2') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_2') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_text(data = data.frame(x = c(rep(c(1.1, 1.2, 1.3), 4)[1:2]),
                              y = c(22.5, 22.5),
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_1') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_1') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_text(data = data.frame(x = c(rep(c(1.1, 1.2, 1.3), 4)[1:8]),
                              y = rep(rev(seq(25.625, 29.375, 1.25)), each = 3)[1:8],
                              label = dea_regulon_celltype_df %>% filter(cluster == 'GC_0') %>% .[['Regulon']],
                              type = dea_regulon_celltype_df %>% filter(cluster == 'GC_0') %>% .[['type']]),
            aes(x = x, y = y, label = label, color = type), size = 3.2, hjust = 1) +
  geom_hline(yintercept = seq(5, 25, 5), linetype = 2, linewidth = .5) +
  scale_color_manual(values = c('#3483ba', '#ff7f0e')) +
  theme(legend.position = "none",
        axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank(), panel.background = element_blank(),
        plot.background = element_blank(), axis.ticks = element_blank(),
        plot.margin = margin(l = 0)) +
  scale_x_continuous(limits = c(1.0, 1.32)) +
  scale_y_discrete(expand = c(0, 0))
ggsave(file.path(dir_figures, '12_heatmap_SCENIC_regulons.pdf'), p1,
       width = 2.5, height = 6, bg = 'white')

# ---- 3.4 Regulon activity correlation (PCOS / Control) -----------------------
cor_obj <- qread(file.path(dir_scenic, '12_SCENIC_regulon_correlations.qs'))
cor_pcos <- cor_obj$cor_pcos
cor_con <- cor_obj$cor_con
regulon_plot <- dea_regulon_celltype %>% map(., function(x) { x[[2]] }) %>%
  Reduce(rbind, .) %>%
  filter(type != 'none') %>%
  mutate(color = ifelse(type == 'Up', '#ff7f0e', '#1f77b4'))

BuRd <- c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#FFFFFF",
          "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061") %>% rev()

pdf(file.path(dir_figures, '12_SCENIC_corrplot_PCOS.pdf'), width = 8, height = 8)
corrplot::corrplot(cor_pcos, method = 'pie', type = 'upper', is.corr = TRUE,
                   diag = TRUE, col = colorRampPalette(BuRd)(200),
                   tl.pos = 'td', tl.col = regulon_plot$color)
dev.off()

pdf(file.path(dir_figures, '12_SCENIC_corrplot_Control.pdf'), width = 8, height = 8)
corrplot::corrplot(cor_con, method = 'pie', type = 'lower', is.corr = TRUE,
                   diag = TRUE, col = colorRampPalette(BuRd)(200),
                   tl.col = regulon_plot$color)
dev.off()

# ---- 3.5 Subtype-specific TF target networks ( GC_3) --------------------
data <- LoadpySCENICOutput(
  regulon.gmt = file.path(dir_scenic_out, 'GC.regulons.gmt'),
  adj.mat.file = file.path(dir_scenic_out, '01-step1_adj.tsv'))
data <- subset(data, importance > 1)


cluster_select <- dea_regulon_celltype_df %>% filter(cluster == 'GC_3') %>% pull(1)
tfs_auto <- gsub('\\(\\+\\)', '', cluster_select)
tfs <- 'Rora'

targets.show <- data %>% filter(TF %in% tfs) %>%
  group_by(TF) %>% arrange(desc(importance)) %>% slice_head(n = 10) %>% pull(2)
p_GC3 <- RegulonGraphVis(data %>% filter(TF %in% tfs) %>%
                           group_by(TF) %>% arrange(desc(importance)),
                         tf.show = tfs, targets.show = targets.show,
                         prop = NULL, colors = c('#6eacd1', '#eb7422'),
                         edge.color = '#c1c1c1', edge.alpha = 1, n = 10) +
  labs(title = 'GC_3 Specific Regulon')

ggsave(file.path(dir_figures, '12_SCENIC_subtype_specific_regulon.png'),
       p_GC3 , width = 4, height = 4)

# ---- 3.7 GC_3: same panel for Rora -------------------------------------------
regulon_plot <- dea_regulon_celltype[[4]][[1]][2]

p1_pcos <- plotUMAP_highlight(seu %>% subset(group == 'PCOS'), celltype.highlight = 'GC_3',
                              celltype_col = 'Gr_subannotation', col_use = '#ff7f0e')
p1_control <- plotUMAP_highlight(seu %>% subset(group != 'PCOS'), celltype.highlight = 'GC_3',
                                 celltype_col = 'Gr_subannotation', col_use = '#1e76b4')
p1_pcos <- p1_pcos + theme_void() + labs(title = 'GC_3 (PCOS)') + NoLegend() +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5, size = 18, face = 'bold'))
p1_control <- p1_control + theme_void() + labs(title = 'GC_3 (Control)') + NoLegend() +
  theme(axis.line = element_blank(),
        plot.title = element_text(hjust = .5, size = 18, face = 'bold'))

p2 <- plotRegulon(seu, binary_mtx = binaryRegulonActivity, regulon_plot, col_use = '#d62728')
p_regulon_GC3 <- ((p1_pcos / p1_control + plot_layout(guides = 'keep')) | p2) +
  plot_layout(guides = 'collect', widths = c(1, 2))
ggsave(file.path(dir_figures, '12_SCENIC_GC_3_Rora_regulon_umap.pdf'),
       p_regulon_GC3, width = 9, height = 6)

Idents(seu) <- seu$Gr_subannotation
vlnplot_GC_3_regulon <- VlnPlot(seu, features = regulon_plot, idents = 'GC_3',
                                group.by = 'group', pt.size = 0.0001,
                                cols = c('#ff7f0e', '#1e76b4') %>% rev()) +
  labs(y = 'Regulon activity', x = NULL, title = paste0('Regulon: ', regulon_plot)) +
  vln_theme
vlnplot_GC_3_expr <- VlnPlot(seu, features = gsub('(+)', '', regulon_plot, fixed = TRUE),
                             idents = 'GC_3', group.by = 'group', pt.size = 0.0001,
                             cols = c('#ff7f0e', '#1e76b4') %>% rev()) +
  labs(x = NULL, title = gsub('(+)', '', regulon_plot, fixed = TRUE)) +
  vln_theme
ggsave(file.path(dir_figures, '12_SCENIC_GC_3_Rora_regulon_vlnplot.pdf'),
       vlnplot_GC_3_regulon | vlnplot_GC_3_expr, width = 8, height = 6)

# ---- 3.8 Fraction of Rora(+)cells per subtype ----------------------
frac_obj <- qread(file.path(dir_scenic, '12_SCENIC_regulon_pos_cell_fractions.qs'))
plot_rora <- frac_obj$plot_rora
plot_prdm1 <- frac_obj$plot_prdm1

bar_theme <- theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

p_Rora <- ggplot(plot_rora, aes(x = Gr_subannotation, y = percent, fill = group)) +
  geom_bar(stat = "identity", width = .5, linetype = 2, linewidth = .6, color = 'black') +
  labs(title = 'GSE268919 (Mus musculus)', x = NULL, y = "% Rora(+)_pos_cells") +
  ggrepel::geom_text_repel(aes(label = label_percent), size = 4) +
  scale_fill_manual(values = alpha(c('#ff7f0e', '#1e76b4') %>% rev(), .7)) +
  scale_y_continuous(expand = c(0, 0)) +
  bar_theme

ggsave(file.path(dir_figures, '12_SCENIC_regulon_binary_cells_percent_barplot.pdf'),
        p_Rora, width = 4 , height = 4)

# ---- 3.9 RORA(+) cells on UMAP -----------------------------------------------
seu <- qread(file.path(dir_scenic_out, 'gc_annotation_seu.qs'))
meta_rora <- data.table::fread(file.path(dir_scenic_out, 'RORA_meta_data.csv'),
                               data.table = FALSE, header = TRUE)
stopifnot(identical(meta_rora$V1, colnames(seu)))
seu$RORA_pos_cells <- meta_rora$Gr_RORA

p3 <- DimPlot(seu, reduction = "X_umap", group.by = "RORA_pos_cells",
              pt.size = 0.000001, split.by = 'group') +
  scale_color_manual(values = c('#1f77b4', '#ff7f0e', '#2c9e6b', '#d62728',
                                '#aa42fc', '#8c564b', 'black')) +
  tidydr::theme_dr() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        plot.title = element_blank(), legend.text = element_text(size = 14),
        axis.title = element_text(size = 14), strip.text = element_text(size = 16))
ggsave(file.path(dir_figures, '12_SCENIC_Rora_umap_plot.pdf'), p3,
       width = 4.5 * 2, height = 4, bg = 'white')

# ---- Session info for reproducibility ----------------------------------------
capture.output(sessionInfo(), file = file.path(dir_scenic, '12_SCENIC_sessionInfo.txt'))
