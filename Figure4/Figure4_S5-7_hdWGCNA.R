# hdWGCNA analysis of granulosa cell subtypes (PCOS vs Control)
# Ref: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
# Part 1: packages | Part 2: data processing | Part 3: visualization

# ============================== Part 1. Packages ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(qs)
  library(reticulate)
  library(ComplexHeatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(scRNAtoolVis)
  library(igraph)
  library(ggrepel)
  library(ggh4x)
  library(viridis)
  library(Hmisc)
  library(IOBR)
  library(yulab.utils)
  library(data.table)
})

# NOTE: 'Matirx' is loaded as in the original script (used as Matirx::t()).
library(Matirx)

# Python environment for reading the scanpy h5ad object
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "celloracle_env")
reticulate::use_condaenv("celloracle_env", required = TRUE)

set.seed(12345)
n_threads <- 6
dir_hd <- './11_hdWGCNA'
dir_figures <- './figures'


# ============================ Part 2. Data processing =========================

# ---- 2.1 Convert the scanpy AnnData to a Seurat object -----------------------
sc <- import('scanpy')
adata_gr <- sc$read_h5ad('./10_cytotrace2/10_after_hb_Gr_subanndata.h5ad')

meta <- adata_gr$obs
var <- adata_gr$var

# PCA / harmony / UMAP embeddings (python stores cells as rows)
X_pca <- adata_gr$obsm['X_pca']
rownames(X_pca) <- rownames(meta)

X_umap <- adata_gr$obsm['X_umap']
rownames(X_umap) <- rownames(meta)
colnames(X_umap) <- c('UMAP-1', 'UMAP-2')

X_harmony <- adata_gr$obsm['X_pca_harmony']
rownames(X_harmony) <- rownames(meta)

# Raw counts; transpose because python stores genes as columns
counts_data <- adata_gr$layers['counts'] %>% Matirx::t()
colnames(counts_data) <- rownames(meta)
rownames(counts_data) <- rownames(var)

seurat_data <- CreateSeuratObject(counts = counts_data, meta.data = meta)
seurat_data[['X_umap']] <- CreateDimReducObject(embeddings = X_umap, key = 'UMAP_', assay = 'RNA')
seurat_data[['X_pca']] <- CreateDimReducObject(embeddings = X_pca, key = 'PCA_', assay = 'RNA')
seurat_data[['harmony']] <- CreateDimReducObject(embeddings = X_harmony, key = 'harmony_', assay = 'RNA')

# KNN graph from scanpy (stored as 'connectivities' in python, used as RNA_snn here)
RNA_snn <- adata_gr$obsp['connectivities']
rownames(RNA_snn) <- rownames(meta)
colnames(RNA_snn) <- rownames(meta)
RNA_snn <- RNA_snn %>% SeuratObject::as.Graph(.)
seurat_data@graphs[['RNA_snn']] <- RNA_snn

# Normalized data matrix (converted to CsparseMatrix)
data_mtx <- adata_gr$raw$X %>% t() %>% as(., 'CsparseMatrix')
colnames(data_mtx) <- rownames(meta)
rownames(data_mtx) <- rownames(var)
seurat_data[['RNA']]@data <- data_mtx

# Highly variable genes from scanpy
highly_variable_genes <- var$highly_variable
VariableFeatures(seurat_data) <- rownames(var)[highly_variable_genes]

qsave(seurat_data, file.path(dir_hd, '11_gr_subannotation_seurat_data.qs'))

# ---- 2.2 hdWGCNA: setup and metacells ----------------------------------------
theme_set(theme_cowplot())
enableWGCNAThreads(nThreads = n_threads)

sce_sub <- qread(file.path(dir_hd, '11_gr_subannotation_seurat_data.qs'))
Idents(sce_sub) <- sce_sub$Gr_subannotation

sce.obj <- SetupForWGCNA(seurat_obj = sce_sub,
                         gene_select = 'fraction',
                         fraction = 0.1,
                         wgcna_name = 'GC_wgcna')

# Metacells per group x subtype on the harmony embedding
sce.obj <- sce.obj %>%
  MetacellsByGroups(seurat_obj = .,
                    group.by = c('group', 'Gr_subannotation'),
                    reduction = 'harmony',
                    k = 30,
                    max_shared = 10,
                    min_cells = 30,
                    ident.group = c('Gr_subannotation'),
                    mode = 'sum',
                    verbose = TRUE) %>%
  NormalizeMetacells(.)

sce.obj <- sce.obj %>%
  SetDatExpr(seurat_obj = ., assay = 'RNA', slot = 'data')

# ---- 2.3 Soft power selection and network construction -----------------------
sce.obj <- sce.obj %>% TestSoftPowers(seurat_obj = ., networkType = 'signed')
power_table <- GetPowerTable(sce.obj)
qsave(power_table, file = file.path(dir_hd, '11_SoftPowers_table.qs'))

# First power reaching scale-free topology fit R^2 >= 0.8
soft_power <- power_table %>% filter(SFT.R.sq >= 0.8) %>% .[1, 1]

sce.obj <- sce.obj %>%
  ConstructNetwork(seurat_obj = .,
                   soft_power = soft_power,
                   deepSplit = 4,
                   detectCutHeight = 0.995,
                   minModuleSize = 100,
                   mergeCutHeight = 0.2,
                   setDatExpr = FALSE,
                   tom_name = '',
                   overwrite_tom = TRUE)

# ---- 2.4 Module eigengenes and connectivity ----------------------------------
sce.obj <- sce.obj %>% ScaleData(., features = VariableFeatures(sce.obj))

# Harmonized MEs, regressing out 'orig.ident'
sce.obj <- sce.obj %>% ModuleEigengenes(., group.by.vars = "orig.ident", exclude_grey = TRUE)
hMEs <- GetMEs(sce.obj)
MEs <- GetMEs(sce.obj, harmonized = FALSE)

sce.obj <- ModuleConnectivity(sce.obj, harmonized = TRUE, wgcna_name = 'GC_wgcna')
sce.obj <- ResetModuleNames(sce.obj, new_name = "M")

modules <- GetModules(sce.obj)
hub_df <- GetHubGenes(sce.obj, n_hubs = 10)

# Module score from top 50 genes (alternative to hMEs)
sce.obj <- ModuleExprScore(sce.obj, n_genes = 50, method = 'Seurat')

# ---- 2.5 Hub gene expression table for the heatmap ----------------------------
dat_ht <- sce.obj %>%
  Seurat::FetchData(object = ., vars = hub_df$gene_name, slot = 'data') %>%
  t()
anno_row <- hub_df[, 1:2]
rownames(anno_row) <- anno_row[, 1]
anno_col <- sce.obj@meta.data %>% dplyr::select(Gr_subannotation, group, Palantir_pseudotime)

qsave(list(dat_ht = dat_ht, anno_row = anno_row, anno_col = anno_col),
      file = file.path(dir_hd, '11_hdWGCNA_hub_heatmap_input.qs'))

qsave(sce.obj, file = file.path(dir_hd, '11_gr_hdWGCNA_object.qs'))

# ---- 2.6 GO enrichment of top-200 hub genes per module ------------------------
top200_hubgenes <- GetHubGenes(sce.obj, n_hubs = 200, wgcna_name = 'GC_wgcna')

eres.list <- lapply(sort(unique(top200_hubgenes$module)), function(clu) {
  genes <- subset(top200_hubgenes, module == clu)$gene_name
  enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL',
           ont = 'BP', pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 1000) %>%
    clusterProfiler::simplify(., cutoff = 0.7, by = "p.adjust", select_fun = min)
})
names(eres.list) <- sort(unique(top200_hubgenes$module))
qsave(eres.list, file = file.path(dir_hd, '11_hdWGCNA_GO_top200_res.qs'))

# Top 5 GO terms per module (ranked by RichFactor) for plotting
go.df.ls <- lapply(names(eres.list), function(xx) {
  res <- eres.list[[xx]] %>% as.data.frame() %>%
    filter(p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  res$Description <- factor(res$Description, levels = unique(res$Description) %>% rev())
  res
})
names(go.df.ls) <- names(eres.list)

go.df.ls1 <- lapply(go.df.ls, function(xx) {
  xx %>% .[seq_len(min(5, nrow(.))), ] %>% arrange(desc(RichFactor))
})

module_colors <- modules %>% dplyr::select(c(module, color)) %>% dplyr::distinct()
mod_colors <- module_colors$color
names(mod_colors) <- module_colors$module
rowcolor <- mod_colors[-1]
names(rowcolor) <- names(go.df.ls)

qsave(list(go.df.ls1 = go.df.ls1, rowcolor = rowcolor),
      file = file.path(dir_hd, '11_hdWGCNA_GO_barplot_input.qs'))

# ---- 2.7 Module-trait correlation (MEs vs subtype / group) --------------------
sce.obj <- qread(file.path(dir_hd, '11_gr_hdWGCNA_object.qs'))

# FIXME: column positions (57:68 = M1-M11 hMEs) depend on the current metadata
# layout and will break if meta.data columns are added upstream; select the
# module columns by name instead when maintaining this script.
meta_hdwgcna <- sce.obj@meta.data %>% dplyr::select(-c(69:80))
meta_hdwgcna <- meta_hdwgcna %>% dplyr::select(57:68, group, Gr_subannotation)

# Dummy-code traits
GC_meta <- model.matrix(~ 0 + Gr_subannotation, data = meta_hdwgcna) %>% as.data.frame()
colnames(GC_meta) <- gsub('Gr_subannotation', '', colnames(GC_meta))
group_meta <- model.matrix(~ 0 + group, data = meta_hdwgcna) %>% as.data.frame()
colnames(group_meta) <- gsub('group', '', colnames(group_meta))
meta_hdwgcna <- meta_hdwgcna %>% cbind(., GC_meta, group_meta)

cor_indexs2 <- c('Control', 'PCOS', paste0('GC_', 0:5))
cor_index1 <- paste0('M', 1:11)

cor_data <- map(seq_along(cor_index1), function(a) {
  IOBR::batch_cor(meta_hdwgcna, target = cor_index1[a], feature = cor_indexs2,
                  method = 'pearson') %>%
    mutate(Cor_index = cor_index1[a]) %>% as.data.frame()
}) %>% Reduce(rbind, .)

cor_data$sig_names <- cor_data$sig_names %>% factor(., levels = c(paste0('GC_', 0:5), 'Control', 'PCOS'))
cor_data$label <- sprintf('%.2f', cor_data$statistic)
cor_data$label <- ifelse(cor_data$p.value < 0.05, cor_data$label, 'X')
cor_data$Cor_index <- factor(cor_data$Cor_index, levels = paste0('M', 1:11))
qsave(cor_data, file = file.path(dir_hd, '11_hdWGCNA_module_traits_cor.qs'))

# ---- 2.8 Differential module eigengene (DME) analyses -------------------------
run_dme <- function(barcodes1, barcodes2, suffix) {
  res <- FindDMEs(sce.obj, barcodes1 = barcodes1, barcodes2 = barcodes2,
                  test.use = 'wilcox', wgcna_name = 'GC_wgcna',
                  harmonized = TRUE, only.pos = FALSE)
  colnames(res)[2] <- paste0('avg_log2FC.', suffix)
  colnames(res)[5] <- paste0('p_val_adj.', suffix)
  res %>% dplyr::select(6, 2, 5)
}

modules <- GetModules(sce.obj)
colors.use <- as.data.frame(table(Modules$module, Modules$color)) %>%
  filter(Freq != 0) %>% dplyr::select(1:2)
Nr.genes <- table(Modules$module) %>% as.data.frame()
colnames(Nr.genes)[2] <- 'Nr.genes'

# Trajectory 1 (Palantir branch 1 vs all other cells)
group1 <- sce.obj@meta.data %>% subset(Branch1_cells == 'TRUE') %>% rownames
group2 <- sce.obj@meta.data %>% subset(Branch1_cells != 'TRUE') %>% rownames
DMEs.branch1 <- run_dme(group1, group2, 'branch')

# Trajectory 2 (Palantir branch 2 vs all other cells)
group1 <- sce.obj@meta.data %>% subset(Branch2_cells == 'TRUE') %>% rownames
group2 <- sce.obj@meta.data %>% subset(Branch2_cells != 'TRUE') %>% rownames
DMEs.branch2 <- run_dme(group1, group2, 'branch')

# PCOS vs Control
group1 <- sce.obj@meta.data %>% subset(group == 'PCOS') %>% rownames
group2 <- sce.obj@meta.data %>% subset(group == 'Control') %>% rownames
DMEs.group <- run_dme(group1, group2, 'group')

# FIXME: the scatter plot below is labeled 'for GC_0' but this contrast is
# GC_5 vs all other subtypes; confirm which subtype the figure should show.
group1 <- sce.obj@meta.data %>% subset(Gr_subannotation == 'GC_5') %>% rownames
group2 <- sce.obj@meta.data %>% subset(Gr_subannotation != 'GC_5') %>% rownames
DMEs.celltype <- run_dme(group1, group2, 'celltype')

DMEs <- DMEs.celltype %>%
  merge(., DMEs.group, by.x = 1, by.y = 1) %>%
  merge(., colors.use, by.x = 1, by.y = 1) %>%
  merge(., Nr.genes, by.x = 1, by.y = 1)
qsave(DMEs, file = file.path(dir_hd, '11_hdWGCNA.DMEs.qs'))

# ---- 2.9 Pseudotime-binned hMEs (trajectory 1 and 2) --------------------------
# Pseudotime per branch and UMAP coordinates
sce.obj$b1_pseudotime <- ifelse(sce.obj$Branch1_cells, sce.obj$Palantir_pseudotime, NA)
sce.obj$b2_pseudotime <- ifelse(sce.obj$Branch2_cells, sce.obj$Palantir_pseudotime, NA)
sce.obj$UMAP1 <- sce.obj@reductions$X_umap@cell.embeddings[, 1]
sce.obj$UMAP2 <- sce.obj@reductions$X_umap@cell.embeddings[, 2]

bin_hMEs <- function(meta_branch, modules_keep, cur_bins = 20) {
  meta_branch$PCOS_pseudotime <- ifelse(meta_branch$group == 'PCOS', meta_branch$Palantir_pseudotime, NA)
  meta_branch$Control_pseudotime <- ifelse(meta_branch$group != 'PCOS', meta_branch$Palantir_pseudotime, NA)

  pseudotime_col <- c('PCOS_pseudotime', 'Control_pseudotime')
  for (ps_col in pseudotime_col) {
    bin_name <- paste0(ps_col, "_bins_", cur_bins)
    meta_branch[[bin_name]] <- as.numeric(meta_branch[[ps_col]]) %>% Hmisc::cut2(g = cur_bins)
  }

  modules <- GetModules(sce.obj, wgcna_name = 'GC_wgcna')
  mods <- levels(modules$module)
  mods <- mods[mods != "grey"]

  avg_list <- list()
  for (ps_col in pseudotime_col) {
    ps_bin_name <- paste0(ps_col, "_bins_", cur_bins)
    avg_scores <- meta_branch %>% group_by(get(ps_bin_name)) %>%
      dplyr::select(all_of(mods)) %>% summarise_all(median)
    colnames(avg_scores)[1] <- "bin"
    avg_df <- reshape2::melt(avg_scores)
    avg_df$bin <- as.numeric(avg_df$bin)
    avg_df$group <- ps_col
    avg_list[[ps_col]] <- avg_df
  }
  plot_df <- do.call(rbind, avg_list)
  plot_df <- plot_df %>% filter(variable %in% modules_keep)
  plot_df$variable <- factor(plot_df$variable, levels = modules_keep)

  list(meta_branch = meta_branch, plot_df = plot_df)
}

res_b1 <- bin_hMEs(sce.obj@meta.data %>% filter(!is.na(b1_pseudotime)), paste0('M', c(9, 8)))
res_b2 <- bin_hMEs(sce.obj@meta.data %>% filter(!is.na(b2_pseudotime)), paste0('M', c(9, 1, 5, 7, 3, 10, 2)))

module_colors <- modules %>% dplyr::select(c(module, color)) %>% dplyr::distinct()
mod_colors <- module_colors$color
names(mod_colors) <- module_colors$module

qsave(list(plot_df_b1 = res_b1$plot_df, plot_df_b2 = res_b2$plot_df,
           mod_colors = mod_colors,
           b1_pseudotime = res_b1$meta_branch, b2_pseudotime = res_b2$meta_branch),
      file = file.path(dir_hd, '11_hdWGCNA_pseudotime_hMEs.qs'))

# ---- 2.10 M8 hub genes export -------------------------------------------------
hub_df <- GetHubGenes(sce.obj, n_hubs = 200, mods = c('M8'))
qsave(hub_df, file = file.path(dir_hd, '11_hdWGCNA_hub_df_M8.qs'))

# 6-column padded matrix for external tool input (Metascape-style)
n_col <- 6
n_row <- ceiling(nrow(hub_df) / n_col)
mat <- matrix(c(hub_df$gene_name, rep(NA, n_row * n_col - nrow(hub_df))),
              nrow = n_row, ncol = n_col, byrow = TRUE)
df <- as.data.frame(mat)
df %>% data.table::fwrite(file.path(dir_hd, '11_gr_hdWGCNA_M8_top200_hub_genes.csv'),
                          row.names = FALSE, quote = FALSE, na = '')

# ---- 2.11 Rora expression along pseudotime ------------------------------------
b_pseudotime <- sce.obj@meta.data
Rora_gene_exp <- sce.obj@assays$RNA@data['Rora', rownames(b_pseudotime)]
b_pseudotime <- b_pseudotime %>% cbind(., Rora_gene_exp)

cur_bins <- 100
b_pseudotime$PCOS_pseudotime <- ifelse(b_pseudotime$group == 'PCOS', b_pseudotime$Palantir_pseudotime, NA)
b_pseudotime$Control_pseudotime <- ifelse(b_pseudotime$group != 'PCOS', b_pseudotime$Palantir_pseudotime, NA)

pseudotime_col <- c('PCOS_pseudotime', 'Control_pseudotime')
for (ps_col in pseudotime_col) {
  bin_name <- paste0(ps_col, "_bins_", cur_bins)
  b_pseudotime[[bin_name]] <- as.numeric(b_pseudotime[[ps_col]]) %>% Hmisc::cut2(g = cur_bins)
}

avg_list <- list()
for (ps_col in pseudotime_col) {
  ps_bin_name <- paste0(ps_col, "_bins_", cur_bins)
  avg_scores <- b_pseudotime %>% group_by(get(ps_bin_name)) %>%
    dplyr::select(Rora_gene_exp) %>% summarise_all(mean)
  colnames(avg_scores)[1] <- "bin"
  avg_df <- reshape2::melt(avg_scores)
  avg_df$bin <- as.numeric(avg_df$bin)
  avg_df$group <- ps_col
  avg_list[[ps_col]] <- avg_df
}
plot_df_rora <- do.call(rbind, avg_list)

qsave(list(plot_df_rora = plot_df_rora, b_pseudotime = b_pseudotime),
      file = file.path(dir_hd, '11_hdWGCNA_Rora_pseudotime.qs'))


# =========================== Part 3. Visualization ============================

sce.obj <- qread(file.path(dir_hd, '11_gr_hdWGCNA_object.qs'))
Idents(sce.obj) <- sce.obj$Gr_subannotation

modules <- GetModules(sce.obj)
module_colors <- modules %>% dplyr::select(c(module, color)) %>% dplyr::distinct()
mod_colors <- module_colors$color
names(mod_colors) <- module_colors$module

subtype_color <- c('#1f77b4', '#ff7f0e', '#2c9e6b', '#d62728', '#aa42fc', '#8c564b')
names(subtype_color) <- levels(Idents(sce.obj))
group_color <- c('#ff7d0b', '#1e76b3')
names(group_color) <- c('PCOS', 'Control')

# ---- 3.1 Soft power plot ------------------------------------------------------
power_table <- qread(file.path(dir_hd, '11_SoftPowers_table.qs'))
sce.obj <- TestSoftPowers(sce.obj, networkType = 'signed')  # restored to hold power plots
plot_list <- PlotSoftPowers(sce.obj)
ggsave(file.path(dir_figures, '11_hdWGCNA_softpower.pdf'),
       wrap_plots(plot_list, ncol = 2), width = 8, height = 7, bg = 'white')

# ---- 3.2 Module dendrogram ----------------------------------------------------
pdf(file.path(dir_figures, '11_hdWGCNA_GC_Dendrogram.pdf'), width = 6, height = 4.5)
PlotDendrogram(sce.obj, main = 'GC hdWGCNA Dendrogram',
               groupLabels = 'Module colors', saveMar = FALSE, cex.colorLabels = 1)
dev.off()

# ---- 3.3 Hub gene kME plots ---------------------------------------------------
p2 <- PlotKMEs(sce.obj, ncol = 4, n_hubs = 10, text_size = 5, wgcna_name = 'GC_wgcna')
ggsave(file.path(dir_figures, '11_hdWGCNA_module_hub_kME_plot.pdf'),
       p2, width = 3.2 * 4, height = 3.8 * 3, bg = 'white')

# ---- 3.4 hMEs on UMAP ---------------------------------------------------------
plot_list <- ModuleFeaturePlot(sce.obj, reduction = 'X_umap', point_size = 0.5,
                               raster = TRUE, raster_dpi = 300,
                               features = 'hMEs', order = 'shuffle',
                               wgcna_name = 'GC_wgcna')
p.list <- wrap_plots(plot_list, ncol = 6)
ggsave(file.path(dir_figures, '11_hdWGCNA_hMEs_featureplot.pdf'),
       p.list, width = 2.4 * 6, height = 2.4 * 2, bg = 'white')

# ---- 3.5 Module ME dot plots (subtype / group) --------------------------------
MEs <- GetMEs(sce.obj, harmonized = TRUE)
mods <- colnames(MEs)
mods <- mods[mods != 'grey']
sce.obj@meta.data <- cbind(sce.obj@meta.data, MEs)

p1 <- DotPlot(sce.obj, features = mods, group.by = 'Gr_subannotation') +
  coord_flip() + RotatedAxis() +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.key.size = unit(0.6, 'cm'))
p2 <- DotPlot(sce.obj, features = mods, group.by = 'group') +
  coord_flip() + RotatedAxis() +
  scale_color_gradient2(high = 'red', mid = 'grey95', low = 'blue') +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.key.size = unit(0.4, 'cm'))

ggsave(file.path(dir_figures, '11_hdWGCNA_module_cell_type_dotplot.pdf'),
       p1, width = 5.5, height = 8, bg = 'white')
ggsave(file.path(dir_figures, '11_hdWGCNA_module_group_dotplot.pdf'),
       p2, width = 4, height = 8, bg = 'white')

# ---- 3.6 Hub gene heatmap -----------------------------------------------------
ht_input <- qread(file.path(dir_hd, '11_hdWGCNA_hub_heatmap_input.qs'))
dat_ht <- ht_input$dat_ht
hub_df <- qread(file.path(dir_hd, '11_hdWGCNA_hub_df_M8.qs'))  # placeholder guard
hub_df <- GetHubGenes(sce.obj, n_hubs = 10)

a <- scRNAtoolVis::averageHeatmap(object = sce.obj, markerGene = hub_df$gene_name)
anno_col <- data.frame(row.names = colnames(a@matrix), Cell_type = colnames(a@matrix))
top_anno <- HeatmapAnnotation(df = anno_col, show_legend = TRUE, which = 'column',
                              border = FALSE, annotation_label = c('Cell Type'),
                              col = list(Cell_type = subtype_color),
                              show_annotation_name = TRUE)
ht.my <- ComplexHeatmap::Heatmap(matrix = a@matrix, name = 'Expression',
                                 col = colorRampPalette(c('#2e92f7', 'white', "#da0612"))(100),
                                 border = TRUE,
                                 show_row_names = TRUE, row_names_gp = gpar(fontsize = 10),
                                 row_names_side = 'right',
                                 show_column_names = TRUE, column_names_rot = 45,
                                 column_names_side = 'top', column_names_gp = gpar(fontsize = 14),
                                 cluster_rows = FALSE, cluster_columns = FALSE,
                                 row_split = hub_df$module, row_title_rot = 0,
                                 column_split = anno_col$Cell_type, column_title = NULL,
                                 row_title_gp = gpar(fontsize = 14),
                                 top_annotation = top_anno)
pdf(file.path(dir_figures, '11_hdWGCNA_Module_celltype_ht.pdf'), width = 5, height = 16)
draw(ht.my, merge_legend = TRUE, heatmap_legend_side = 'right',
     annotation_legend_side = 'right')
dev.off()

# ---- 3.7 GO enrichment bar plots per module -----------------------------------
go_input <- qread(file.path(dir_hd, '11_hdWGCNA_GO_barplot_input.qs'))
go.df.ls1 <- go_input$go.df.ls1
rowcolor <- go_input$rowcolor

enrich.bar.ls <- lapply(names(go.df.ls1), function(M.name) {
  mytheme <- theme(legend.position = 'none',
                   plot.title = element_text(size = 14, face = 'bold', hjust = .5),
                   axis.title = element_text(size = 13),
                   axis.text = element_text(size = 11),
                   axis.ticks = element_blank(),
                   axis.text.x = element_blank(),
                   axis.line.y.left = element_blank(),
                   axis.ticks.y = element_blank())
  go.df <- go.df.ls1[[M.name]]
  Description <- go.df$Description %>% as.character() %>%
    map_chr(., function(char.use) {
      ifelse(nchar(char.use) > 60, yulab.utils::str_wrap(char.use, 60), char.use)
    })
  go.df$Description <- Description
  go.df$Description <- factor(go.df$Description, levels = unique(go.df$Description) %>% rev())

  ggplot() +
    geom_bar(data = go.df, aes(x = RichFactor, y = Description),
             color = rowcolor[M.name], width = 0.8, stat = 'identity', fill = 'transparent') +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.y = element_blank()) +
    geom_text(data = go.df, aes(x = 0, y = Description, label = Description),
              size = 5, color = 'black', hjust = -0.01) +
    labs(x = 'GeneRatio', y = NULL, title = M.name) +
    mytheme
})
qsave(enrich.bar.ls, file = file.path(dir_hd, '11_hdWGCNA_all_module_enrich_plot_ls.qs'))

p.comb <- cowplot::plot_grid(plotlist = enrich.bar.ls, ncol = 3, byrow = TRUE, nrow = 4)
ggsave(file.path(dir_figures, '11_hdWGCNA_module_GO_function.pdf'),
       p.comb, width = 6 * 3, height = 3 * 4)

# ---- 3.8 Module network plots -------------------------------------------------
dir.create(file.path(dir_figures, 'hdWGCNA'), showWarnings = FALSE)
ModuleNetworkPlot(sce.obj, mods = 'all',
                  outdir = file.path(dir_figures, 'hdWGCNA'))

options(future.globals.maxSize = 100 * 1024^3)
HubGeneNetworkPlot(sce.obj, n_hubs = 5, n_other = 5,
                   edge_prop = 0.75, mods = 'all')

# ---- 3.9 Hub gene UMAP --------------------------------------------------------
sce.obj <- RunModuleUMAP(sce.obj, n_hubs = 10, n_neighbors = 15, min_dist = 0.2,
                         exclude_grey = TRUE)
pdf(file.path(dir_figures, '11_hdWGCNA_ModuleUMAP.pdf'), width = 10, height = 10, bg = 'white')
ModuleUMAPPlot(sce.obj, edge.alpha = 0.25, sample_edges = TRUE,
               edge_prop = 0.2, label_hubs = 5, keep_grey_edges = FALSE,
               vertex.label.cex = 0.5)
dev.off()
qsave(sce.obj, file.path(dir_hd, '11_gr_hdWGCNA_object.qs'))

# ---- 3.10 Module-trait correlation heatmap ------------------------------------
cor_data <- qread(file.path(dir_hd, '11_hdWGCNA_module_traits_cor.qs'))
p_cor <- ggplot(data = cor_data, aes(x = sig_names, y = Cor_index, fill = statistic)) +
  geom_tile(color = 'black', width = 0.9, height = 0.9, linewidth = 0.4) +
  geom_text(aes(label = label), size = 4.5) +
  scale_fill_gradient2(low = '#4177b9', high = '#d76364') +
  labs(x = 'Traits', y = 'hdWGCNA Modules', fill = 'Cor') +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16, hjust = 1, vjust = 1, angle = 45),
        axis.title = element_text(size = 16), axis.ticks = element_blank(),
        legend.key.size = unit(.3, 'cm'))
ggsave(file.path(dir_figures, '11_hdWGCNA_module_traits_cor.pdf'),
       p_cor, width = 6, height = 8)

# ---- 3.11 Branch DME bar plots ------------------------------------------------
plot_dme <- function(res, xlabel, outfile) {
  DMEs_plot <- res %>%
    merge(., colors.use, by.x = 1, by.y = 1) %>%
    merge(., Nr.genes, by.x = 1, by.y = 1)

  logFC_cutoff <- 0.5
  pval_cutoff <- 0.05
  DMEs_plot$branch <- ifelse(DMEs_plot$avg_log2FC.branch > logFC_cutoff &
                               DMEs_plot$p_val_adj.branch < pval_cutoff, '*', '')
  DMEs_plot <- DMEs_plot %>%
    arrange(desc(avg_log2FC.branch)) %>%
    dplyr::mutate(module = factor(as.character(.$module), levels = rev(as.character(.$module))))

  p <- ggplot(DMEs_plot, aes(avg_log2FC.branch, module, fill = module)) +
    geom_col(color = 'black', linetype = 2, linewidth = 0.5) +
    geom_text(aes(label = branch, color = module), size = 6, nudge_x = 0.1) +
    geom_vline(xintercept = c(logFC_cutoff), linetype = 5) +
    scale_color_manual(values = structure(as.character(DMEs_plot$Var2), names = DMEs_plot$module %>% as.character())) +
    scale_fill_manual(values = structure(as.character(DMEs_plot$Var2), names = DMEs_plot$module %>% as.character())) +
    labs(x = xlabel, y = NULL, title = 'Differential module eigengene analysis') +
    theme_bw(base_rect_size = 1.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.text = element_text(size = 10, colour = 'black'),
          legend.title = element_text(face = 'bold', size = 12),
          axis.text.x = element_text(size = 12, colour = 'black'),
          axis.text.y = element_text(size = 12, colour = 'black'),
          axis.title = element_text(size = 14, colour = 'black'),
          panel.grid = element_blank()) +
    Seurat::NoLegend()
  ggsave(outfile, p, width = 5, height = 4)
}

# NOTE: recomputes FindDMEs to restore objects saved only in Part 2 workflow;
# rerun Part 2.8 first if DME results are missing.
group1 <- sce.obj@meta.data %>% subset(Branch1_cells == 'TRUE') %>% rownames
group2 <- sce.obj@meta.data %>% subset(Branch1_cells != 'TRUE') %>% rownames
DMEs.branch1 <- run_dme(group1, group2, 'branch')
plot_dme(DMEs.branch1, 'Average log2FC for Trajectory 1',
         file.path(dir_figures, '11_hdWGCNA_Branch1_DME.pdf'))

group1 <- sce.obj@meta.data %>% subset(Branch2_cells == 'TRUE') %>% rownames
group2 <- sce.obj@meta.data %>% subset(Branch2_cells != 'TRUE') %>% rownames
DMEs.branch2 <- run_dme(group1, group2, 'branch')
plot_dme(DMEs.branch2, 'Average log2FC for Trajectory 2',
         file.path(dir_figures, '11_hdWGCNA_Branch2_DME.pdf'))

# ---- 3.12 Group vs subtype DME scatter ----------------------------------------
DMEs <- qread(file.path(dir_hd, '11_hdWGCNA.DMEs.qs'))
p <- ggplot(data = DMEs, aes(x = avg_log2FC.group, y = avg_log2FC.celltype)) +
  geom_point(aes(color = module, fill = module, size = Nr.genes)) +
  ggrepel::geom_text_repel(aes(label = module, color = module), data = DMEs, size = 3) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = 5) +
  geom_hline(yintercept = c(-1.5, 1.5), linetype = 5) +
  scale_color_manual(values = structure(as.character(DMEs$Var2), names = DMEs$module)) +
  labs(x = 'Average log2 (Fold Change) for group',
       y = 'Average log2 (Fold Change) for GC_0',
       title = 'Differential module eigengene analysis') +
  scale_size_continuous(breaks = seq(100, 500, 100)) +
  guides(color = 'none', fill = 'none') +
  annotate("rect", xmin = -Inf, xmax = -0.25, ymin = 1.5, ymax = Inf, alpha = .2, fill = '#1e76b3') +
  annotate("rect", xmin = 0.25, xmax = Inf, ymin = 1.5, ymax = Inf, alpha = .2, fill = '#ff7d0b') +
  annotate("rect", xmin = -Inf, xmax = -0.25, ymin = -Inf, ymax = -1.5, alpha = .2, fill = '#1e76b3') +
  annotate("rect", xmin = 0.25, xmax = Inf, ymin = -Inf, ymax = -1.5, alpha = .2, fill = '#ff7d0b') +
  theme_bw(base_rect_size = 1.5) +
  theme(plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 10, colour = 'black'),
        legend.title = element_text(face = 'bold', size = 14),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.title = element_text(size = 14, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'right', panel.background = element_rect(fill = 'white'))
ggsave(file.path(dir_figures, '11_hdWGCNA_DME_scatter.pdf'), p, width = 5, height = 4)

# ---- 3.13 Pseudotime hME line plots + heatmaps (trajectory 1 and 2) -----------
pt_obj <- qread(file.path(dir_hd, '11_hdWGCNA_pseudotime_hMEs.qs'))
mod_colors <- pt_obj$mod_colors

plot_pseudotime_lines <- function(plot_df, outfile, width, height) {
  group_cp <- group_color
  names(group_cp) <- c('PCOS_pseudotime', 'Control_pseudotime')
  line_size <- 1

  selected_mod_colors <- mod_colors[levels(plot_df$variable)]
  class_strip <- ggh4x::strip_themed(background_x = ggh4x::elem_list_rect(fill = selected_mod_colors))

  plot_df$group <- factor(plot_df$group, levels = names(group_cp))
  p <- ggplot(plot_df, aes(x = as.numeric(bin), y = value, color = group)) +
    geom_smooth(se = FALSE, size = line_size, method = 'gam') +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    xlab("Pseudotime") + ylab("hMEs") +
    theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.line.x = element_blank(), axis.line.y = element_blank(),
          axis.text.y = element_text(size = 12), axis.title = element_text(size = 14),
          legend.text = element_text(size = 10, colour = 'black'),
          legend.title = element_text(face = 'bold', size = 14),
          panel.border = element_rect(size = 1, fill = NA, color = "black")) +
    labs(color = "Trajectory") +
    scale_color_manual(values = group_cp, labels = gsub('_pseudotime', '', names(group_cp)))
  p + ggh4x::facet_wrap2(. ~ variable, ncol = 4, nrow = 2, scales = "free", strip = class_strip) +
    theme(panel.background = element_blank(), strip.background = element_blank(),
          strip.placement = 'outside', panel.grid = element_blank(),
          strip.text.x.top = element_text(face = 'bold', color = 'white', size = 12))
}

plot_pseudotime_lines(pt_obj$plot_df_b1,
                      file.path(dir_figures, '11_hdWGCNA_pseudotime_hME_trajectory1_lineplot.pdf'),
                      4 * 2, 3.2)
ggsave(file.path(dir_figures, '11_hdWGCNA_pseudotime_hME_trajectory1_lineplot.pdf'),
       last_plot(), width = 4 * 2, height = 3.2, bg = 'white')

plot_pseudotime_lines(pt_obj$plot_df_b2,
                      file.path(dir_figures, '11_hdWGCNA_pseudotime_hME_trajectory2_lineplot.pdf'),
                      3 * 4, 2.6 * 2)
ggsave(file.path(dir_figures, '11_hdWGCNA_pseudotime_hME_trajectory2_lineplot.pdf'),
       last_plot(), width = 3 * 4, height = 2.6 * 2, bg = 'white')

cols_use <- c('#1f77b4', '#ff7f0e', '#2c9e6b', '#d62728', '#aa42fc', '#8c564b',
              '#ff7d0b', '#1e76b3')
vect_cols <- structure(cols_use, names = c(paste0('GC_', 0:5), 'PCOS', 'Control'))

plot_branch_heatmap <- function(meta_branch, modules_keep, height_cm, outfile) {
  plot_dat <- meta_branch %>% rownames_to_column('Cell_name') %>%
    arrange(Palantir_pseudotime) %>%
    mutate(order = 1:nrow(.))

  module_df <- plot_dat %>% dplyr::select(any_of(modules_keep)) %>% scale() %>% t()
  annotation_col <- plot_dat %>%
    dplyr::select(Palantir_pseudotime, dpt_pseudotime, monocle3_pseudotime,
                  Gr_subannotation, group)
  top_annotation <- HeatmapAnnotation(df = annotation_col, show_legend = TRUE,
                                      which = 'column', border = FALSE,
                                      annotation_label = c('Palantir pseudotime', 'DPT pseudotime',
                                                           'Monocle3 pseudotime', 'Cell Type', 'Group'),
                                      col = list(Gr_subannotation = subtype_color,
                                                 group = group_color),
                                      annotation_legend_param = list(grid_height = unit(3, "mm"),
                                                                     grid_width = unit(3, "mm")),
                                      show_annotation_name = TRUE)
  set.seed(1234)
  ht_my <- Heatmap(module_df, name = 'hMEs', height = unit(height_cm, 'cm'),
                   cluster_rows = FALSE, cluster_columns = FALSE,
                   show_column_names = FALSE, show_row_names = TRUE,
                   row_names_gp = gpar(fontsize = 10),
                   top_annotation = top_annotation,
                   heatmap_legend_param = list(direction = 'vertical',
                                               title_position = 'topcenter',
                                               grid_height = unit(3, 'mm'),
                                               grid_width = unit(3, 'mm')))
  pdf(outfile, width = 12, height = 6)
  draw(ht_my)
  dev.off()
}

plot_branch_heatmap(pt_obj$b1_pseudotime, paste0('M', c(8, 9)), 0.8 * 2,
                    file.path(dir_figures, '11_hdWGCNA_Module_pseudotime1_ht.pdf'))
plot_branch_heatmap(pt_obj$b2_pseudotime, paste0('M', c(9, 1, 5, 7, 3, 10, 2)), 0.8 * 7,
                    file.path(dir_figures, '11_hdWGCNA_Module_pseudotime2_ht.pdf'))

# ---- 3.14 Rora expression along pseudotime ------------------------------------
rora_obj <- qread(file.path(dir_hd, '11_hdWGCNA_Rora_pseudotime.qs'))
plot_df_rora <- rora_obj$plot_df_rora
b_pseudotime <- rora_obj$b_pseudotime

group_cp <- group_color
names(group_cp) <- c('PCOS_pseudotime', 'Control_pseudotime')
plot_df_rora$group <- factor(plot_df_rora$group, levels = names(group_cp))
p <- ggplot(plot_df_rora, aes(x = as.numeric(bin), y = value, color = group)) +
  geom_smooth(se = FALSE, size = 1) +
  xlab("Pseudotime") + ylab("Expression level of Rora") +
  theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 10, colour = 'black'),
        legend.title = element_text(size = 12),
        legend.background = element_blank(),
        panel.border = element_rect(size = 1, fill = NA, color = "black")) +
  scale_color_manual(values = group_cp, labels = gsub('_pseudotime', '', names(group_cp)))
patch <- p + theme(panel.background = element_blank(), strip.background = element_blank(),
                   plot.title = element_text(size = 16, hjust = .5),
                   legend.position = c(0.19, 0.88), legend.key.size = unit(.3, 'cm'))

plot_dat2 <- b_pseudotime %>% rownames_to_column('Cell_name') %>%
  arrange(Palantir_pseudotime) %>%
  mutate(order = 1:nrow(.))
ppp2 <- ggplot() +
  geom_tile(aes(x = order, y = 'Palantir_pseudotime', fill = Palantir_pseudotime),
            linewidth = 0.6, height = 0.9, data = plot_dat2, show.legend = FALSE) +
  geom_tile(aes(x = order, y = 'Group', color = group), data = plot_dat2,
            linewidth = 0.6, height = 0.9, show.legend = FALSE) +
  geom_tile(aes(x = order, y = 'Cell Type', color = Gr_subannotation),
            linewidth = 0.6, height = 0.9, data = plot_dat2, show.legend = FALSE) +
  scale_fill_gradientn(colors = viridis::viridis(100)) +
  scale_x_continuous(name = NULL, expand = c(0, 0), position = 'top') +
  scale_y_discrete(name = NULL, expand = c(0, 0), position = 'left') +
  scale_color_manual(values = vect_cols) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.title = element_text(size = 16))

p_rora <- ppp2 + patch + plot_layout(nrow = 2, heights = c(0.8, 3))
ggsave(file.path(dir_figures, '11_hdWGCNA_pseudotime_Rora_expr_lineplot.pdf'),
       p_rora, width = 4.5, height = 4, bg = 'white')

