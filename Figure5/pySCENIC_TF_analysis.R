# SCENIC workflow: metacell construction, pySCENIC input preparation, and
# regulon activity score (RAS) analysis.
# Part 1: packages | Part 2: data processing | Part 3: visualization
# Ref: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html

# ============================== Part 1. Packages ==============================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(qs)
  library(SCopeLoomR)
  library(data.table)
  library(arrow)
  library(clusterProfiler)
  library(viridis)
  library(ggsci)
  library(tidydr)
  library(cowplot)
})

source('./source/makeMetaCells.R')
source('./source/compute_module_score.R')

set.seed(123)
n_cores <- 10
dir_scenic <- './12_SCENIC'
dir_scenic_out <- './12_SCENIC/output'
dir_cistarget <- './cisTarget_db'
dir_figures <- './figures'


# ============================ Part 2. Data processing =========================

# ---- 2.1 Build metacells per sample (30000 cells -> ~1800 metacells) ---------
seu <- qread('./11_hdWGCNA/11_gr_subannotation_seurat_data.qs')

seu.list <- SplitObject(seu, split.by = "orig.ident")
metacells.list <- lapply(seq_along(seu.list), function(ii) {
  makeMetaCells(seu = seu.list[[ii]], min.cells = 10,
                reduction = "X_umap", dims = 1:2, k.param = 10, cores = n_cores)
})
mc.mat <- lapply(metacells.list, function(mc) mc$mat) %>% Reduce(cbind, .)
mc.cellmeta <- lapply(metacells.list, function(mc) mc$metadata) %>% Reduce(rbind, .)

seu2 <- CreateSeuratObject(mc.mat)
seu2 <- NormalizeData(seu2)

dir.create(dir_scenic_out, showWarnings = FALSE, recursive = TRUE)
qsave(mc.mat, file.path(dir_scenic_out, '00-1.mc.mat.qs'))

# ---- 2.2 Prepare pySCENIC input files ----------------------------------------
# (1) TF list from the cisTarget motif2tf table
motif2tfs <- data.table::fread(file.path(dir_cistarget, 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'))
TFs <- sort(unique(motif2tfs$gene_name))
dir.create(dir_cistarget, showWarnings = FALSE)
writeLines(TFs, file.path(dir_cistarget, 'mgi_tfs.motifs-v10.txt'))

# (2) Metacell expression matrix as loom
mc.mat <- qread(file.path(dir_scenic_out, '00-1.mc.mat.qs'))
# Drop lowly expressed genes (detected in < 5 metacells)
expr.in.cells <- rowSums(mc.mat > 0)
mc.mat <- mc.mat[expr.in.cells >= 5, ]
# Keep only genes present in the cisTarget ranking database
cisdb <- arrow::read_feather(file.path(dir_cistarget, 'mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'))
genes.use <- intersect(colnames(cisdb), rownames(mc.mat))
mc.mat <- mc.mat[genes.use, ]

loom <- SCopeLoomR::build_loom(
  file.name = file.path(dir_scenic_out, '00-2.mc_mat_for_step1.loom'),
  dgem = mc.mat,
  default.embedding = NULL)
loom$close()
rm(loom)
gc()

##### run pySCENIC in shell #####

# ---- 2.3 Parse pySCENIC regulons ---------------------------------------------
seu <- qread('./11_hdWGCNA/11_gr_subannotation_seurat_data.qs')
regulons <- clusterProfiler::read.gmt(file.path(dir_scenic_out, 'GC.regulons.gmt'))
rg.names <- unique(regulons$term)
regulon.list <- lapply(rg.names, function(rg) {
  subset(regulons, term == rg)$gene
})
names(regulon.list) <- sub("[0-9]+g", "\\+", rg.names)
qsave(regulon.list, file.path(dir_scenic_out, '03-1.GC.regulons.qs'))

# ---- 2.4 Regulon activity scores (RAS) via AUCell ----------------------------
regulon.list <- qread(file.path(dir_scenic_out, '03-1.GC.regulons.qs'))
seu <- ComputeModuleScore(seu, gene.sets = regulon.list, min.size = 10, cores = 6)
DefaultAssay(seu) <- "AUCell"

# UMAP and PCA on the RAS matrix (correlation metric works best here)
seu <- RunUMAP(object = seu, features = rownames(seu), metric = "correlation",
               reduction.name = "umapRAS", reduction.key = "umapRAS_")
DefaultAssay(seu) <- "AUCell"
seu <- ScaleData(seu)
seu <- RunPCA(object = seu, features = rownames(seu),
              reduction.name = "pcaRAS", reduction.key = "pcaRAS_")

# Cluster cells on the RAS PCA
seu <- seu %>%
  FindNeighbors(., reduction = 'pcaRAS') %>%
  FindClusters(., resolution = 1.0)

qsave(seu, file.path(dir_scenic_out, 'gc_annotation_seu.qs'))

# ---- 2.5 Export Rora regulon genes -------------------------------------------
regulon.list <- qread(file.path(dir_scenic_out, '03-1.GC.regulons.qs'))
dat_rora <- regulon.list[['Rora(+)']]

# 6-column padded matrix for external tool input
n_col <- 6
n_row <- ceiling(length(dat_rora) / n_col)
mat <- matrix(c(dat_rora, rep(NA, n_row * n_col - length(dat_rora))),
              nrow = n_row, ncol = n_col, byrow = TRUE)
df <- as.data.frame(mat)
df %>% data.table::fwrite(file.path(dir_scenic, '12_rora_regulon_genes.csv'),
                          row.names = FALSE, quote = FALSE, na = '')


# =========================== Part 3. Visualization ============================

seu <- qread(file.path(dir_scenic_out, 'gc_annotation_seu.qs'))

dimplot_theme <- theme(panel.grid = element_blank(),
                       panel.background = element_blank(),
                       plot.title = element_text(hjust = .5, size = 16),
                       legend.text = element_text(size = 14),
                       axis.title = element_text(size = 14))

# ---- 3.2 PCA on RAS ----------------------------------------------------------
# PC3 separates GC_3, PC4/PC6 separate GC_0
p3 <- DimPlot(seu, reduction = "pcaRAS", group.by = "Gr_subannotation",
              dims = c(3, 8), pt.size = 0.0001) + ggsci::scale_color_d3("category20")
p4 <- DimPlot(seu, reduction = "pcaRAS", group.by = "group",
              cols = rev(c("#da0612", '#2e92f7')), dims = c(3, 8),
              split.by = 'group', pt.size = 0.0001)
p_pca <- cowplot::plot_grid(plotlist = list(p3, p4), rel_widths = c(1.2, 2))
ggsave(file.path(dir_figures, '12_SCENIC_PCAras_plot.pdf'), p_pca,
       width = 12, height = 4, bg = 'white')

# ---- 3.3 UMAP on RAS ---------------------------------------------------------
p3 <- DimPlot(seu, reduction = "umapRAS", group.by = "Gr_subannotation",
              pt.size = 0.000001) + ggsci::scale_color_d3("category20") +
  tidydr::theme_dr() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  dimplot_theme
p4 <- DimPlot(seu, reduction = "umapRAS", group.by = "seurat_clusters",
              pt.size = 0.000001) + ggsci::scale_color_d3("category20") +
  tidydr::theme_dr() +
  guides(color = guide_legend(override.aes = list(size = 3), ncol = 2)) +
  dimplot_theme
p5 <- DimPlot(seu, reduction = "umapRAS", group.by = "group",
              cols = c("#629fca", '#ffa657'), pt.size = 0.000001) +
  tidydr::theme_dr() +
  guides(color = guide_legend(override.aes = list(size = 3))) +
  dimplot_theme

ggsave(file.path(dir_figures, '12_SCENIC_UMAPras_plot.pdf'), p3 + p4 + p5,
       width = 4.5 * 3, height = 4, bg = 'white')
