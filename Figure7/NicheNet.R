# NicheNet: RORA(+) GC -> Macrophage ligand-receptor analysis in PCOS

library(qs)
library(tidyverse)
library(data.table)
library(ggvenn)
library(nichenetr)
library(Seurat)
library(SeuratObject)
library(ggview)
library(scales)
library(patchwork)
library(cowplot)

dir_comm <- './17_cell_communication'
dir_figures <- './figures'

# ---- Select the macrophage target gene set -------------------------------------

# 1. PCOS-up markers of the macrophage cluster (single-cell DEG)
deg_res_ls <- qread(file.path(dir_comm, '17_2_celltype_group_DEG_results.qs'))
macro_deg_res <- deg_res_ls[[5]]
macro_deg_res <- macro_deg_res %>%
  mutate(diff_pct = (pct.1 - pct.2) / 100) %>%
  filter(pct.1 > 0.25, p_val_adj < 0.05, avg_log2FC > 0.25, diff_pct > 0.1) %>%
  arrange(desc(avg_log2FC))

# Exclude mitochondrial, ribosomal, and hemoglobin genes
genes_all <- macro_deg_res[['gene']]
MT_genes <- genes_all[str_starts(genes_all, 'mt-')]
ribo_genes <- genes_all[grep("^(RP[SL]|Rps|Rpl)", genes_all)]
hb_genes <- genes_all[grep("^Hb[^(p)]", genes_all)]
macro_deg_res <- macro_deg_res %>%
  filter(!gene %in% c(MT_genes, ribo_genes, hb_genes)) %>%
  arrange(desc(avg_log2FC))

# 2. Macrophage markers from the major cell type annotation (pseudobulk DEG)
macro_cluster_genes <- readxl::read_xlsx('./07_major_celltype_DEG/07_Major_annotation_d_dea_results.xlsx',
                                         sheet = 5)
macro_cluster_genes_filter <- macro_cluster_genes %>%
  mutate(diff_pct = pct_nz_group - pct_nz_reference) %>%
  filter(pct_nz_group > 0.25, pvals_adj < 0.05, logfoldchanges > 0.25, diff_pct > 0.1) %>%
  filter(!names %in% c(MT_genes, ribo_genes, hb_genes)) %>%
  arrange(desc(logfoldchanges))

# Intersection of the two marker sets
genes_common <- intersect(macro_cluster_genes_filter$names, macro_deg_res$gene)

a <- macro_cluster_genes_filter %>% filter(names %in% genes_common) %>%
  dplyr::select(names, logfoldchanges, pvals_adj, diff_pct) %>%
  rename(gene = names,
         avg_log2FC_cluster = logfoldchanges,
         p_val_adj_cluster = pvals_adj,
         diff_pct_cluster = diff_pct)
b <- macro_deg_res %>% filter(gene %in% genes_common) %>%
  dplyr::select(gene, avg_log2FC, p_val_adj, diff_pct) %>%
  rename(gene = gene,
         avg_log2FC_PCOS = avg_log2FC,
         p_val_adj_PCOS = p_val_adj,
         diff_pct_PCOS = diff_pct)

gene_common_df <- a %>% merge(., b, by = 'gene') %>%
  arrange(desc(avg_log2FC_cluster), desc(avg_log2FC_PCOS))

gene_common_df %>% data.table::fwrite(., file.path(dir_comm, '17_Macrophage_cluster_group_DEG_df.csv'))

# ---- Fig. A: venn of the two macrophage marker sets ------------------------------

venn_data <- list(
  Macrophage_cluster_markers = macro_cluster_genes_filter$names,
  Macrophage_PCOS_markers = macro_deg_res$gene)

p7 <- ggvenn(venn_data,
             c("Macrophage_cluster_markers", "Macrophage_PCOS_markers"),
             text_size = 6, auto_scale = FALSE,
             set_name_size = 5, show_elements = FALSE,
             show_percentage = TRUE, stroke_size = .4, stroke_linetype = 1,
             stroke_color = "black",
             fill_color = c("#b057fd", "#f7a454"),
             set_name_color = c("#b057fd", "#f7a454"))
p7


ggsave(file.path(dir_figures, '17_Macrophage_common_genes_vennplot.pdf'),
       p7, width = 7, height = 7)

# ---- NicheNet --------------------------------------------------------------------

options(timeout = 600)
organism <- "mouse"

# Prior models
lr_network <- read_rds('./source/NicheNet/lr_network_mouse_21122021.rds')
ligand_target_matrix <- read_rds('./source/NicheNet/ligand_target_matrix_nsga2r_final_mouse.rds')
weighted_networks <- readRDS('./source/NicheNet/weighted_networks_nsga2r_final_mouse.rds')

# Deduplicate
lr_network <- lr_network %>% distinct(from, to)

# Seurat object (RORA(+) labels already merged)
sce_dat <- qread(file.path(dir_comm, 'RORA_sce_dat_PCOS_for_cellchat.qs'))
seuratObj <- alias_to_symbol_seurat(sce_dat, "mouse")  # map names to symbols

Idents(seuratObj) <- seuratObj$Cluster_PCOS

receiver <- "Macrophage"          # receiving cell type
sender_celltypes <- c("RORApos_GCs")  # sending cell type

# Subset to sender and receiver
seuratObj <- seuratObj %>% subset(Cluster_PCOS %in% c(receiver, sender_celltypes))

# a. Filter expressed receptors in the receiver
expressed_genes_receiver <- get_expressed_genes(receiver, seuratObj, pct = 0.1)
length(expressed_genes_receiver)
all_receptors <- unique(lr_network$to)
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)
expressed_receptors %>% length()

# b. Filter ligands matching the expressed receptors (step 1: prior knowledge)
potential_ligands <- lr_network %>%
  filter(to %in% expressed_receptors) %>%
  pull(from) %>%
  unique()

# step 2: keep ligands expressed in the sender
list_expressed_genes_sender <- sender_celltypes %>% unique() %>%
  lapply(get_expressed_genes, seuratObj, 0.1)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()
potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender)

potential_ligands_focused %>% length()
expressed_receptors %>% length()

# Weighted ligand-receptor network
weighted_networks_lr <- weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from", "to"))

# Background gene set
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

qsave(background_expressed_genes, file = file.path(dir_comm, '17_1_nichenet_backgroundgenes.qs'))

# Target geneset: PCOS-up macrophage markers (intersection set above)
gene_common_df <- data.table::fread(file.path(dir_comm, '17_Macrophage_cluster_group_DEG_df.csv'),
                                    data.table = FALSE)
geneset_oi <- gene_common_df$gene
geneset_oi %>% length()

# step 3: ligand activity prediction
nichenet_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  potential_ligands = potential_ligands_focused,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix)

nichenet_activities <- nichenet_activities %>%
  dplyr::arrange(-aupr_corrected) %>%
  mutate(rank = rank(dplyr::desc(aupr_corrected)))
nichenet_activities

qsave(nichenet_activities, file = file.path(dir_comm, '17_1_nichenet_ligand_activities.qs'))

nichenet_activities <- qread(file.path(dir_comm, '17_1_nichenet_ligand_activities.qs'))

# ---- Fig. B: prioritize ligands also upregulated in RORA(+) cells ----------------

# DEGs of RORA(+) vs other cells
res_deg <- qread('./12_SCENIC/12_DEA_for_RORA_cells.qs')
res_deg_filter <- res_deg %>%
  mutate(Ratio = round(pct_in / pct_out, 3),
         pct.fc = pct_in - pct_out) %>%
  filter(padj < 0.05, pct.fc > 25, logFC > 0.25) %>%
  filter(group == 1)

# Intersection of RORA(+) cell markers and active ligands
ligand_common <- res_deg_filter$feature %>% intersect(., nichenet_activities$test_ligand)
nichenet_activities$test_ligand %>% length()

# Four-way venn
genes_A <- macro_cluster_genes_filter$names          # macrophage cluster markers
genes_B <- macro_deg_res$gene                        # macrophage PCOS-up markers
genes_C <- nichenet_activities$test_ligand           # inferred ligands
genes_D <- res_deg_filter$feature                    # RORA(+) cell markers

venn_data <- list(
  Macrophage_cluster_markers = genes_A,
  Macrophage_PCOS_markers = genes_B,
  NicheNet_potential_ligands = genes_C,
  `RORA(+)_pos_cells_markers` = genes_D)

p7 <- ggvenn(venn_data,
             c(names(venn_data)), text_size = 6, auto_scale = FALSE,
             set_name_size = 5, show_elements = FALSE,
             show_percentage = TRUE, stroke_size = .4, stroke_linetype = 1,
             stroke_color = "black",
             fill_color = c("#8e3119", "#00649d", "#50b688", "#ffe071"),
             set_name_color = c("#8e3119", "#00649d", "#50b688", "#ffe071")) +
  theme(plot.margin = margin(r = 0.4, l = 0.4, unit = 'cm'),
        panel.background = element_blank())
p7


ggsave(file.path(dir_figures, '17_NicheNet_RORA_ligands_genes_vennplot.pdf'),
       p7, width = 7, height = 7)

# Ligands in the RORA marker set but not macrophage markers (final 16 ligands)
genes_aaa <- ligand_common %>% intersect(., macro_cluster_genes_filter$names)
genes_bbb <- ligand_common %>% intersect(., macro_deg_res$gene)
genes_ccc <- c(genes_aaa, genes_bbb) %>% unique()

cat(paste0(setdiff(ligand_common, genes_ccc), '\n'))

ligand_intersect <- setdiff(ligand_common, genes_ccc)

nichenet_activities <- nichenet_activities %>%
  filter(test_ligand %in% ligand_intersect) %>%
  dplyr::arrange(-aupr_corrected)

cat(paste0(nichenet_activities$test_ligand %>% as.character(), '\n'))

# ---- Fig. C left: ligand expression dot plot -------------------------------------

receiver <- "Macrophage"
sender_celltypes <- c("RORApos_GCs")

p_dotplot <- DotPlot(seuratObj,
                     assay = 'RNA',
                     group.by = 'Cluster_PCOS',
                     features = as.character(nichenet_activities$test_ligand) %>% rev(),
                     idents = c(sender_celltypes, receiver),
                     cols = "RdYlBu") +
  RotatedAxis() +
  coord_flip() +
  theme(axis.title.x = element_blank(),
        legend.position = 'left')

p_dotplot

# ---- Active target gene inference -------------------------------------------------

best_upstream_ligands <- nichenet_activities$test_ligand %>% as.character()

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 250) %>% bind_rows() %>% drop_na()

# ---- Fig. D right: ligand-target heatmap ------------------------------------------

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>%
  rev() %>% make.names()
order_targets <- active_ligand_target_links_df$target %>% unique() %>%
  intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) <- rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) <- colnames(active_ligand_target_links) %>% make.names()

qsave(list(order_ligands = order_ligands, order_targets = order_targets),
      file = file.path(dir_comm, '17_1_nichenet_ligands_targets.qs'))

ligands_targets <- qread(file.path(dir_comm, '17_1_nichenet_ligands_targets.qs'))
order_targets <- ligands_targets$order_targets
order_ligands <- ligands_targets$order_ligands

vis_ligand_target <- active_ligand_target_links[order_targets %>% rev(), order_ligands] %>% t()

p_ligand_target_network <- vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands in Rora(+)_pos_cells",
                      "Predicted target genes in Macrophage cells",
                      legend_position = "right",
                      x_axis_position = "top",
                      legend_title = "Regulatory potential") +
  theme(axis.text = element_text(size = 12, hjust = .5, vjust = .5),
        legend.key.size = unit(.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 10)) +
  scale_fill_gradient2(low = c("whitesmoke"), high = '#3f77ab')
p_ligand_target_network


ggsave(file.path(dir_figures, '17_1_NicheNet_ligand_target_activity.pdf'),
       p_ligand_target_network, width = 5, height = 5, bg = 'white')

# ---- Fig. D left: ligand activity (AUPR) bars -------------------------------------

nichenet_activities <- nichenet_activities[match(order_ligands, nichenet_activities$test_ligand), ]
nichenet_activities$test_ligand <- factor(nichenet_activities$test_ligand,
                                          levels = nichenet_activities$test_ligand,
                                          ordered = TRUE)

p_ligand_aupr <- nichenet_activities %>% ggplot() +
  geom_bar(aes(y = test_ligand, x = aupr, fill = aupr), stat = "identity") +
  scale_fill_gradient(low = "whitesmoke", high = '#785ca7') +
  scale_x_continuous(expand = c(0, 0), position = 'top') +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = 'Ligand activity', y = NULL) +
  theme_classic() +
  theme(plot.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        legend.key.size = unit(.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.position = 'top',
        legend.background = element_blank(),
        plot.margin = margin(l = -2))

# ---- Fig. E: target gene expression dot plot --------------------------------------

p_dotplot2 <- DotPlot(seuratObj,
                      assay = 'RNA',
                      group.by = 'Cluster_PCOS',
                      features = as.character(order_targets),
                      idents = c(sender_celltypes, receiver),
                      cols = "RdYlBu") +
  RotatedAxis() +
  theme(axis.title.x = element_blank(),
        axis.title.y.left = element_blank())

p_dotplot2


ggsave(file.path(dir_figures, '17_1_NicheNet_target_genes_expression_dotplot.pdf'),
       p_dotplot2, width = 5.5, height = 5)

# ---- Ligand-receptor activity heatmap ---------------------------------------------

best_upstream_ligands <- nichenet_activities$test_ligand %>% as.character()

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands,
  expressed_receptors,
  lr_network, weighted_networks$lr_sig)

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both")

vis_ligand_receptor_network %>% colnames()
vis_ligand_receptor_network %>% rownames(.)

p_ligand_receptor <- (make_heatmap_ggplot(t(vis_ligand_receptor_network[, best_upstream_ligands]),
                                          y_name = "Prioritized Rora(+)_pos_cells-ligands",
                                          x_name = "Receptors expressed by Macrophage cells",
                                          color = '#D1806B',
                                          legend_title = "Prior interaction potential")) +
  theme(axis.text.x = element_text(hjust = .5, vjust = .5, size = 12),
        legend.key.size = unit(.3, 'cm'),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 10))

p_ligand_receptor

receiver <- "Macrophage"
sender_celltypes <- c("RORApos_GCs")
p_dotplot1 <- DotPlot(seuratObj,
                      assay = 'RNA',
                      group.by = 'Cluster_PCOS',
                      features = as.character(nichenet_activities$test_ligand) %>% rev(),
                      idents = c(sender_celltypes, receiver),
                      cols = "RdYlBu") +
  RotatedAxis() +
  coord_flip() +
  labs(x = 'Prioritized Rora(+)_pos_cells-ligands') +
  theme(axis.title.x = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

p_dotplot1

# ---- Combined panel ---------------------------------------------------------------

p_comb <- plot_grid(
  p_dotplot1 + theme(legend.position = "none"),
  p_ligand_receptor + theme(legend.position = "none",
                            axis.text.y = element_blank(),
                            axis.title.y = element_blank()),
  p_ligand_aupr + theme(legend.position = "none"),
  align = "hv",
  ncol = 3,
  rel_widths = c(1.5, 4, 1.5))


ggsave(file.path(dir_figures, '17_1_NicheNet_ligand_receptor_activity.pdf'),
       p_comb, width = 10, height = 6, bg = 'white')

# ---- Extracted legends ------------------------------------------------------------

legends <- plot_grid(
  as_ggplot(get_legend(p_dotplot1)),
  as_ggplot(get_legend(p_ligand_receptor)),
  as_ggplot(get_legend(p_ligand_aupr)),
  nrow = 2,
  align = "v")

plot_grid(legends,
          rel_heights = c(10, 2), nrow = 2, align = "hv")

ggview:::ggview(width = 12, height = 10)

plot_grid(legends)

ggview:::ggview(width = 6, height = 6)

ggsave(file.path(dir_figures, '17_1_NicheNet_ligand_receptor_activity_legends.pdf'),
       plot_grid(legends), width = 6, height = 6, bg = 'white')
