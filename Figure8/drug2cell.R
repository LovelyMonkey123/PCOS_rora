# drug2cell drug activity prediction
rm(list = ls())
options(stringsAsFactors = FALSE)

library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)

# ---- Cell annotation for drug2cell metacells ----------------------------------

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

rownames(meta_use) <- meta_use$Row.names
identical(meta_use$Row.names, seurat_dat@meta.data %>% rownames(.))
# merge() reorders rows; restore the original cell order
meta_use <- meta_use[seurat_dat@meta.data %>% rownames(.), ]
identical(meta_use$Row.names, seurat_dat@meta.data %>% rownames(.))

seurat_dat@meta.data$Cluster_all <- meta_use$Cluster_all
# Remove GGGG cells
seurat_dat <- seurat_dat %>% subset(Cluster_all != 'GGGG')
seurat_dat$Cluster_all <- factor(seurat_dat$Cluster_all)

metadata <- seurat_dat@meta.data %>% dplyr::select(Cluster_all, group, Major_annotation_d)

# Export cell annotation
data.table::fwrite(metadata, './14_drug/14_metacell_cluster_all.csv', quote = FALSE, row.names = TRUE)

# ---- drug2cell: human-to-mouse drug target conversion -------------------------

d2c <- import('drug2cell')

# ChEMBL drug target data
a <- d2c$data$chembl()

source('./source/sc_functions.R')

if (TRUE) {
  bbb <- list()
  for (category_num in 65:100) {
    drug_list <- a[[category_num]]
    length_list <- seq_along(drug_list)
    names_use <- names(drug_list)
    drug_ls_res <- list()
    for (i in length_list) {
      drug_ls_res[[i]] <- convertHumanGeneList(drug_list[[i]])
    }
    names(drug_ls_res) <- names_use
    bbb[[names(a)[category_num]]] <- drug_ls_res
  }

  qsave(bbb, file = './12_drug2cell/drug2cell_drug_files_2.qs')

  drug_files <- qread('./12_drug2cell/drug2cell_drug_files.qs')
  drug_files <- c(drug_files, bbb)

  drug_files_ls <- drug_files %>% Reduce(c, .)
  # Drop empty target sets
  drug_files_ls_f <- drug_files_ls[map(drug_files_ls, length) %>% Reduce(c, .) != 0]

  qsave(drug_files_ls_f, file = './12_drug2cell/drug2cell_mouse_ls.qs')
  jsonlite::write_json(drug_files_ls_f, "./12_drug2cell/drug2cell_mouse_ls.json")
}

# ---- Drug DEG results for RORA(+) cells ---------------------------------------

# Per-subtype drug differential results (RORApos_GCs vs each GC subtype)
rora_meta_ls <- map(1:6, function(x) {
  rora_meta <- readxl::read_excel('./14_drug/14_drug_dea_RORA_group_results.xlsx', sheet = x) %>%
    as.data.frame() %>%
    mutate(name_a = gsub('\\.[0-9]{1,2}$', '', .$names)) %>%
    distinct(., name_a, .keep_all = TRUE)
  return(rora_meta)
})
names(rora_meta_ls) <- paste0('GC_', 0:5)

# Filtering cutoffs
pct_group <- 0.75
pct_ref <- 0.75
log_fc <- log2(1.5)
adj_p_cutoff <- 0.05

rora_drugs_deg <- map(names(rora_meta_ls), function(a) {
  b <- rora_meta_ls[[a]] %>%
    filter(pct_nz_group >= pct_group, logfoldchanges > log_fc,
           pct_nz_reference < pct_ref, pvals_adj < adj_p_cutoff) %>%
    mutate(Cluster_compare = paste0('Rora_pos_GCs_vs_', a),
           diff_pct = pct_nz_group - pct_nz_reference,
           drug_num = nrow(.)) %>%
    rename(avg_log2FC = logfoldchanges)
})

# Number of significant drugs per subtype: 171 66 286 163 420 133
map(rora_drugs_deg, nrow) %>% Reduce(c, .)

plot_dat <- rora_drugs_deg %>%
  Reduce(rbind, .)

paste0('Rora_pos_GCs_vs_', paste0('GC_', 0:5))

# Scatter of drug logFC (y) vs diff_pct (x)
p <- ggplot(plot_dat) +
  ggrastr::rasterise(geom_point(aes(x = diff_pct, y = avg_log2FC, color = Cluster_compare)),
                     dpi = 500, scale = 0.75) +
  geom_hline(yintercept = log_fc, linetype = 2, linewidth = 0.8, color = 'black') +
  scale_color_manual(values = c('Rora_pos_GCs_vs_GC_0' = '#1f77b4',
                                'Rora_pos_GCs_vs_GC_1' = '#ff7f0e',
                                'Rora_pos_GCs_vs_GC_2' = '#279e68',
                                'Rora_pos_GCs_vs_GC_3' = '#d62728',
                                'Rora_pos_GCs_vs_GC_4' = '#aa40fc',
                                'Rora_pos_GCs_vs_GC_5' = '#8c564b')) +
  labs(y = 'log2(Fold change)', x = NULL) +
  facet_wrap(~ Cluster_compare, nrow = 2) +
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

ggsave('./figures/14_drug2cell_DEG_Rora_point_plot.pdf', p, width = 4 * 2.8, height = 8)

# ---- Drugs shared across all subtype comparisons (upset plot) -----------------

common_drug_names <- rora_drugs_deg %>% map(., function(x) {
  x %>% dplyr::select(name_a, Cluster_compare)
})

names_list <- lapply(common_drug_names, `[[`, "name_a")
final_intersection <- Reduce(intersect, names_list)

library(ggupset)

dat_upset <- Reduce(rbind, common_drug_names)
dat_upset$Cluster_compare <- factor(dat_upset$Cluster_compare,
                                    levels = paste0('Rora_pos_GCs_vs_', paste0('GC_', 0:5)))
table(dat_upset$Cluster_compare)

dd <- dat_upset %>%
  group_by(name_a) %>%
  summarize(Cluster_compare = list(Cluster_compare))

p2_8 <- ggplot(data = dd, aes(x = Cluster_compare)) +
  geom_bar(fill = alpha('#1f77b4', .6), color = '#1f77b4', alpha = 0.66, linetype = 1) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -1) +
  scale_x_upset(order_by = 'degree') +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 76)) +
  theme_classic() +
  theme_combmatrix(combmatrix.panel.point.color.fill = c('#1f77b4'),
                   combmatrix.panel.line.color = '#00000f',
                   combmatrix.panel.striped_background = TRUE,
                   combmatrix.panel.line.size = 0.8,
                   combmatrix.label.text = element_text(color = c("#43A2CA",
                                                                  '#ff7f0e', '#279e68', '#d62728',
                                                                  '#aa40fc', '#8c564b')),
                   combmatrix.panel.point.size = 5,
                   combmatrix.panel.point.color.empty = 'grey')

p2_8

ggview:::ggview(p2_8, width = 6, height = 6)

ggsave('./figures/14_drug2cell_DEG_Rora_upset_last.pdf', p2_8, width = 6, height = 6)

# ---- ATC classification of the shared drugs -----------------------------------

final_df_ls <- rora_drugs_deg %>% map(., function(x) {
  x %>% filter(name_a %in% final_intersection)
})

# Drug category correspondence table
drug_corresponding_dat <- readxl::read_excel('./14_drug/chembl_drugs_corresponding.xlsx', sheet = 1) %>%
  as.data.frame()
colnames(drug_corresponding_dat) <- gsub(' ', '_', colnames(drug_corresponding_dat))

is.element(final_intersection, drug_corresponding_dat$Drug_Name) %>% table()

# ATC level 1 code dictionary
# Ref: https://www.atccode.com/
atc_level1 <- data.frame(
  code = c("A", "B", "C", "D", "G", "H", "J", "L", "M", "N", "P", "R", "S", "V"),
  name = c(
    "Alimentary tract and metabolism",
    "Blood and blood forming organs",
    "Cardiovascular system",
    "Dermatological drugs",
    "Genitourinary system and reproductive hormones",
    "Systemic hormonal preparations, excluding reproductive hormones and insulins",
    "Antiinfectives for systemic use",
    "Antineoplastic and immunomodulating agents",
    "Musculoskeletal system",
    "Nervous system",
    "Antiparasitic products, insecticides and repellents",
    "Respiratory system",
    "Sensory organs",
    "Various ATC structures"),
  stringsAsFactors = FALSE)

qsave(atc_level1, './14_drug/ATC_code_level1.qs')

# ATC level 2 codes present in the drug list
level2_data <- data.frame(code = c('G03', 'C10', 'L03', 'A09'),
                          level_2_name = c('Sex hormones and modulators of the genital system',
                                           'Lipid modifying agents',
                                           'Immunostimulants drugs',
                                           'Digestives, including enzymes'))

drug_corresponding_dat_level2 <- drug_corresponding_dat %>%
  filter(Drug_Name %in% final_intersection) %>%
  rename(code = ATC_Category) %>%
  .[c(20:nrow(.)), ] %>%
  merge(., level2_data, by = 'code') %>%
  dplyr::select(Drug_Name1, level_2_name)

drug_corresponding_dat_f <- drug_corresponding_dat %>%
  filter(Drug_Name %in% final_intersection) %>%
  rename(code = ATC_Category) %>%
  head(19) %>%
  merge(., atc_level1, by = 'code', all.x = TRUE) %>%
  merge(., drug_corresponding_dat_level2, by = 'Drug_Name1', all.x = TRUE)

drug_corresponding_dat_f[is.na(drug_corresponding_dat_f)] <- 'No-category'

drug_corresponding_dat_f_wr <- drug_corresponding_dat_f %>%
  dplyr::select(1, 2, 3, 4, 5, 6)
colnames(drug_corresponding_dat_f_wr) <- c('Drug_name', 'ATC_Code', 'full_Drug_name',
                                           'ChEMBL_ID', 'ATC_level1', 'ATC_level2')
drug_corresponding_dat_f_wr <- drug_corresponding_dat_f_wr %>%
  arrange(ATC_Code)

# Save final drug table
data.table::fwrite(drug_corresponding_dat_f_wr,
                   file = './14_drug/14_drug2cell_filtered_drug_for_RORA_final.csv',
                   quote = FALSE, row.names = TRUE)

# ---- Drug activity matrix -> Seurat object for plotting -----------------------

# Python environment (used upstream to extract the drug matrix)
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "celloracle_env")
reticulate::use_condaenv("celloracle_env", required = TRUE)

cat(paste0(drug_corresponding_dat_f_wr[, 3], '\n'))

# Drug matrix extracted in python; rows = cells, columns = drugs
mtx_drug <- data.table::fread('./14_drug/14_drug2cell_filtered_drug_for_rora_matrix.csv.gz', data.table = FALSE)
rownames(mtx_drug) <- mtx_drug[, 1]
mtx_drug <- mtx_drug[, -1]

mtx_drug_t <- mtx_drug %>% t()
mtx_drug_t <- mtx_drug_t[match(drug_corresponding_dat_f_wr[, 3], rownames(mtx_drug_t)), ]
rownames(mtx_drug_t) <- drug_corresponding_dat_f_wr[, 1]

# Cell annotation including RORA(+) label
meta_with_rora <- data.table::fread('./14_drug/14_drug2cell_Cluster_all_with_RORA.csv',
                                    data.table = FALSE, header = TRUE)
rownames(meta_with_rora) <- meta_with_rora[, 1]
meta_with_rora <- meta_with_rora[, -1, drop = FALSE]

# UMAP coordinates
umap_coords <- data.table::fread('./14_drug/14_drug2cell_Cluster_all_UMAP_coords.csv', data.table = FALSE)
rownames(umap_coords) <- umap_coords[, 1]
umap_coords <- umap_coords[, -1, drop = FALSE]

# Seurat object mainly for DotPlot-based visualization
library(Seurat)

seu_drug <- CreateSeuratObject(counts = mtx_drug_t,
                               meta.data = meta_with_rora,
                               assay = 'RNA')
seu_drug[['UMAP']] <- CreateDimReducObject(embeddings = umap_coords %>% as.matrix(),
                                           key = 'UMAP_',
                                           assay = 'RNA')

DimPlot(seu_drug, group.by = 'Cluster_all_with_RORA',
        cols = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
                 '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
                 '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d'))

Idents(seu_drug) <- seu_drug$Cluster_all_with_RORA

level_use <- c('RORApos_GCs', paste0('GC_', 0:5), 'Oocyte', 'Theca',
               'Epithelial', 'Endothelial', 'Perivascular', 'Stromal',
               'Macrophage', 'NK T cell', 'B', 'Neutrophil')
seu_drug$Cluster_all_with_RORA <- factor(seu_drug$Cluster_all_with_RORA, levels = level_use)

# Final drug activity dot plot
p_dot <- DotPlot(seu_drug,
                 features = rownames(seu_drug) %>% rev(),
                 group.by = 'Cluster_all_with_RORA',
                 cols = c("#ffffff", "#f781bf"),
                 cluster.idents = FALSE) +
  labs(y = 'Cell types', x = NULL, title = 'Drug activity score (drug2cell)') +
  scale_x_discrete(position = 'top') +
  coord_flip() +
  theme_bw() +
  RotatedAxis() +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

ggview:::ggview(width = 12, height = 8, bg = 'white')

drug_corresponding_dat_f_wr$Drug_name <- factor(drug_corresponding_dat_f_wr$Drug_name,
                                                  levels = drug_corresponding_dat_f_wr$Drug_name %>% rev())

# ATC level 2 labels flanking the dot plot
library(ggplot2)
library(cowplot)

p_text <- ggplot(data = drug_corresponding_dat_f_wr) +
  geom_text(aes(x = 1, y = Drug_name, label = ATC_level2),
            hjust = 1, size = 5) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_nothing() +
  theme(plot.margin = margin(0, 0, 0, 0))

ggview:::ggview(width = 12, height = 10, bg = 'white')
ggview:::ggview(p_text + p_dot, width = 17, height = 8, bg = 'white')

ggsave('./figures/14_drug2cell_drug_Rora_final_dotplot.pdf', p_text + p_dot,
       width = 17, height = 8, bg = 'white')

# Violin plots of drug scores by group
p_vln <- VlnPlot(seu_drug, features = rownames(seu_drug)[c(1:5, 18:19)],
                 group.by = 'group', pt.size = 0.00001, raster = TRUE)

ggview:::ggview(width = 12, height = 12, bg = 'white')

# Feature plots of selected drugs
p_ls <- FeaturePlot(seu_drug, features = rownames(seu_drug)[c(2:5, 18:19)], combine = FALSE)
p_ls <- p_ls %>% map(., function(x) {
  x + scale_color_viridis_c() +
    tidydr::theme_dr() +
    theme(panel.grid = element_blank(),
          legend.key.size = unit(.4, 'cm'))
}) %>% patchwork::wrap_plots(., ncol = 6)

p_ls %>% ggview:::ggview(width = 3 * 6, height = 3)

ggsave('./figures/14_drug2cell_drug_Rora_final_featureplot.pdf', p_ls,
       width = 3 * 6, height = 3, bg = 'white')
ggsave('./figures/14_drug2cell_drug_Rora_final_featureplot.png', p_ls,
       width = 3 * 6, height = 3, bg = 'white', dpi = 300, limitsize = FALSE)
