# Augur: cell type prioritization of perturbation response (GC subtypes)
# Ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC11798837 (Fig. 2)

rm(list = ls())
options(stringsAsFactors = FALSE)

library(reticulate)
library(Augur)
library(Matrix)
library(tidyverse)
library(data.table)
library(qs)
library(ggplot2)
library(ggforce)

# Python environment for reading the h5ad object
options(reticulate.conda_binary = "./miniconda3/bin/conda", sc_env_name = "sc_ss_ov")
reticulate::use_condaenv("sc_ss_ov", required = TRUE)

dir_augur <- './09_augur'
dir_figures <- './figures'

sc <- import('scanpy')

adata_gr <- sc$read_h5ad('./08_GC_subclustering/08_after_hb_Gr_subanndata.h5ad')

# Preprocessed expression data (raw.X; genes as columns in python, transpose)
mtx <- adata_gr$raw$X %>% t()

# Spot/gene metadata
meta <- adata_gr$obs
var <- adata_gr$var

# Add row/column names to the count matrix
colnames(mtx) <- rownames(meta)
rownames(mtx) <- rownames(var)

# Phenotype metadata
pheno.data <- meta

# ---- Augur: calculate AUC -----------------------------------------------------

augur_res <- Augur::calculate_auc(
  input = mtx %>% as.matrix(),
  meta = pheno.data,
  label_col = 'group',
  cell_type_col = 'Gr_subannotation',
  n_subsamples = 50,
  subsample_size = 200,
  feature_perc = 0.8,
  var_quantile = 0.9,
  n_threads = 8,
  show_progress = TRUE,
  classifier = 'rf',
  rf_params = list(trees = 100, mtry = 2, min_n = NULL, importance = "accuracy"))

augur_res_ls <- augur_res[-1]

qsave(augur_res_ls, file = file.path(dir_augur, '09_gr_augur_res.qs'))

# ---- Plot AUC per cell type ---------------------------------------------------

augur_res_ls <- qread(file.path(dir_augur, '09_gr_augur_res.qs'))

plot_data <- augur_res_ls$AUC %>%
  arrange(desc(auc)) %>%
  mutate(cell_type = factor(cell_type, levels = rev(cell_type)))

p <- ggplot(plot_data) +
  geom_link(aes(x = 0.5, y = cell_type,
                xend = auc, yend = cell_type,
                alpha = after_stat(index),
                color = cell_type,
                size = after_stat(index)),
            n = 500, show.legend = FALSE)

p1 <- p +
  geom_point(aes(x = auc, y = cell_type), color = "black", fill = "white",
             size = 6, shape = 21) +
  geom_text(aes(x = auc, y = cell_type), label = plot_data[['auc']] %>% round(., 3),
            size = 3, nudge_x = 0.05) +
  geom_vline(xintercept = 0.5, linetype = 2, color = 'grey70') +
  xlab("AUC") + ylab("") +
  labs(title = 'Augur') +
  scale_x_continuous(limits = c(0.5, 1)) +
  scale_color_manual(values = c(
    "GC_0" = "#1f77b4",
    "GC_1" = "#ff7f0e",
    "GC_2" = "#2ca02c",
    "GC_3" = "#d62728",
    "GC_4" = "#9467bd",
    "GC_5" = "#8c564b")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = .5, size = 12))

ggsave(file.path(dir_figures, '09_augur_plot.pdf'), p1,
       width = 3.2 * 1.2, height = 4 * 1.2)
