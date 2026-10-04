# Palantir pseudotime visualization and correlation between pseudotime methods

rm(list = ls())
options(stringsAsFactors = FALSE)

library(reticulate)
library(Matrix)
library(tidyverse)
library(data.table)
library(qs)
library(Seurat)
library(corrplot)
library(ggraph)
library(igraph)

dir_pt <- './10_cytotrace2'
dir_figures <- './figures'

# ---- Load pseudotime results and UMAP coordinates -----------------------------

pr.res <- read.table(file.path(dir_pt, 'pseudotime_res_all.csv'),
                     sep = ",", header = TRUE, row.names = 1)

umap_coord <- data.table::fread(file.path(dir_pt, '10_UMAP_coords.csv'),
                                sep = ",", header = TRUE, data.table = FALSE)
rownames(umap_coord) <- umap_coord[, 1]
umap_coord <- umap_coord[, -1]
umap_coord <- umap_coord %>% as.data.frame()

plot_data <- umap_coord %>% cbind(., pr.res)

# ---- Branch assignments on UMAP ------------------------------------------------

p_branch1 <- ggplot() +
  ggrastr::rasterise(geom_point(size = 0.01,
                                data = plot_data %>% filter(Branch1_cells == 'False'),
                                mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Branch1_cells')),
                     dpi = 500, scale = 1) +
  ggrastr::rasterise(geom_point(size = 0.5,
                                data = plot_data %>% filter(Branch1_cells == 'True'),
                                mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Branch1_cells')),
                     dpi = 500, scale = 1) +
  scale_color_manual(values = rev(c('#5784db', 'grey90'))) +
  labs(title = 'GC cells in Trajectory 1', x = 'UMAP 1', y = 'UMAP 2') +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = .5),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  Seurat::NoLegend()

p_branch2 <- ggplot() +
  ggrastr::rasterise(geom_point(size = 0.01,
                                data = plot_data %>% filter(Branch2_cells == 'False'),
                                mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Branch2_cells')),
                     dpi = 500, scale = 1) +
  ggrastr::rasterise(geom_point(size = 0.5,
                                data = plot_data %>% filter(Branch2_cells == 'True'),
                                mapping = aes_string(x = 'UMAP_1', y = 'UMAP_2', color = 'Branch2_cells')),
                     dpi = 500, scale = 1) +
  scale_color_manual(values = rev(c('#5784db', 'grey90'))) +
  labs(title = 'GC cells in Trajectory 2', x = 'UMAP 1', y = 'UMAP 2') +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = .5),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  Seurat::NoLegend()

# ---- Differential potential along pseudotime -----------------------------------

gc_colors <- c("GC_0" = "#1f77b4",
               "GC_1" = "#ff7f0e",
               "GC_2" = "#2ca02c",
               "GC_3" = "#d62728",
               "GC_4" = "#9467bd",
               "GC_5" = "#8c564b")

p_potential_path1 <- ggplot() +
  ggrastr::rasterise(geom_jitter(mapping = aes_string(x = 'Palantir_pseudotime', y = 'Branch1',
                                                      color = 'Gr_subannotation'),
                                 data = plot_data %>% filter(Branch1_cells == 'True'),
                                 size = .01), dpi = 500, scale = 0.75) +
  ggrastr::rasterise(geom_jitter(mapping = aes_string(x = 'Palantir_pseudotime', y = 'Branch1',
                                                      color = 'Branch1_cells'),
                                 data = plot_data %>% filter(Branch1_cells == 'False'),
                                 size = .01), dpi = 500, scale = 0.75) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4), reverse = TRUE)) +
  scale_color_manual(values = c('False' = 'grey90', gc_colors)) +
  labs(x = 'Pseudotime', y = 'Differential Potential', title = 'Trajectory 1') +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5, size = 14)) +
  Seurat::NoLegend()

p_potential_path2 <- ggplot() +
  ggrastr::rasterise(geom_jitter(mapping = aes_string(x = 'Palantir_pseudotime', y = 'Branch2',
                                                      color = 'Gr_subannotation'),
                                 data = plot_data %>% filter(Branch2_cells == 'True'),
                                 size = .01), dpi = 500, scale = 0.75) +
  ggrastr::rasterise(geom_jitter(mapping = aes_string(x = 'Palantir_pseudotime', y = 'Branch2',
                                                      color = 'Branch2_cells'),
                                 data = plot_data %>% filter(Branch2_cells == 'False'),
                                 size = .01), dpi = 500, scale = 0.75) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4), reverse = TRUE)) +
  scale_color_manual(values = c('False' = 'grey90', gc_colors)) +
  labs(x = 'Pseudotime', y = 'Differential Potential', title = 'Trajectory 2') +
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'right',
        plot.title = element_text(hjust = .5, size = 14))

ggsave(file.path(dir_figures, '10_Palantir_potential.pdf'),
       p_potential_path1 + p_potential_path2,
       width = 4.5 * 2, height = 4)
ggsave(file.path(dir_figures, '10_Palantir_trajectory.pdf'),
       p_branch1 + p_branch2,
       width = 3.8 * 2, height = 4)

# ---- Correlation between pseudotime methods ------------------------------------

dat <- data.table::fread(file.path(dir_pt, 'pseudotime_res_all.csv'), data.table = FALSE)

cor.dat <- dat[, c('dpt_pseudotime', 'monocle3_pseudotime', 'Palantir_pseudotime')]
plot_dat <- cor(cor.dat, method = 'spearman')

# Correlation network; keep only the upper triangle (exclude the diagonal)
plot_dat[lower.tri(plot_dat, diag = TRUE)] <- NA

plot_dat_df <- plot_dat %>% as.data.frame() %>%
  rownames_to_column('index_a') %>%
  pivot_longer(., !index_a, names_to = 'index_b', values_to = 'Cor') %>%
  na.omit()

graph <- graph_from_data_frame(plot_dat_df, directed = FALSE)
V(graph)$color <- names(V(graph))

set.seed(12346)

p_cor <- ggraph(graph, layout = 'fr') +
  geom_edge_link(aes(width = Cor, label = round(Cor, 2)),
                 angle_calc = 'along',
                 color = 'grey66',
                 show.legend = FALSE,
                 label_dodge = unit(.3, 'cm')) +
  geom_node_point(aes(color = color), size = 10, show.legend = FALSE) +
  geom_node_text(aes(label = name), vjust = 1, hjust = 1, check_overlap = TRUE) +
  scale_edge_width(name = 'Cor', range = c(0.6, 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_graph(base_family = "sans") +
  theme(legend.position = "right",
        plot.margin = margin(l = 2, unit = 'cm')) +
  labs(edge_width = "Correlation")

ggsave(file.path(dir_figures, '10_cor_pseudotime.pdf'), p_cor, width = 6, height = 6)
