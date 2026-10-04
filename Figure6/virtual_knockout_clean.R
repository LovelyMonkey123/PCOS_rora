# In-silico knockdown of Rora in granulosa cells
# scTenifoldKnk ref: https://github.com/cailab-tamu/scTenifoldKnk

rm(list = ls())
options(stringsAsFactors = FALSE)

library(scTenifoldKnk)
library(tidyverse)
library(qs)
library(Seurat)

# ---- scTenifoldKnk: virtual Rora knockdown (PCOS cells) -----------------------

# Load hdWGCNA object; drop the module columns appended to metadata
sce.obj <- qread('./11_hdWGCNA/11_gr_hdWGCNA_object.qs')
sce.obj@meta.data <- sce.obj@meta.data %>% dplyr::select(-c(57:ncol(.)))

countdata <- GetAssayData(sce.obj, slot = 'counts')

# PCOS cells only; restrict to the hdWGCNA variable genes
meta.df <- sce.obj@meta.data %>% filter(group == 'PCOS')
variable_genes <- sce.obj@misc$GC_wgcna$datExpr %>% colnames(.)
dat <- countdata[variable_genes, rownames(meta.df)] %>% as.data.frame()

result <- scTenifoldKnk(dat,
                        qc = FALSE,
                        gKO = 'Rora',
                        nc_nNet = 10,
                        # ,nc_nCells = 500
                        nc_nComp = 3,
                        nCores = 12)

qsave(result, file = './13_virtual_knockout/13_scTenifoldKnk_result.qs')

# ---- Significantly perturbed genes after Rora KO ------------------------------

result <- qread('./13_virtual_knockout/13_scTenifoldKnk_result.qs')

perturbed_genes <- result$diffRegulation[order(-result$diffRegulation$FC), ] %>%
  filter(p.adj < 0.05)

fwrite(perturbed_genes, './13_virtual_knockout/13_scTenifoldKnk_result.csv', quote = FALSE,
       row.names = TRUE)

# ---- Top 20 perturbed genes ---------------------------------------------------

top_genes <- head(result$diffRegulation[order(-result$diffRegulation$FC), ], 20)

ggplot(top_genes, aes(x = reorder(gene, Z), y = Z)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  coord_flip() +
  labs(title = 'Top 20 Differentitally Regulated Genes\n by Rora KO',
       x = 'Gene', y = 'Z statistic') +
  theme_bw() +
  theme(plot.title = element_text(hjust = .5),
        panel.grid = element_blank())

p_top20 <- ggplot(top_genes, aes(x = Z, y = reorder(gene, Z))) +
  geom_segment(aes(x = 0, xend = Z, y = reorder(gene, Z), yend = reorder(gene, Z)),
               color = 'black', linetype = 2) +
  geom_point(size = 4, pch = 21, fill = '#facd6d', color = 'black') +
  labs(title = 'Top 20 Differentitally Regulated Genes\n by Rora KO',
       y = 'Gene', x = 'Z statistic') +
  scale_x_continuous(expand = c(0, 0), limits = c(0, max(top_genes[['Z']]) + 0.2)) +
  theme_classic() +
  theme(axis.title = element_text(size = 18),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.grid = element_blank(),
        legend.position = c(0.8, 0.9),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.5, 'cm'),
        legend.title = element_text(size = 14),
        plot.title = element_text(hjust = .5, size = 16))

ggsave('./figures/13_scTenifoldKnk_Rora_ko_top20_genes.pdf', p_top20,
       width = 6, height = 4.6)

# KO network plot helpers
source('./11_scTenifoldKnk.R')

pdf('./13_virtual_knockout/13_scTenifoldKnk_ko_Rora.pdf', width = 15, height = 10)
plotKO(result, 'Rora')
dev.off()

plotDR(result)

# ---- CellOracle: perturbation score (PS) on UMAP ------------------------------

library(ggpubr)
library(tidyverse)

# CellOracle grid results
coords_celloracle <- data.table::fread('./13_virtual_knockout/13_celloracle_gridpoints_coordinates.csv', data.table = FALSE)
inner_score_celloracle <- data.table::fread('./13_virtual_knockout/13_celloracle_inner_score.csv', data.table = FALSE)

# Keep grid points passing the mass filter and attach PS
coords_celloracle_filter <- coords_celloracle %>% filter(mass_filter != 'TRUE') %>%
  cbind(., inner_score_celloracle)

# PS map with shift arrows (Rora KO shift direction)
scale <- 2
p_flow <- ggplot(coords_celloracle_filter, aes(x = Dim1, y = Dim2)) +
  geom_point(aes(color = score), size = 2, alpha = 1) +
  geom_segment(aes(x = Dim1,
                   y = Dim2,
                   xend = Dim1 + flow_x * scale,
                   yend = Dim2 + flow_y * scale),
               arrow = arrow(length = unit(0.06, "cm"), type = "open"),
               linewidth = 0.25) +
  scale_color_gradient2(high = '#4d8eb8', low = '#df4e57', name = 'PS', mid = '#f7f7f6') +
  theme_void() +
  theme(legend.position = "right",
        plot.margin = margin(5, 5, 5, 5),
        legend.key.size = unit(.4, 'cm')) +
  coord_fixed(ratio = 1)

ggview:::ggview(p_flow, width = 5, height = 5, bg = 'white')

ggsave('./figures/13_celloracle_PS_PS_UMAP.pdf', p_flow, width = 5, height = 5, bg = 'white')

# Reference (control) flow vectors
ref_scale <- 0.3
p_flow_ref <- ggplot(coords_celloracle_filter, aes(x = Dim1, y = Dim2)) +
  geom_segment(aes(x = Dim1,
                   y = Dim2,
                   xend = Dim1 + ref_flow_x * ref_scale,
                   yend = Dim2 + ref_flow_y * ref_scale),
               arrow = arrow(length = unit(0.08, "cm"), type = "open"),
               linewidth = 0.6,
               alpha = 0.8) +
  scale_color_viridis_c(option = "B", name = "Flow Speed") +
  theme_void() +
  theme(legend.position = "right",
        plot.margin = margin(5, 5, 5, 5),
        legend.key.size = unit(.4, 'cm')) +
  coord_fixed(ratio = 1)

# ---- PS vs pseudotime ---------------------------------------------------------

colorbar_use <- colorRampPalette(colors = c('#df4e57', 'white', '#4d8eb8'))(100)

p_pseudotime_psscore <- ggplot(data = coords_celloracle_filter, aes(x = pseudotime, y = score, fill = score)) +
  geom_point(shape = 21, color = 'black') +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, color = 'red') +
  scale_x_continuous(breaks = seq(0.2, 0.8, 0.2)) +
  scale_fill_gradientn(colors = colorbar_use, name = 'PS') +
  labs(x = 'pseudotime', y = 'Inner product score') +
  theme_bw(base_rect_size = 1.2) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.1, 0.9),
        plot.margin = margin(5, 5, 5, 5),
        legend.key.size = unit(.3, 'cm'),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggview:::ggview(width = 5, height = 5)

ggsave('./figures/13_celloracle_pseudotime_PS_scatterplot.pdf', p_pseudotime_psscore,
       width = 5, height = 5)

# PS distribution across digitized pseudotime bins
coords_celloracle_filter$pseudotime_id <- factor(coords_celloracle_filter$pseudotime_id)
p_boxplot <- ggplot(data = coords_celloracle_filter, aes(x = pseudotime_id, y = score)) +
  geom_jitter(aes(color = score)) +
  geom_boxplot(fill = 'transparent', outliers = FALSE) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.6, color = 'red') +
  scale_color_gradientn(colors = colorbar_use, name = 'PS') +
  labs(x = 'Digitized_pseudotime', y = 'Inner product score') +
  theme_bw(base_rect_size = 1.2) +
  theme(panel.grid = element_blank(),
        panel.background = element_blank(),
        legend.position = c(0.1, 0.9),
        plot.margin = margin(5, 5, 5, 5),
        legend.key.size = unit(.3, 'cm'),
        legend.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16))

ggview:::ggview(width = 6, height = 5)

ggsave('./figures/13_celloracle_pseudotime_PS_boxplot.pdf', p_boxplot, width = 6, height = 5)

# Digitized pseudotime on UMAP
scale <- 2
p_flow <- ggplot(coords_celloracle_filter, aes(x = Dim1, y = Dim2, color = pseudotime_id)) +
  geom_point(size = 2, alpha = 1) +
  ggsci::scale_color_d3(palette = 'category20') +
  labs(color = 'Digitized_pseudotime') +
  theme_void() +
  theme(legend.position = "right",
        plot.margin = margin(5, 5, 5, 5),
        legend.key.size = unit(.4, 'cm')) +
  coord_fixed(ratio = 1)

p_flow

ggview:::ggview(p_flow, width = 6, height = 6, bg = 'white')

ggsave('./figures/13_celloracle_digitized_pseudotime.pdf', p_flow, width = 6, height = 6)
