# Venn diagram: intersection of M8 hub genes, GC_3 markers, and GC_3-specific
# TFs (supporting Rora as the key regulator)

library(qs)
library(tidyverse)
library(data.table)
library(ggvenn)

dir_figures <- './figures'

# hdWGCNA module M8 top-200 hub genes
hub_df <- qread('./11_hdWGCNA/11_hdWGCNA_hub_df_M8.qs')
M8_genes <- hub_df$gene_name

# GC subtype marker genes (PCOS vs Control)
pct_group <- 0.25
log_fc <- 1
adj_p_cutoff <- 0.05
diff_pct_cutoff <- 0.25
GC_dea_ls <- map(1:6, function(a) {
  dea_GC <- readxl::read_excel('./08_GC_subclustering/after_hb_Gr_subannotation_dea_results.xlsx',
                               sheet = a) %>%
    as.data.frame() %>%
    filter(pct_nz_group >= pct_group, logfoldchanges > log_fc, pvals_adj < adj_p_cutoff) %>%
    mutate(diff_pct = pct_nz_group - pct_nz_reference) %>%
    filter(diff_pct < diff_pct_cutoff) %>%
    dplyr::rename(avg_log2FC = logfoldchanges)
  return(dea_GC)
})
names(GC_dea_ls) <- paste0('GC_', 0:5)

GC3_de_genes <- GC_dea_ls[['GC_3']]$names

# Check Rora presence in the M8 / GC_3 overlap
intersect(GC3_de_genes, M8_genes) %>% grep('Rora', ., value = TRUE)

# GC_3-specific regulon TFs from SCENIC
dea_regulon_celltype <- qread('./12_SCENIC/12_SCENIC_celltype_DEA.qs')
GC3_specfic_tf <- dea_regulon_celltype[[4]][[1]] %>% gsub('\\(\\+\\)', '', .)

# ---- Venn diagram -------------------------------------------------------------

venn_data <- list(
  M8_hub_genes = M8_genes,
  GC_3_markers = GC3_de_genes,
  GC3_specfic_TF = GC3_specfic_tf)

p7 <- ggvenn(venn_data,
             names(venn_data), text_size = 6, auto_scale = FALSE,
             set_name_size = 5, show_elements = FALSE,
             show_percentage = TRUE, stroke_size = .4, stroke_linetype = 1,
             stroke_color = "black",
             fill_color = c("#1b9e77", "#9eb9f3", '#f1ce63'),
             set_name_color = c("#1b9e77", "#9eb9f3", '#f1ce63'))
p7

ggview:::ggview(p7, width = 7, height = 7)

ggsave(file.path(dir_figures, 'supplementary_common_genes_vennplot.pdf'),
       p7, width = 7, height = 7)

# Genes shared by all three sets
genes_common <- intersect(GC3_de_genes, M8_genes) %>%
  intersect(., GC3_specfic_tf)

cat(paste0(genes_common, '\n'))
