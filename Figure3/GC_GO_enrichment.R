# GO enrichment of marker genes for each granulosa cell subtype

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(cowplot)
library(qs)
library(yulab.utils)

dir_gc <- './08_GC_subclustering'
dir_figures <- './figures'

# ---- Marker genes per GC subtype ----------------------------------------------
pct_group <- 0.25
log_fc <- 0.585
adj_p_cutoff <- 0.05
diff_pct_cutoff <- 0.25
GC_dea_ls <- map(1:6, function(a) {
  dea_GC <- readxl::read_excel(file.path(dir_gc, 'after_hb_Gr_subannotation_dea_results.xlsx'),
                               sheet = a) %>%
    as.data.frame() %>%
    filter(pct_nz_group >= pct_group, logfoldchanges > log_fc, pvals_adj < adj_p_cutoff) %>%
    mutate(diff_pct = pct_nz_group - pct_nz_reference) %>%
    filter(diff_pct < diff_pct_cutoff) %>%
    dplyr::rename(avg_log2FC = logfoldchanges)
  return(dea_GC)
})
names(GC_dea_ls) <- paste0('GC_', 0:5)

# ---- GO-BP enrichment per subtype ---------------------------------------------

eres.list <- lapply(seq_along(GC_dea_ls), function(clu) {
  genes <- GC_dea_ls[[clu]][['names']]
  enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = 'SYMBOL', ont = 'BP',
           pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 1000) %>%
    clusterProfiler::simplify(., cutoff = 0.7, by = "p.adjust", select_fun = min)
})
names(eres.list) <- names(GC_dea_ls)

qsave(eres.list, file = file.path(dir_gc, '08_Gr_subannotation_GO_res.qs'))

eres.list <- qread(file.path(dir_gc, '08_Gr_subannotation_GO_res.qs'))

# Top 20 terms per subtype (ranked by p.adjust)
go.df.ls <- names(eres.list) %>%
  lapply(., function(xx) {
    res <- eres.list[[xx]] %>% as.data.frame() %>%
      filter(p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      slice_head(n = 20)
    Description.levels <- unique(res$Description)
    res$Description <- factor(res$Description, levels = Description.levels %>% rev())
    res
  })
names(go.df.ls) <- names(eres.list)

# Top 5 terms per subtype (ranked by RichFactor) for plotting
go.df.ls1 <- go.df.ls %>%
  lapply(., function(xx) {
    xx <- xx %>%
      .[1:5, ] %>%
      arrange(desc(RichFactor))
    res <- xx
    return(res)
  })

# Subtype colors
rowcolor <- c('#1f77b4', '#ff7f0e', '#2c9e6b', '#d62728', '#aa42fc', '#8c564b')
names(rowcolor) <- names(go.df.ls1)

# GO bar plot per subtype
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

  Description <- go.df$Description %>%
    as.character() %>%
    map_chr(., function(char.use) {
      ifelse(nchar(char.use) > 60, yulab.utils::str_wrap(char.use, 60), char.use)
    })
  go.df$Description <- Description
  Description.levels <- unique(go.df$Description)
  go.df$Description <- factor(go.df$Description, levels = Description.levels %>% rev())

  ggplot() +
    geom_bar(data = go.df,
             aes(x = RichFactor, y = Description),
             color = rowcolor[M.name],
             width = 0.8,
             stat = 'identity',
             fill = 'transparent') +
    scale_x_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.y = element_blank()) +
    geom_text(data = go.df,
              aes(x = 0, y = Description, label = Description),
              size = 5, color = 'black', hjust = -0.01) +
    labs(x = 'GeneRatio', y = NULL, title = M.name) +
    mytheme
})

qsave(enrich.bar.ls, file = file.path(dir_gc, '08_Gr_subannotation_enrich_plot_ls.qs'))

enrich.bar.ls_m <- enrich.bar.ls

# Combine all subtype panels (2 rows x 3 columns)
p.comb <- cowplot::plot_grid(plotlist = enrich.bar.ls_m,
                             ncol = 3, byrow = TRUE, nrow = 2)

ggsave(file.path(dir_figures, '08_Gr_subannotation_GO_function.pdf'), p.comb,
       width = 6 * 3, height = 3 * 2)
