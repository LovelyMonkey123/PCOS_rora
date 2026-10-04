# ssGSEA validation of the RORA(+) cell signature in microarray cohort GSE98595

library(tidyverse)
library(ggplot2)
library(qs)
library(GSVA)
library(ggpubr)

dir_val <- './19_validation'
dir_figures <- './figures'

# Bulk microarray cohort (PCOS vs control)
cohort_ls <- read_rds(file.path(dir_val, 'PCOS_cohort_overall.RDS'))

group_dat <- cohort_ls[["array_ls"]][["GSE98595"]][["clin_dat"]][["diagnosis"]]
expr_dat <- cohort_ls[["array_ls"]][["GSE98595"]][["expr_dat"]]

# Human-converted RORA(+) cell signature (from 19_external_validation)
genesets <- qread(file.path(dir_val, 'RORA_pos_cells_signature_human.qs'))

gsva_exp <- expr_dat %>% as.matrix()

# GSVA parameter object
ssgseapar <- ssgseaParam(exprData = gsva_exp,
                         geneSets = genesets,
                         minSize = 3,
                         maxSize = 500,
                         normalize = TRUE,
                         use = 'everything')

result <- gsva(ssgseapar)

result <- result %>% as.data.frame() %>% t() %>% as.data.frame()

# ---- Group comparison (control samples only) -----------------------------------

# Sample indices of the untreated group
untreated_indexs <- c(2, 4, 6, 8, 11, 12, 15, 16)

plot_df <- expr_dat %>% t() %>%
  as.data.frame() %>%
  cbind(., group = group_dat, result) %>%
  .[untreated_indexs, ]

plot_df$group <- factor(plot_df$group,
                        levels = names(plot_df$group %>% table()),
                        labels = c('Normal', 'PCOS'))

cols <- c('#1f77b4', '#ff7f0e')
p.GSE98595 <- ggplot(data = plot_df, aes_string(x = 'group', y = 'RORA_pos_cells')) +
  geom_boxplot(outlier.colour = NA, aes_string(fill = "group"), color = 'black',
               size = 0.6, alpha = 0.65) +
  geom_violin(alpha = 0.5, aes_string(fill = "group"), color = NA, trim = TRUE) +
  geom_jitter(shape = 21, size = 1.4, width = 0.3,
              aes_string(fill = "group", color = "group")) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  labs(x = 'GSE98595 \n (lutein GCs from human)', y = NULL,
       title = 'RORA(+)_pos_cells_ssGSEA') +
  stat_compare_means(size = 4) +
  theme_classic(base_rect_size = 1.5) +
  theme(plot.title = element_text(face = 'bold', size = 14, hjust = 0.5),
        legend.text = element_text(size = 10, colour = 'black'),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14, colour = 'black'),
        axis.title.x = element_text(size = 14, colour = 'black'),
        panel.grid = element_blank(),
        legend.position = 'none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.key.size = unit(.3, 'cm'),
        panel.background = element_blank(),
        plot.margin = margin(t = 0.2, r = 0.2, b = 0.1, l = 0.2, unit = 'cm'))

p.GSE98595

ggsave(file.path(dir_figures, '19_GSE98595_RORA_cells_group.pdf'),
       p.GSE98595, width = 3.2, height = 4)
