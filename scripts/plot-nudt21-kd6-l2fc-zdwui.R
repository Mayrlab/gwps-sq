library(tidyverse)
library(magrittr)
library(ggbeeswarm)


df_wui <- readRDS("tbl/df_wui_kd6_essential_ui10.Rds")

df_nudt21 <- filter(df_wui, target_gene == "NUDT21") %>%
    filter(!is.na(wui))


df_ntp <- filter(df_wui, target_gene == "non-targeting")

df_means <- df_ntp %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(is.na(wui)) < 0.1) %>%
    summarize(tpm_mean=Hmisc::wtd.mean(tpm_gene, n_cells),
              tpm_sd=sqrt(Hmisc::wtd.var(tpm_gene, n_cells)),
              wui_mean=Hmisc::wtd.mean(wui, n_cells),
              wui_sd=sqrt(Hmisc::wtd.var(wui, n_cells)),
              .groups='drop') %>%
    filter(!is.na(wui_mean))

df_means %>%
    ggplot(aes(x=tpm_mean, y=wui_mean)) +
    geom_point() +
    scale_x_log10()

df_dwui <- df_nudt21 %>%
    inner_join(df_means, by=c("gene_id", "gene_name")) %>%
    mutate(z_dtpm=(tpm_gene-tpm_mean)/tpm_sd,
           l2fc=log2(tpm_gene/tpm_mean),
           dwui=wui-wui_mean,
           z_dwui=dwui/wui_sd,
           z_dwui_adjust=ifelse(abs(z_dwui) > 10, sign(z_dwui)*10, z_dwui),
           l2fc_adjust=ifelse(abs(l2fc) > 3, sign(l2fc)*3, l2fc))

df_dwui %>%
    ggplot(aes(x="NUDT21", y=l2fc)) +
    geom_boxplot()

df_dwui %>%
    mutate(outlier=abs(l2fc) >= 1) %>%
    ggplot(aes(x="NUDT21", y=l2fc_adjust)) +
    geom_hline(yintercept=c(-1,1), linetype='dashed') +
    geom_violin(draw_quantiles=c(0.25,0.5,0.75), color='grey20', fill='grey90', alpha=0.8) +
    geom_quasirandom(aes(color=outlier, size=outlier)) +
    scale_color_manual(values=c('black', 'plum1'), guide=NULL) +
    scale_size_manual(values=c(0.05, 0.3), guide=NULL) +
    labs(x=NULL, y="log2 fold-change") +
    cowplot::theme_minimal_hgrid()

ggsave("img/target-readouts/nudt21-l2fc-kd6.pdf", width=3, height=4)

df_dwui %>%
    mutate(outlier=abs(l2fc) >= 1 | abs(z_dwui) >= 3) %>%
    ggplot(aes(x=z_dwui_adjust, y=l2fc_adjust)) +
    geom_vline(xintercept=0, color='grey80') +
    geom_hline(yintercept=c(-1,1), linetype='dashed') +
    geom_vline(xintercept=c(-3,3), linetype='dashed') +
    geom_point(aes(color=outlier, size=outlier)) +
    scale_color_manual(values=c('black', 'plum3'), guide=NULL) +
    scale_size_manual(values=c(0.03, 0.35), guide=NULL) +
    labs(x="Z-score dWUI", y="log2 fold-change") +
    cowplot::theme_minimal_hgrid()

ggsave("img/target-readouts/nudt21-zdwui-l2fc-kd6.pdf", width=4, height=4)

