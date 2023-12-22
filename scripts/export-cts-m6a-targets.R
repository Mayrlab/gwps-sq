library(tidyverse)
library(magrittr)
library(ggbeeswarm)

## genes of interest
GOIS=c("DCP1A", "TGFBR1", "PURA", "BAG3")

## SummarizedExperiment object
se_kd6 <- readRDS("data/se/kd6_essential_bulk_expressed.Rds")

## subset to txs of interest
idx_gois <- rowData(se_kd6)$gene_name %in% GOIS
se_select <- se_kd6[idx_gois,]

## extract counts and pivot to long form
df_cts <- assay(se_select, "counts") %>%
    as.matrix() %>%
    as_tibble(rownames="transcript_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="umi_count")

## attach annotation inform
df_rowdata <- rowData(se_select) %>% as_tibble %>%
    dplyr::select(transcript_id, utr_name, utr_rank, gene_id, gene_name, is_novel)

df_coldata <- colData(se_select) %>% as_tibble

df_final <- right_join(df_rowdata, df_cts, by="transcript_id") %>%
    left_join(df_coldata, by="sgID_AB")

## export
saveRDS(df_final, "data/counts/df_m6a_targets.kd6.Rds")

## Example Plot
df_final %>%
    filter(gene_name=="BAG3",
           target_gene %in% c("CBLL1", "YTHDC1", "non-targeting")) %>%
    group_by(sgID_AB) %>%
    mutate(ui=umi_count/sum(umi_count)) %>%
    ungroup() %>%
    mutate(is_targeted=target_gene != 'non-targeting') %>%
    ggplot(aes(x=utr_name, y=ui, color=is_targeted)) +
    geom_quasirandom() +
    scale_color_manual(values=c("FALSE"="black", "TRUE"="red")) +
    theme_bw()

