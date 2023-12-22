library(magrittr)
library(tidyverse)
library(writexl)
library(SummarizedExperiment)
library(Matrix)

FILE_CLUSTERS="tbl/df_target_nnclusters_kd6_essential_ui10.csv"
FILE_CLUSTERS_NEW="metadata/kd6-cluster-info-sm.csv"
FILE_SE="data/se/kd6_essential_bulk_expressed.Rds"

FILE_OUT_XLSX="tbl/counts_tx_cluster_kd6.xlsx"
FILE_OUT_RDS="tbl/df_counts_tx_cluster_kd6.Rds"

## ====================
## LOAD DATA
## ====================

df_clusters <- read_csv(FILE_CLUSTERS) %>%
    mutate(cluster_id=as.character(cluster_id))

old2new <- read_csv(FILE_CLUSTERS_NEW, col_types='c___c__') %>%
    deframe %>% `[<-`("NTP", "NTP")

se <- readRDS(FILE_SE)

idx_ntp <- se$target_gene == 'non-targeting'
idx_clusters <- se$sgID_AB %in% df_clusters$sgID_AB

## filter targets of interest
se %<>% `[`(,idx_clusters | idx_ntp)
## filter expressed transcripts
se %<>% `[`(rowSums(assay(., 'counts')) > 0,)

## attach cluster IDs
colData(se) %<>% as_tibble %>%
    left_join(df_clusters, by=c("target_gene", "target_gene_id", "sgID_AB")) %>%
    mutate(cluster_id=ifelse(is.na(cluster_id), "NTP", cluster_id)) %>%
    DataFrame(row.names=.$sgID_AB)

## ====================
## SUMMARIZE CLUSTERS
## ====================

idx_clusters <- c("NTP", 1:18)
M_target_cluster <- fac2sparse(se$cluster_id) %>% t %>%
    `colnames<-`(old2new[colnames(.)]) %>%
    `[`(,idx_clusters)

cts_tx_cluster <- assay(se, 'counts') %*% M_target_cluster

## ====================
## CONVERT DATAFRAME
## ====================

df_rd <- rowData(se) %>% as_tibble %>%
    dplyr::select(any_of(c('transcript_id', 'gene_name', 'gene_id', 'utr_rank')))

df_txs <- as.matrix(cts_tx_cluster) %>%
    as_tibble(rownames='transcript_id') %>%
    inner_join(x=df_rd, by='transcript_id') %>%
    arrange(gene_name, utr_rank)

df_txs_long <- df_txs %>%
    pivot_longer(cols=all_of(idx_clusters), names_to="cluster", values_to="counts")

## ====================
## EXPORT
## ====================

df_txs %>% write_xlsx(path=FILE_OUT_XLSX)
saveRDS(df_txs_long, file=FILE_OUT_RDS)


