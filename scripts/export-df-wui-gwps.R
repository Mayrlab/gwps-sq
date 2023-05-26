library(magrittr)
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)

FILE_SE_TX="data/se/kd8_gwps_bulk_expressed.Rds"
FILE_OUT="tbl/df_wui_kd8_gwps.Rds"
MIN_UMIS=5000

## Load data
se_tx_target <- readRDS(FILE_SE_TX)

cts_tx_target <- assay(se_tx_target, "counts")

M_gene_tx <- rowData(se_tx_target) %>% as_tibble %>%
    dplyr::select(transcript_id, gene_id) %>%
    deframe %>%
    fac2sparse

cts_gene_target <- M_gene_tx %*% cts_tx_target

tpm_gene_target <- cts_gene_target %>%
    { . %*% Diagonal(ncol(.), 1e6/colSums(.)) }

##################
## IDENTIFY MULTI
##################

df_multi <- rowData(se_tx_target) %>% as_tibble %>%
    filter(utr_type_raw == 'multi') %>%
    filter(rowSums(cts_tx_target[transcript_id,]) >= MIN_UMIS) %>%
    group_by(gene_id) %>%
    filter(dplyr::n() > 1) %>%
    mutate(utr_rank_expr=rank(utr_rank),
           utr_wt=(utr_rank_expr - 1)/(max(utr_rank_expr) - 1)) %>%
    ungroup() %>%
    dplyr::select(transcript_id, transcript_name, gene_id, gene_name,
                  utr_rank, utr_rank_expr, utr_wt) %>%
    arrange(gene_id, utr_rank)

idx_tx <- df_multi$transcript_id

if (interactive()) {
    df_multi %>% distinct(gene_id) %>% nrow %>%
        sprintf(fmt="Found %d multi-UTR genes.") %>%
        cat()

    df_multi %>% nrow %>%
        sprintf(fmt=" These consist of %d 3' UTR isoforms.") %>%
        cat(fill=TRUE)
}

################
## COMPUTE WUI
################
wt_tx <- df_multi %>% dplyr::select(transcript_id, utr_wt) %>% deframe

M_gene_tx <- df_multi %>%
    dplyr::select(transcript_id, gene_id) %>%
    deframe %>%
    fac2sparse

cts_gene_target <- M_gene_tx %*% cts_tx_target[idx_tx,]

ui_tx_target <- cts_tx_target[idx_tx,] / (t(M_gene_tx) %*% cts_gene_target)
#ui_tx_target[is.na(ui_tx_target)] <- 0

wui_gene_target <- M_gene_tx %*% Diagonal(length(wt_tx), wt_tx) %*% ui_tx_target


####################
## CONVERT TO TIBBLE
####################

df_wui <- wui_gene_target %>% as.matrix %>%
    as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="wui", values_drop_na=TRUE)

df_tpm <- tpm_gene_target[rownames(wui_gene_target),] %>%
    as.matrix %>% as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="tpm")

df_gene <- df_multi %>%
    group_by(gene_id, gene_name) %>%
    summarize(n_utrs=dplyr::n(), .groups='drop')

df_tbl <- left_join(df_wui, df_tpm, by=c("gene_id", "sgID_AB")) %>%
    left_join(df_gene, by="gene_id") %>%
    left_join(as_tibble(colData(se_tx_target)), by="sgID_AB") %>%
    dplyr::select(target_gene, target_gene_id, sgID_AB, n_cells,
           gene_id, gene_name, n_utrs, wui, tpm) %>%
    mutate(across(where(is.character), factor))

####################
## EXPORT
####################

saveRDS(df_tbl, FILE_OUT)
