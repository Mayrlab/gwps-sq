library(magrittr)
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)
library(plyranges)

MIN_MEAN_UI=0.10
FILE_SE_TX="data/se/kd6_essential_bulk_expressed.Rds"
FILE_GTF="../hcl/data/gff/utrome.e30.t5.gc39.pas3.f0.9999.w500.gtf.gz"

FILE_WUI_OUT=sprintf("tbl/df_wui_kd6_essential_ui%2.0f.Rds", 100*MIN_MEAN_UI)
FILE_TPM_OUT="tbl/df_tpm_kd6_essential.Rds"
FILE_UTRS_OUT=sprintf("tbl/df_utrs_kd6_essential_ui%2.0f.csv", 100*MIN_MEAN_UI)

## Load data
se_tx_target <- readRDS(FILE_SE_TX)

df_utr3 <- read_gff(FILE_GTF, col_names=c("type", "transcript_id", "gene_id")) %>%
    filter(type == 'transcript') %>%
    select("gene_id", "transcript_id") %>%
    anchor_3p() %>%
    mutate(width=0) %>%
    as_tibble() %>%
    select("gene_id", "transcript_id", "seqnames", "end", "strand")

###########################################
## DETERIME MULTI-UTR GENES
###########################################

idx_ntp <- colnames(se_tx_target) %>% str_detect("non-targeting")
idx_multi_raw <- which(rowData(se_tx_target)$utr_type_no_ipa == 'multi' & !rowData(se_tx_target)$is_ipa)

cts_tx_ntp <- assay(se_tx_target[idx_multi_raw, idx_ntp], 'counts')

M_gene_tx <- rowData(se_tx_target[idx_multi_raw,]) %>%
    as_tibble %>%
    dplyr::select(transcript_id, gene_id) %>%
    deframe %>%
    fac2sparse

cts_gene_ntp <- M_gene_tx %*% cts_tx_ntp

ui_tx_ntp <- as.matrix(cts_tx_ntp / (t(M_gene_tx) %*% cts_gene_ntp))

ncells_ntp <- se_tx_target$n_cells[idx_ntp]

mui_tx <- rowWeightedMeans(ui_tx_ntp, ncells_ntp, na.rm=TRUE)
mui_tx[is.na(mui_tx)] <- 0

idx_tx <- M_gene_tx[,mui_tx >= MIN_MEAN_UI] %>% {
    idx_i <- rowSums(.) > 1;
    idx_j <- colSums(.[idx_i,]) > 0;
    .[idx_i, idx_j]
} %>% colnames

df_multi <- rowData(se_tx_target[idx_tx,]) %>% as_tibble %>%
    mutate(ensembl_id=str_remove(gene_id, "\\.[0-9]+$"),
           mui=mui_tx[transcript_id]) %>%
    group_by(ensembl_id) %>%
    mutate(utr_rank_expr=rank(utr_rank),
           utr_wt=(utr_rank_expr - 1)/(max(utr_rank_expr) - 1),
           is_top2=rank(-mui) <= 2,
           is_su=utr_rank == min(utr_rank[is_top2]),
           is_lu=utr_rank == max(utr_rank[is_top2])) %>%
    ungroup() %>%
    dplyr::select(transcript_id, transcript_name, gene_id, gene_name, ensembl_id,
                  utr_rank, utr_rank_expr, utr_wt, mui, is_su, is_lu) %>%
    arrange(gene_id, utr_rank)

## use ordered version
idx_tx <- df_multi$transcript_id

if (interactive()) {
    df_multi %>% distinct(gene_id) %>% nrow %>%
        sprintf(fmt="Found %d multi-UTR genes.") %>%
        cat()

    df_multi %>% nrow %>%
        sprintf(fmt=" These consist of %d 3' UTR isoforms.") %>%
        cat(fill=TRUE)
}

## export transcripts
left_join(df_multi, df_utr3, by=c("transcript_id", "gene_id")) %>%
    write_csv(FILE_UTRS_OUT)

###########################################
## COMPUTE TPMs
###########################################

cts_tx_target <- assay(se_tx_target, "counts")

tpm_tx_target <- cts_tx_target %>%
    { . %*% Diagonal(ncol(.), 1e6/colSums(.)) }

M_gene_tx <- rowData(se_tx_target) %>% as_tibble %>%
    dplyr::select(transcript_id, gene_id) %>%
    deframe %>%
    fac2sparse

cts_gene_target <- M_gene_tx %*% cts_tx_target

tpm_gene_target <- M_gene_tx %*% tpm_tx_target

###########################################
## COMPUTE SU & LU TPMs
###########################################

idx_su <- filter(df_multi, is_su) %>%
    dplyr::select(gene_id, transcript_id) %>%
    deframe()

idx_lu <- filter(df_multi, is_lu) %>%
    dplyr::select(gene_id, transcript_id) %>%
    deframe()

tpm_su_gene_target <- tpm_tx_target[idx_su,] %>%
    `rownames<-`(names(idx_su))

tpm_lu_gene_target <- tpm_tx_target[idx_lu,] %>%
    `rownames<-`(names(idx_lu))

###########################################
## COMPUTE IPA
###########################################

idx_ipa_txs <- rowData(se_tx_target)$is_ipa
idx_ipa_genes <- rowSums(M_gene_tx[, idx_ipa_txs]) > 0

ipa_gene_target <- (M_gene_tx[idx_ipa_genes, idx_ipa_txs] %*% cts_tx_target[idx_ipa_txs, ]) / cts_gene_target[idx_ipa_genes,]

################
## COMPUTE WUI
################
wt_tx <- df_multi %>% dplyr::select(transcript_id, utr_wt) %>% deframe

M_gene_tx <- df_multi %>%
    dplyr::select(transcript_id, gene_id) %>%
    deframe %>%
    fac2sparse

cts_tx_target <- assay(se_tx_target[idx_tx,], "counts")

cts_gene_target <- M_gene_tx %*% cts_tx_target

ui_tx_target <- cts_tx_target / (t(M_gene_tx) %*% cts_gene_target)
#ui_tx_target[is.na(ui_tx_target)] <- 0

wui_gene_target <- M_gene_tx %*% Diagonal(length(wt_tx), wt_tx) %*% ui_tx_target

########################
## CONVERT TO TIBBLES ##
########################

df_ipa <- ipa_gene_target %>% as.matrix %>%
    as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="pct_ipa", values_drop_na=TRUE)

df_wui <- wui_gene_target %>% as.matrix %>%
    as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="wui", values_drop_na=TRUE)

df_tpm_gene <- tpm_gene_target[union(rownames(wui_gene_target), rownames(ipa_gene_target)),] %>%
    as.matrix %>% as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="tpm_gene")

df_tpm_su <- tpm_su_gene_target %>%
    as.matrix %>% as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="tpm_su")

df_tpm_lu <- tpm_lu_gene_target %>%
    as.matrix %>% as_tibble(rownames="gene_id") %>%
    pivot_longer(cols=-1, names_to="sgID_AB", values_to="tpm_lu")

df_gene <- rowData(se_tx_target) %>%
    as_tibble %>%
    distinct(gene_id, gene_name)

df_tbl <- left_join(df_tpm_gene, df_wui, by=c("gene_id", "sgID_AB")) %>%
    left_join(df_tpm_su, by=c("gene_id", "sgID_AB")) %>%
    left_join(df_tpm_lu, by=c("gene_id", "sgID_AB")) %>%
    left_join(df_ipa, by=c("gene_id", "sgID_AB")) %>%
    left_join(df_gene, by="gene_id") %>%
    left_join(as_tibble(colData(se_tx_target)), by="sgID_AB") %>%
    dplyr::select(target_gene, target_gene_id, sgID_AB, n_cells,
           gene_id, gene_name, tpm_gene, tpm_su, tpm_lu, wui, pct_ipa) %>%
    mutate(across(where(is.character), factor))

df_tpm_gene_all <- tpm_gene_target %>%
    as.matrix %>% as_tibble(rownames="gene_id") %>%
    left_join(df_gene, by="gene_id") %>%
    dplyr::select(gene_id, gene_name, everything()) %>%
    pivot_longer(cols=-c(1,2), names_to="sgID_AB", values_to="tpm_gene") %>%
    dplyr::select(gene_id, gene_name, sgID_AB, tpm_gene) %>%
    mutate(across(where(is.character), factor))

####################
## EXPORT
####################

saveRDS(df_tbl, FILE_WUI_OUT)
saveRDS(df_tpm_gene_all, FILE_TPM_OUT)
