library(Matrix)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(tidyverse)
library(magrittr)

cts_tx_group <- read_rds("data/counts/kd8_gwps.cts_tx_group.Rds") %>%
    # only detected transcripts
    `[`(rowSums(.) > 0,)

sce <- loadHDF5SummarizedExperiment("data/sce", "kd8_gwps.annot.txs.")

df_groups <- colData(sce) %>% as_tibble %>%
    group_by(target_gene, target_gene_id, sgID_AB) %>%
    summarize(n_cells=dplyr::n(), .groups='drop') %>%
    DataFrame(row.names=.$sgID_AB) %>%
    `[`(colnames(cts_tx_group),)

df_txs <- rowData(sce) %>% `[`(rownames(cts_tx_group),)

se <- SummarizedExperiment(assays=list(counts=cts_tx_group),
                           rowData=df_txs, colData=df_groups)

saveRDS(se, "data/se/kd8_gwps_bulk_expressed.Rds")
