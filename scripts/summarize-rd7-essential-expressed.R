library(magrittr)
library(tidyverse)
library(Matrix)
library(SingleCellExperiment)

FILE_SCE_TX="data/sce/rd7_essential.annot.txs.Rds"

## Load data
sce <- readRDS(FILE_SCE_TX)

M_cell_group <- fac2sparse(sce$sgID_AB) %>% t

cts_tx_group <- counts(sce) %*% M_cell_group

cts_tx_group %<>% `[`(rowSums(.) > 0, )

df_rowdata <- rowData(sce[rownames(cts_tx_group),])

df_coldata <- colData(sce) %>% as_tibble %>%
    group_by(target_gene, target_gene_id, sgID_AB) %>%
    summarize(n_cells=dplyr::n(), .groups='drop') %>%
    DataFrame(row.names=.$sgID_AB) %>%
    `[`(colnames(cts_tx_group),)

se_tx_group <- SummarizedExperiment(assays=list(counts=cts_tx_group),
                                    rowData=df_rowdata,
                                    colData=df_coldata)

## export
saveRDS(se_tx_group, "data/se/rd7_essential_bulk_expressed.Rds")

