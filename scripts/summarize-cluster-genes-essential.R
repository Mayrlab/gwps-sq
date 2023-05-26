library(tidyverse)
library(writexl)

df_targets_raw <- read_csv("tbl/df_target_raw_nnclusters_kd6_essential_ui10.csv")
df_targets_sig <- read_csv("tbl/df_target_nnclusters_kd6_essential_ui10.csv")

targets_all <- unique(df_targets_raw$target_gene)
targets_sig <- unique(df_targets_sig$target_gene)

df_genes_raw <- read_csv("tbl/df_gene_nnclusters_kd6_essential_ui10.csv")

genes_all <- unique(df_genes_raw$gene_name)
genes_sig <- filter(df_genes_raw, !(cluster_id %in% c(9,12,13))) %$% unique(gene_name)

write_xlsx(list(perturbations_in=data.frame(gene_symbol=targets_all),
                perturbations_out=data.frame(gene_symbol=targets_sig),
                multi_in=data.frame(gene_symbol=genes_all),
                multi_out=data.frame(gene_symbol=genes_sig)),
           path="tbl/genes_kd6_essential.xlsx")
