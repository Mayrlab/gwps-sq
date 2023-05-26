library(tidyverse)
library(magrittr)

MAX_NTP_NA=0.20
## min avg TPM in NTPs
MIN_MEAN_TPM=20

df_wui <- readRDS("tbl/df_wui_kd8_gwps_ui10.Rds")

df_ntp <- filter(df_wui, target_gene == 'non-targeting') %>%
    group_by(gene_id, gene_name) %>%
    filter(mean(is.na(wui)) <= MAX_NTP_NA) %>%
    filter(mean(tpm) >= MIN_MEAN_TPM) %>%
    filter(!is.na(wui)) %>%
    summarize(mean_wui=weighted.mean(wui, n_cells), .groups='drop',
              sd_wui=sqrt(sum((wui-mean_wui)^2)/n()))

plot_zdwui_ma <- function (df_zdwui) {
    df_zdwui %>%
        ggplot(aes(x=tpm, zdwui)) +
        geom_hline(yintercept=0, linetype="dashed", color="magenta") +
        geom_point() +
        scale_x_log10() +
        scale_y_continuous(limits=c(-4,4)) +
        labs(x="TPM", y="dWUI [Z-scale]") +
        facet_wrap(facets=vars(n_cells)) +
        theme_bw()
}

### PRDM9
df_prdm9 <- df_wui %>%
    filter(target_gene == "PRDM9") %>%
    inner_join(df_ntp, by=c("gene_id", "gene_name")) %>%
    mutate(dwui=wui-mean_wui, zdwui=dwui/sd_wui) %>%
    dplyr::select(gene_id, sgID_AB, n_cells, wui, dwui, zdwui, tpm)

plot_zdwui_ma(df_prdm9)

### RAD52
df_rad52 <- df_wui %>%
    filter(target_gene == "RAD52") %>%
    inner_join(df_ntp, by=c("gene_id", "gene_name")) %>%
    mutate(dwui=wui-mean_wui, zdwui=dwui/sd_wui) %>%
    dplyr::select(gene_id, sgID_AB, n_cells, wui, dwui, zdwui, tpm)

plot_zdwui_ma(df_rad52)

df_rad52 %>%
    dplyr::select(gene_id, n_cells, zdwui) %>%
    pivot_wider(names_from=n_cells, names_prefix="zdwui_", values_from=zdwui) %>%
    ggplot(aes(x=zdwui_242, y=zdwui_190)) +
    geom_point()

### ADAMTS3 (actually would like to see ADAMTS9)
df_adamts3 <- df_wui %>%
    filter(target_gene == "ADAMTSL4") %>%
    inner_join(df_ntp, by=c("gene_id", "gene_name")) %>%
    mutate(dwui=wui-mean_wui, zdwui=dwui/sd_wui) %>%
    dplyr::select(gene_id, sgID_AB, n_cells, wui, dwui, zdwui, tpm)

plot_zdwui_ma(df_adamts3)

