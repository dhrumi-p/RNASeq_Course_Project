#!/usr/bin/Rscript

## Libraries
library(dplyr)
library(ggplot2)
library(tidyverse)

## ---------------------------------------------------------------------- ##
genes_of_interest <- c(
  "ENSMUSG00000055170", # Ifng
  "ENSMUSG00000000386", # Mx1
  "ENSMUSG00000029417", # Cxcl9
  "ENSMUSG00000037321"  # Tap1
)

gene_symbols <- c("Ifng", "Mx1", "Cxcl9", "Tap1")

## ---------------------------------------------------------------------- ##

extract_lfc <- function(res, contrast_name) {
  res %>%
    as.data.frame() %>%
    rownames_to_column("ENSEMBL") %>%
    filter(ENSEMBL %in% genes_of_interest) %>%
    mutate(
      Gene = gene_symbols[match(ENSEMBL, genes_of_interest)],
      Contrast = contrast_name,
      negLog10Padj = -log10(padj)
    )
}


lfc_all <- bind_rows(
  extract_lfc(res_WT, "WT: Case vs Control"),
  extract_lfc(res_DKO, "DKO: Case vs Control"),
  extract_lfc(res_WTvsDKO_infected, "WT vs DKO (Case)")
)

lfc_all

## ---------------------------------------------------------------------- ##
ggplot(lfc_all,
       aes(x = Gene,
           y = log2FoldChange,
           fill = log2FoldChange > 0)) +
  geom_col(
    width = 0.6,
    color = "black",
    linewidth = 0.3
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_text(
    aes(label = paste0("adj p = ", formatC(padj, format = "e", digits = 1))),
    vjust = ifelse(lfc_all$log2FoldChange > 0, -0.8, 1.6),
    size = 3,
    color = "black"
  ) +
  facet_wrap(~ Contrast, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = c("TRUE" = "tomato", "FALSE" = "cyan2"),
    guide = "none"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.15, 0.25))
  ) +
  labs(
    title = "Gene expression across all conditions",
    subtitle = "Log2 fold change with adjusted p-values",
    y = "Log2 fold change",
    x = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "italic", size = 11),
    axis.title.y = element_text(size = 13)
  )
## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##
norm_counts <- counts(dds, normalized = TRUE)

norm_df <- norm_counts[genes_of_interest, ] %>%
  as.data.frame() %>%
  rownames_to_column("ENSEMBL") %>%
  pivot_longer(
    cols = -ENSEMBL,
    names_to = "Sample",
    values_to = "NormalisedCounts"
  ) %>%
  mutate(
    Gene = gene_symbols[match(ENSEMBL, genes_of_interest)]
  ) %>%
  left_join(
    meta_data %>%
      rownames_to_column("Sample"),
    by = "Sample"
  )

norm_df <- norm_df %>%
  mutate(
    Condition = factor(
      paste(Genotype, Infection),
      levels = c(
        "WT Control",
        "WT Infected",
        "DKO Control",
        "DKO Infected"
      )
    )
  )

## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##
ggplot(norm_df,
       aes(x = Condition,
           y = NormalizedCounts,
           fill = Infection)) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.75,
    width = 0.65
  ) +
  geom_jitter(
    width = 0.15,
    size = 2,
    alpha = 0.8
  ) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(
    values = c("Control" = "grey",
               "Infected" = "red")
  ) +
  labs(
    title = "Normalised expression of genes",
    subtitle = "DESeq2 normalised counts",
    y = "Normalised counts",
    x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.title = element_text(face = "bold")
  )+ scale_y_log10()


## ---------------------------------------------------------------------- ##
## ---------------------------------------------------------------------- ##

