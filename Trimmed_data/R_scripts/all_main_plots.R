#!/usr/bin/Rscript

# ---------------------------------------------------------------------------- #
# Main plots used in the study
# ---------------------------------------------------------------------------- #

## Libraries
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(patchwork)
library(UpSetR)
# ---------------------------------------------------------------------------- #
# MA plots
# ---------------------------------------------------------------------------- #

par(mfrow = c(2,2))

plotMA(res_WT, main = "WT Case vs Control", colSig = "red", colNonSig = "grey")
plotMA(res_DKO, main = "DKO Case vs Control", colSig = "red", colNonSig = "grey")
plotMA(res_WTvsDKO_control, main = "WT vs DKO Control", colSig = "red", colNonSig = "grey")
plotMA(res_WTvsDKO_infected, main = "WT vs DKO Case", colSig = "red", colNonSig = "grey")

par(mfrow = c(1,1))

# ---------------------------------------------------------------------------- #
#  Summary of DE genes
# ---------------------------------------------------------------------------- #

sig_WT  <- res_WT  %>% as.data.frame() %>% filter(padj < 0.05)
sig_DKO <- res_DKO %>% as.data.frame() %>% filter(padj < 0.05)


data.frame(Comparison   = c("WT Infected vs Control", "DKO Infected vs Control"),
           DE_genes     = c(nrow(sig_WT), nrow(sig_DKO)),
           Upregulated  = c(sum(sig_WT$log2FoldChange >  0),sum(sig_DKO$log2FoldChange > 0)),
           Downregulated = c(sum(sig_WT$log2FoldChange <  0),sum(sig_DKO$log2FoldChange < 0)))

#                Comparison DE_genes Upregulated Downregulated
# 1  WT Infected vs Control     9898        4186          5712
# 2 DKO Infected vs Control     8072        3800          4272
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #

# WT infection response
WT_up <- rownames(res_WT)[
  res_WT$padj < padj_cutoff &
    res_WT$log2FoldChange > lfc_cutoff]

WT_down <- rownames(res_WT)[
  res_WT$padj < padj_cutoff &
    res_WT$log2FoldChange < -lfc_cutoff]

# DKO infection response
DKO_up <- rownames(res_DKO)[
  res_DKO$padj < padj_cutoff &
    res_DKO$log2FoldChange > lfc_cutoff]

DKO_down <- rownames(res_DKO)[
  res_DKO$padj < padj_cutoff &
    res_DKO$log2FoldChange < -lfc_cutoff]

# Interaction-significant genes
INT_sig <- rownames(res_interaction)[
  res_interaction$padj < padj_cutoff]

# cleaning all the sets we obtained up
WT_up    <- unique(na.omit(WT_up))
WT_down  <- unique(na.omit(WT_down))
DKO_up   <- unique(na.omit(DKO_up))
DKO_down <- unique(na.omit(DKO_down))
INT_sig  <- unique(na.omit(INT_sig))

# WT-specific infection response
WT_specific_up   <- setdiff(WT_up, DKO_up)
WT_specific_down <- setdiff(WT_down, DKO_down)

# DKO-specific response
DKO_specific_up   <- setdiff(DKO_up, WT_up)
DKO_specific_down <- setdiff(DKO_down, WT_down)

# Common responses
Common_up   <- intersect(WT_up, DKO_up)
Common_down <- intersect(WT_down, DKO_down)

# Genes whose induction is lost in DKO
Lost_in_DKO <- intersect(
  WT_up,
  rownames(res_interaction)[
    res_interaction$padj < padj_cutoff &
      res_interaction$log2FoldChange < 0
  ]
)

# Genes gained or exaggerated in DKO
Gained_in_DKO <- intersect(DKO_up,rownames(res_interaction)[
  res_interaction$padj < padj_cutoff &
    res_interaction$log2FoldChange > 0])

gene_sets <- list(
  WT_Up       = WT_up,
  WT_Down     = WT_down,
  DKO_Up      = DKO_up,
  DKO_Down    = DKO_down,
  Interaction = INT_sig)

upset(fromList(gene_sets),
      order.by = "freq",
      mainbar.y.label = "Gene intersections",
      sets.x.label    = "Genes per condition")


# ---------------------------------------------------------------------------- #
# Volcano plot
# ---------------------------------------------------------------------------- #

plot_volcano <- function(
    res,
    title,
    pCutoff = 0.05,
    FCcutoff = 1,
    pointSize = 2.5,
    labSize = 3
) {
  
  res <- as.data.frame(res)
  
  EnhancedVolcano(
    res,
    lab = rownames(res),
    x = "log2FoldChange",
    y = "padj",
    title = title,
    subtitle = NULL,
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    pointSize = pointSize,
    labSize = labSize
  ) +
    theme(
      plot.title   = element_text(hjust = 0.5),
      legend.text  = element_text(size = 9),
      legend.title = element_text(size = 10)
    )
}

# ---------------------------------------------------------------------------- #
volcano_WT <- plot_volcano(
  res_WT,
  title = "WT Infected vs Control"
)

volcano_DKO <- plot_volcano(
  res_DKO,
  title = "DKO Infected vs Control"
)

volcano_WTvsDKO_ctrl <- plot_volcano(
  res_WTvsDKO_control,
  title = "WT vs DKO (Control)"
)

volcano_WTvsDKO_inf <- plot_volcano(
  res_WTvsDKO_infected,
  title = "WT vs DKO (Infected)"
)

# ---------------------------------------------------------------------------- #
# Combining plots 
# ---------------------------------------------------------------------------- #

volcano_panel <-
  (volcano_WT | volcano_DKO) /
  (volcano_WTvsDKO_ctrl | volcano_WTvsDKO_inf) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Display plot
volcano_panel
