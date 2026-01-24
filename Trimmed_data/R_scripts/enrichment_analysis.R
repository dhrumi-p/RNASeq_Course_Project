#!/usr/bin/Rscript


# ---------------------------------------------------------------------------- #
# GO enrichment code 
# ---------------------------------------------------------------------------- # 
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
# ---------------------------------------------------------------------------- # 

# WT samples

sig_WT <- as.data.frame(res_WT) %>%
  filter(!is.na(padj), padj < 0.05)

WT_up <- rownames(sig_WT[sig_WT$log2FoldChange > 1, ])
WT_down <- rownames(sig_WT[sig_WT$log2FoldChange < 1, ])

ego_up <- enrichGO(
  gene = WT_up,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  readable = TRUE
)

a1 <- dotplot(ego_up, showCategory = 8) +
  ggtitle("BP: WT infected vs control (upregulated genes)")

barplot(ego_up, showCategory = 8)+
  ggtitle("GO BP enrichment: WT infected vs control (upregulated genes)")


ego_down <- enrichGO(
  gene = WT_down,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  readable = TRUE
)

a2 <- dotplot(ego_down, showCategory = 8) +
  ggtitle("BP: WT infected vs control (downregulated genes)")

barplot(ego_down, showCategory = 8)


ego_up_MF <- enrichGO(
  gene = WT_up,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  readable = TRUE
)

a3 <- dotplot(ego_up_MF, showCategory = 8) +
  ggtitle("MF: WT infected vs control (upregulated genes)")

barplot(ego_up_MF, showCategory = 8)

a1 + a2 + a3

ego_down_MF <- enrichGO(
  gene = WT_down,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  readable = TRUE
)

dotplot(ego_down_MF, showCategory = 8) +
  ggtitle("MF: WT infected vs control (downregulated genes)")

barplot(ego_down_MF, showCategory = 8)


ego_up_CC <- enrichGO(
  gene = WT_up,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  readable = TRUE
)

dotplot(ego_up_CC, showCategory = 8) +
  ggtitle("CC: WT infected vs control (upregulated genes)")

barplot(ego_up_CC, showCategory = 8)

ego_down_CC <- enrichGO(
  gene = WT_down,
  universe = sig_WT,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  readable = TRUE
)

dotplot(ego_down_CC, showCategory = 8) +
  ggtitle("CC: WT infected vs control (downregulated genes)")

barplot(ego_down_CC, showCategory = 8)

# ---------------------------------------------------------------- #
# DKO samples
sig_DKO <- as.data.frame(res_DKO) %>%
  filter(!is.na(padj), padj < 0.05)

DKO_up <- rownames(sig_DKO[sig_DKO$log2FoldChange > 1, ])
DKO_down <- rownames(sig_DKO[sig_DKO$log2FoldChange < 1, ])

ego_up_DKO <- enrichGO(
  gene = DKO_up,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  readable = TRUE
)

ego_down_DKO <- enrichGO(
  gene = DKO_down,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  readable = TRUE
)

dotplot(ego_up_DKO, showCategory = 8) +
  ggtitle("BP: DKO infected vs control (Upregulated genes)")

barplot(ego_up_DKO, showCategory = 8)


dotplot(ego_down_DKO, showCategory = 8) +
  ggtitle("BP: DKO infected vs control (downregulated genes)")

barplot(ego_down_DKO, showCategory = 8)


ego_up_DKO_MF <- enrichGO(
  gene = DKO_up,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  readable = TRUE
)

ego_down_DKO_MF <- enrichGO(
  gene = DKO_down,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  readable = TRUE
)

dotplot(ego_up_DKO_MF, showCategory = 8) +
  ggtitle("MF: DKO infected vs control (Upregulated genes)")

barplot(ego_up_DKO_MF, showCategory = 8)


dotplot(ego_down_DKO_MF, showCategory = 8) +
  ggtitle("MF: DKO infected vs control (downregulated genes)")

barplot(ego_down_DKO_MF, showCategory = 8)


ego_up_DKO_CC <- enrichGO(
  gene = DKO_up,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  readable = TRUE
)

ego_down_DKO_CC <- enrichGO(
  gene = DKO_down,
  universe = sig_DKO,
  OrgDb = org.Mm.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  readable = TRUE
)

dotplot(ego_up_DKO_CC, showCategory = 8) +
  ggtitle("CC: DKO infected vs control (Upregulated genes)")

barplot(ego_up_DKO_CC, showCategory = 8)


dotplot(ego_down_DKO_CC, showCategory = 8) +
  ggtitle("CC: DKO infected vs control (downregulated genes)")

barplot(ego_down_DKO_CC, showCategory = 8)
