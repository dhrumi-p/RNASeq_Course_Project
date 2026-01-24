#!/usr/bin/Rscript


# RNA-Sequencing class, DEseq2 (differential gene expression) (BLOOD samples) 
# The preprocessing of samples was done including: quality control, trmming and alignment
# The counts were obtanined which were then further used for the downstream analysis.


# ---------------------------------------------------------------------------- #
#  libraries
# ---------------------------------------------------------------------------- #
#Used For differential expression analysis
library(DESeq2)
library(tidyverse)

#Used For Plotting and Visualization
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(EnhancedVolcano)

#Used For stats
library(vsn)
library(matrixStats)

#Used For data structuring, annotations and analysing
library(readr)
library(tidyverse)
library(dplyr)
library(tibble)
library(org.Mm.eg.db)
library(patchwork)
library(clusterProfiler)
library(AnnotationDbi)

# ---------------------------------------------------------------------------- #
#  counts and metadata
# ---------------------------------------------------------------------------- #

counts <- read_csv("gene_counts_trimmed.csv")

counts <- counts %>%
  column_to_rownames(var = "Geneid")

head(counts)
dim(counts) # 78334    15

# ---------------------------------------------------------------------------- #
meta_data <- read_csv("metaData_blooodSamples.csv")

meta_data <- meta_data %>%
  select(2,3) %>%
  column_to_rownames("Sample_id")

head(meta_data)

counts <- counts[, rownames(meta_data)]


meta_data <- meta_data %>%
  mutate(Genotype = ifelse(grepl("WT", Sample_type), "WT", "DKO"),
         Infection = ifelse(grepl("Case", Sample_type), "Infected", "Control")) # factors according to condition and genotype

meta_data$Genotype  <- factor(meta_data$Genotype,  levels = c("WT", "DKO"))
meta_data$Infection <- factor(meta_data$Infection, levels = c("Control", "Infected"))


meta_data$Sample_type <- factor(meta_data$Sample_type,
                                levels = c("Blood_WT_Control","Blood_DKO_Control",
                                           "Blood_WT_Case","Blood_DKO_Case")) # Convert Sample_type to different levels of conditions used in experiment

meta_data

# ---------------------------------------------------------------------------- #
#  Creating the DESeq2 dataset
# ---------------------------------------------------------------------------- #

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = meta_data,
  design    = ~ Genotype + Infection + Genotype:Infection
)

# To Filter lowly expressed genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

dds
# dim: 30031 15 

# ---------------------------------------------------------------------------- #
#  Runing the DESeq2
# ---------------------------------------------------------------------------- #

dds <- DESeq(dds)
resultsNames(dds)

# ---------------------------------------------------------------------------- #
#  Variance stabilizing transformation (for QC)
# ---------------------------------------------------------------------------- #

vsd <- vst(dds, blind = TRUE)

# ---------------------------------------------------------------------------- #
# Quality Control plots
# ---------------------------------------------------------------------------- #
# PCA plots
pcaData <- plotPCA(vsd, intgroup = c("Genotype","Infection"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Infection, shape = Genotype)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal() 

## Sample distance heatmap
sampleDists <- dist(t(assay(vsd)))

pheatmap(as.matrix(sampleDists),
         annotation_col = meta_data,
         main = "Sample-to-sample distances")

# ---------------------------------------------------------------------------- #
# Gene expression heatmaps
# ---------------------------------------------------------------------------- #
# Top 500 most variable genes
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)

pheatmap(assay(vsd)[topVarGenes, ],
         scale = "row",
         annotation_col = meta_data,
         show_rownames = FALSE,
         main = "Top 500 most variable genes")



# ---------------------------------------------------------------------------- #
# Differential expression contrasts
# ---------------------------------------------------------------------------- #
resultsNames(dds)
levels(dds$Genotype)
levels(dds$Infection)

## WT: Infected vs Control
res_WT <- results(dds,
                  name = "Infection_Infected_vs_Control")

res_WT <- res_WT[order(res_WT$padj), ]
summary(res_WT)


## DKO: Infected vs Control
res_DKO <- results(dds,
                   contrast = list(c("Infection_Infected_vs_Control",
                                     "GenotypeDKO.InfectionInfected")))

res_DKO <- res_DKO[order(res_DKO$padj), ]
summary(res_DKO)
# LFC > 0 (up)       : 4482, 15%
# LFC < 0 (down)     : 5158, 17%

## WT vs DKO under infected conditions (interaction-aware)
res_WTvsDKO_infected <- results(dds,contrast = list(c("Genotype_DKO_vs_WT",
                                                      "GenotypeDKO.InfectionInfected")))

res_WTvsDKO_infected <- res_WTvsDKO_infected[order(res_WTvsDKO_infected$padj), ]
summary(res_WTvsDKO_infected)
# LFC > 0 (up)       : 5758, 19%
# LFC < 0 (down)     : 4637, 15%

# WT vs DKO in controls
res_WTvsDKO_control <- results(dds,name = "Genotype_DKO_vs_WT")
summary(res_WTvsDKO_control)
# LFC > 0 (up)       : 96, 0.32%
# LFC < 0 (down)     : 296, 0.99%

res_interaction <- results(dds, name = "GenotypeDKO.InfectionInfected")
summary(res_interaction)
# LFC > 0 (up)       : 3681, 12%
# LFC < 0 (down)     : 2674, 8.9%

# ---------------------------------------------------------------------------- #
# Export DEG tables
# ---------------------------------------------------------------------------- #

write.csv(as.data.frame(res_WT),
          file = "DEG_WT_Infected_vs_Control.csv")

write.csv(as.data.frame(res_DKO),
          file = "DEG_DKO_Infected_vs_Control.csv")

write.csv(as.data.frame(res_WTvsDKO_infected),
          file = "DEG_WT_vs_DKO_Infected.csv")

write.csv(as.data.frame(res_WTvsDKO_control),
          file = "DEG_WT_vs_DKO_Control.csv")

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
# TOP 3 UP / DOWN GENE PLOTS (SIGNIFICANT ONLY)
get_top_genes <- function(res, padj_cutoff = 0.05, n = 3) {
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj), ]
  res <- res[res$padj < padj_cutoff, ]
  
  up <- res %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    head(n)
  
  down <- res %>%
    dplyr::filter(log2FoldChange < 0) %>%
    dplyr::arrange(log2FoldChange) %>%
    head(n)
  
  rbind(up, down)
}

# Barplot function for top genes
library(ggplot2)

plot_top_genes <- function(res, title) {
  top_genes <- get_top_genes(res)
  
  top_genes$gene <- rownames(top_genes)
  top_genes$direction <- ifelse(top_genes$log2FoldChange > 0, "Up", "Down")
  
  ggplot(top_genes,
         aes(x = reorder(gene, log2FoldChange),
             y = log2FoldChange,
             fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "firebrick", "Down" = "steelblue")) +
    theme_minimal(base_size = 12) +
    labs(
      title = title,
      x = NULL,
      y = "log2 Fold Change"
    )
}

p1 <- plot_top_genes(res_WT, "WT: Infected vs Control")
p2 <- plot_top_genes(res_DKO, "DKO: Infected vs Control")
p3 <- plot_top_genes(res_WTvsDKO_control, "WT vs DKO (Control)")
p4 <- plot_top_genes(res_WTvsDKO_infected, "WT vs DKO (Infected)")


(p1 | p2) /
  (p3 | p4)



# --------------------------------------------------------------- #
sessionInfo()
