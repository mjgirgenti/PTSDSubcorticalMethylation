library(dplyr)
library(biomaRt)
library(ggplot2)
library(reshape2)
library(DESeq2)
library(ggpubr)
library(tidyr)
library(ggrepel)
library(data.table)

## Metadata ----
demograph <- readxl::read_xlsx("~/lab/reference/NPBB.AllSubject Table.export20191003.xlsx") %>%
  dplyr::select(BrNum, AgeDeath, Sex, Race, Smoking, PMI, `Manner Of Death`, PrimaryDx)

metadata <- read.csv("Samples.csv") %>%
  merge(., demograph, by.x = "Brain", by.y = "BrNum")

rownames(metadata) <- metadata$LibName
metadata <- metadata %>%
  mutate(Region = factor(Region, levels = c("S", "CA", "DG", "BLA", "MeA", "CeA")),
         Brain = as.character(Brain))

## Filtering ----
# Gene Count
geneCount <- read.csv("GeneExpressionCount.csv", row.names = 1)
# Gene Annotation
ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://feb2023.archive.ensembl.org") # version 109

geneAnnot <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'gene_biotype'), 
                   filters = 'ensembl_gene_id', 
                   values = rownames(geneCount), 
                   mart = ensembl)
colnames(geneAnnot) <- c("GeneID", "Chr", "GeneSymbol", "GeneBiotype")
# Filter by count
filter_count <- function(cts, metadata, threshold = 1) { # mean/group > 1
  # Split batch into biological groups
  groups <- split(metadata$LibName, paste(metadata$Diagnosis, metadata$Region), sep ="")
  # Calculate row means
  for (group in groups) {
    cts %>%
      as.data.frame() %>%
      dplyr::select(all_of(group)) %>%
      base::rowMeans(.) > threshold -> x
    cts <- cts[x, ]
  }
  return(cts)
} # Mean count per region per diagnosis > threshold

geneCount_H <- geneCount %>%
  filter_count(., metadata, threshold = 1)

dds <- DESeqDataSetFromMatrix(geneCount_H[,metadata$LibName], 
                              colData = metadata,
                              design = ~ Batch + Sex + Diagnosis + Race)
normed_df <- vst(dds)
PCA_highCount <- ggarrange(plotPCA(normed_df, intgroup = "Batch") + theme_light() + labs(col = "Batch"),
                           plotPCA(normed_df, intgroup = "Region") + theme_light() + labs(col = "Region"),
                           plotPCA(normed_df, intgroup = "Parameter") + theme_light() + labs(col = "Parameter"),
                           plotPCA(normed_df, intgroup = "Sex") + theme_light() + labs(col = "Sex"))

## DEG ----
# Run DESeq2 for each region
region_list = c("S","CA","DG","BLA","MeA","CeA")

deg_results <- list()
for (region in region_list) {
  sub_meta <- subset(metadata, Region == region)
  sub_meta[,c("PMI", "AgeDeath", "RIN")] <- scale(sub_meta[,c("PMI", "AgeDeath", "RIN")], center = T)
  sub_meta$BatchPara = paste0(sub_meta$Batch, sub_meta$Parameter)
  sub_dds <- DESeqDataSetFromMatrix(countData = geneCount_H[,sub_meta$LibName],
                                    colData = sub_meta,
                                    design = ~ PMI + BatchPara + Sex + Race + AgeDeath + AgeDeath^2 + RIN + RIN^2 + Diagnosis)
  sub_dds <- DESeq(sub_dds)
  for (comparison in list(c("PTSD", "Control"), c("MDD", "Control"), c("PTSD", "MDD"))) {
    comparison_name <- paste0(comparison[2], "_", comparison[1])
    deg_results[[comparison_name]][[region]] <- results(sub_dds, contrast = c("Diagnosis", comparison[1], comparison[2])) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("GeneID") %>%
      mutate(Region = region, Comparison = comparison_name) # Easier for combining results later
  }
}

rm(sub_meta, sub_dds)
saveRDS(deg_results, "DEG_results.rds") # Save the output


## Figures ----
# Identify DEGs
pdf("Results/DEGs.pdf")
rbindlist(deg_results$Control_PTSD) %>%
  mutate(Region = factor(Region, levels = region_list),
         padj = replace_na(padj, 1)) %>%
  merge(., geneAnnot, by = "GeneID") %>%
  mutate(fdr_sig = padj < 0.05 & abs(log2FoldChange) > log2(1.2)) %>%
  mutate(name = ifelse(pvalue < 0.005, GeneSymbol, NA),
         color = ifelse(fdr_sig, ifelse(Region %in% c("S", "CA", "DG"), "DEhpc", "DEamg"), "NotDE")) %>%
  group_by(Region) %>%
  ggplot(., aes(x = log2FoldChange, y = -log10(pvalue), fill = color)) + 
  geom_hline(yintercept = -log10(0.001), linetype = "dashed") +
  geom_point(shape = 21, col = "black", size = 2, show.legend = F) + 
  facet_wrap(~Region) + 
  geom_text_repel(aes(label = name), show.legend = F) + 
  scale_fill_manual(values = c('#F8766D', '#619CFF', 'gray')) +
  labs(title = 'Using p-value',
       x = expression("Log"[2]*"(fold change)"),
       y = expression("-Log"[10]*"p-value")) +
  theme_bw()
# PCA
annotate_figure(PCA_highCount, top = text_grob("PCA highCount"))
dev.off()
