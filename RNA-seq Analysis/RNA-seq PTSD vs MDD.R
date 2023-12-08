library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)
library(ggpubr)
library(clusterProfiler)

# Get MDDvPTSD differential methylation result
res_anno_merged <- get(load("../../FinalResults/PTSDvMDD_DSS_Filtered_NoSmooth_Ancestry_Prop_PTSDvMDD_FinalResults.RData"))

# Get MDDvPTSD differential expression result
deg_result <- readRDS("../DEG_results.rds")
names(deg_result$MDD_PTSD) <- c("Sub", "CA", "DG", "BLA", "MNA", "CNA")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://feb2023.archive.ensembl.org") # version 109
# listEnsemblArchives() 
# listAttributes(ensembl, what = "name") 
geneAnnot <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'gene_biotype'), 
                   filters = 'ensembl_gene_id', 
                   values = deg_result$MDD_PTSD$S$GeneID, 
                   mart = ensembl)

# Combine methylation and expression stat
methyl_stat_df <- data.frame(chr_pos = rep(res_anno_merged$chr_pos, 6),
                             nearestGeneName = rep(res_anno_merged$nearestGeneName, 6),
                             region = rep(c("BLA", "MNA", "CNA", "Sub", "CA", "DG"), each = nrow(res_anno_merged)),
                             methyl_stat = c(res_anno_merged$BLA_stat, res_anno_merged$MNA_stat, res_anno_merged$CNA_stat,
                                             res_anno_merged$Sub_stat, res_anno_merged$CA_stat, res_anno_merged$DG_stat),
                             methyl_pval = c(res_anno_merged$BLA_pval, res_anno_merged$MNA_pval, res_anno_merged$CNA_pval,
                                             res_anno_merged$Sub_pval, res_anno_merged$CA_pval, res_anno_merged$DG_pval))

exp_stat_df <- data.frame(GeneID = rep(deg_result$MDD_PTSD$Sub$GeneID, 6),
                          region = rep(c("BLA", "MNA", "CNA", "Sub", "CA", "DG"), each = nrow(deg_result$MDD_PTSD$Sub)),
                          exp_stat = c(deg_result$MDD_PTSD$BLA$stat, deg_result$MDD_PTSD$MNA$stat, deg_result$MDD_PTSD$CNA$stat,
                                       deg_result$MDD_PTSD$Sub$stat, deg_result$MDD_PTSD$CA$stat, deg_result$MDD_PTSD$DG$stat),
                          exp_pval = c(deg_result$MDD_PTSD$BLA$pvalue, deg_result$MDD_PTSD$MNA$pvalue, deg_result$MDD_PTSD$CNA$pvalue,
                                       deg_result$MDD_PTSD$Sub$pvalue, deg_result$MDD_PTSD$CA$pvalue, deg_result$MDD_PTSD$DG$pvalue)) %>%
  merge(., geneAnnot, by.x = "GeneID", by.y = "ensembl_gene_id")
combined_stat_df <- merge(methyl_stat_df, exp_stat_df, by.x = c("nearestGeneName", "region"), by.y = c("external_gene_name", "region"))

# Extract DE and DM genes in agreement
upDE_hypoMet <- combined_stat_df %>%
  subset(methyl_pval < 0.05 & exp_pval < 0.05) %>%
  subset(methyl_stat < 0 & exp_stat > 0)
upDE_hypoMet_gene <- split(upDE_hypoMet$GeneID, f = upDE_hypoMet$region) %>%
  lapply(., unique)

downDE_hyperMet <- combined_stat_df %>%
  subset(methyl_pval < 0.05 & exp_pval < 0.05) %>%
  subset(methyl_stat > 0 & exp_stat < 0)
downDE_hyperMet_gene <- split(downDE_hyperMet$GeneID, f = downDE_hyperMet$region) %>%
  lapply(., unique)

# PTSD_dom <- compareCluster(upDE_hypoMet_gene, fun = "enrichGO", 
#                            universe = unique(combined_stat_df$GeneID),
#                            OrgDb = 'org.Hs.eg.db', 
#                            keyType = "ENSEMBL", ont = "BP", 
#                            pvalueCutoff = 0.1)
# 
# MDD_dom <- compareCluster(downDE_hyperMet_gene, fun = "enrichGO", 
#                           universe = unique(combined_stat_df$GeneID),
#                           OrgDb = 'org.Hs.eg.db', 
#                           keyType = "ENSEMBL", ont = "BP", 
#                           pvalueCutoff = 0.1)
# dotplot(PTSD_dom) # Preview to pick out pathways
# dotplot(MDD_dom) # Preview to pick out pathways

# Enrichment for PTSD > MDD and PTSD < MDD genes
PTSD_dom_full <- compareCluster(upDE_hypoMet_gene, fun = "enrichGO", 
                           universe = unique(combined_stat_df$GeneID),
                           OrgDb = 'org.Hs.eg.db', 
                           keyType = "ENSEMBL", ont = "BP", 
                           pvalueCutoff = 1, qvalueCutoff = 1)
MDD_dom_full <- compareCluster(downDE_hyperMet_gene, fun = "enrichGO", 
                          universe = unique(combined_stat_df$GeneID),
                          OrgDb = 'org.Hs.eg.db', 
                          keyType = "ENSEMBL", ont = "BP", 
                          pvalueCutoff = 1, qvalueCutoff = 1)

# Visualization
pdf("MDDvPTSD_Enrichment_heatmap.pdf", width =  8, height = 4)
PTSD_dom_full@compareClusterResult %>%
  subset(Description %in% c("modulation of chemical synaptic transmission",
                            "regulation of trans-synaptic signaling",
                            "learning",
                            "signal release",
                            "regulation of neuronal synaptic plasticity",
                            "epithelial cell development")) %>%
  mutate(Cluster = factor(Cluster, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub")),
         significant_sign = ifelse(p.adjust < 0.1, ifelse(p.adjust < 0.01, "**","*"), NA)) %>%
  ggplot(., aes(x = Cluster, y = Description, fill = -log10(p.adjust), label = significant_sign)) + 
  geom_tile() +
  geom_text() +
  geom_tile(data = data.frame(X = rep(1:6, 6), Y = rep(1:6, each = 6)), aes(x = X, y = Y), fill = NA, color = "black", inherit.aes = F) +
  scale_fill_gradient(low = "lavenderblush", high = "red") +
  labs(x = "", y = "", fill = expression("-Log"[10]*"(p-adjusted)"),
       title = "PTSD > MDD") + coord_fixed() +
  theme_bw() + theme(text = element_text(size = 14), panel.grid = element_blank(),
                     plot.title = element_text(hjust = 0.5), 
                     legend.position = "bottom") 

MDD_dom_full@compareClusterResult %>%
  subset(Description %in% c("chronic inflammatory response",
                            "myeloid leukocyte activation",
                            "regulation of leukocyte activation",
                            "response to lipopolysaccharide",
                            "regulation of epithelial cell proliferation",
                            "negative regulation of leukocyte differentiation")) %>%
  mutate(Cluster = factor(Cluster, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub")),
         significant_sign = ifelse(p.adjust < 0.1, ifelse(p.adjust < 0.01, "**","*"), NA)) %>%
  ggplot(., aes(x = Cluster, y = Description, fill = -log10(p.adjust), label = significant_sign)) + 
  geom_tile() +
  geom_text() +
  geom_tile(data = data.frame(X = rep(1:6, 6), Y = rep(1:6, each = 6)), aes(x = X, y = Y), fill = NA, color = "black", inherit.aes = F) +
  scale_fill_gradient(low = "mintcream", high = "blue") +
  labs(x = "", y = "", fill = expression("-Log"[10]*"(p-adjusted)"),
       title = "PTSD < MDD") + coord_fixed() +
  theme_bw() + theme(text = element_text(size = 14), panel.grid = element_blank(),
                     plot.title = element_text(hjust = 0.5), 
                     legend.position = "bottom")
dev.off()
