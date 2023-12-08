library(dplyr)
library(tidyr)
library(ggplot2)
library(biomaRt)
library(ggpubr)

# Load differential methylation result
met_res_PTSD <- get(load("../FinalResults/PTSD_DSS_Filtered_NoSmooth_Ancestry_Prop_FinalResults.RData"))
met_res_MDD <- get(load("../FinalResults/MDD_DSS_Filtered_NoSmooth_Ancestry_Prop_FinalResults.RData"))

# Load differential expression result
deg_result <- readRDS("DEG_results.rds")
names(deg_result$Control_PTSD) <- c("Sub", "CA", "DG", "BLA", "MNA", "CNA")
names(deg_result$Control_MDD) <- c("Sub", "CA", "DG", "BLA", "MNA", "CNA")

ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://feb2023.archive.ensembl.org") # version 109
# listEnsemblArchives() 
# listAttributes(ensembl, what = "name") 
geneAnnot <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'gene_biotype'), 
                   filters = 'ensembl_gene_id', 
                   values = deg_result$Control_PTSD$S$GeneID, 
                   mart = ensembl)

# Combine methylation and expression results
combine_met_exp <- function(met_res, exp_res) {
  methyl_stat_df <- data.frame(chr_pos = rep(met_res$chr_pos, 6),
             nearestGeneName = rep(met_res$nearestGeneName, 6),
             region = rep(c("BLA", "MNA", "CNA", "Sub", "CA", "DG"), each = nrow(met_res)),
             methyl_stat = c(met_res$BLA_stat, met_res$MNA_stat, met_res$CNA_stat,
                             met_res$Sub_stat, met_res$CA_stat, met_res$DG_stat),
             methyl_pval = c(met_res$BLA_pval, met_res$MNA_pval, met_res$CNA_pval,
                             met_res$Sub_pval, met_res$CA_pval, met_res$DG_pval))
  # Reformat expression result
  exp_stat_df <- data.frame(GeneID = rep(exp_res$Sub$GeneID, 6),
                            region = rep(c("BLA", "MNA", "CNA", "Sub", "CA", "DG"), each = nrow(exp_res$Sub)),
                            exp_stat = c(exp_res$BLA$stat, exp_res$MNA$stat, exp_res$CNA$stat,
                                         exp_res$Sub$stat, exp_res$CA$stat, exp_res$DG$stat),
                            exp_pval = c(exp_res$BLA$pvalue, exp_res$MNA$pvalue, exp_res$CNA$pvalue,
                                         exp_res$Sub$pvalue, exp_res$CA$pvalue, exp_res$DG$pvalue)) %>%
    merge(., geneAnnot, by.x = "GeneID", by.y = "ensembl_gene_id")
  combined_stat_df <- merge(methyl_stat_df, exp_stat_df, by.x = c("nearestGeneName", "region"), by.y = c("external_gene_name", "region"))
  return(combined_stat_df)
}

combined_stat_df_PTSD <- combine_met_exp(met_res_PTSD, exp_res = deg_result$Control_PTSD)
combined_stat_df_MDD <- combine_met_exp(met_res_MDD, exp_res = deg_result$Control_MDD)

# Extract proportion of significant DE + DM that align
extract_stat_number <- function(combined_stat_df) {
  stat_number <- combined_stat_df %>%
    subset(methyl_pval < 0.05 & exp_pval < 0.05) %>%
    group_by(region, sign(methyl_stat), sign(exp_stat)) %>%
    summarise(count = n()) %>% ungroup() %>%
    group_by(region) %>%
    mutate(percentage = count/sum(count)) %>%
    as.data.frame()
  colnames(stat_number) <- c("region", "sign_methyl", "sign_exp", "count", "percentage")
  stat_number$region <- factor(stat_number$region, 
                               levels = c("BLA", "MNA", "CNA", "Sub", "CA", "DG"),
                               labels = c("BLA", "MeA", "CeA", "Sub", "CA", "DG"))
  stat_number <- stat_number %>%
    subset(sign_methyl*sign_exp < 0)
  return(stat_number)
}

stat_number_PTSD <- extract_stat_number(combined_stat_df_PTSD)
stat_number_MDD <- extract_stat_number(combined_stat_df_MDD)

# Get the table of DE+DM genes
significant_met_exp <- combined_stat_df_PTSD %>%
  subset(methyl_pval < 0.05 & exp_pval < 0.05) %>%
  mutate(met_signi = ifelse(methyl_stat > 0, "HyperMet","HypoMet"),
         exp_signi = ifelse(exp_stat > 0, "UpExp","DownExp")) %>%
  group_by(region, exp_signi, nearestGeneName) %>%
  summarise(MethylationState = ifelse(length(table(met_signi)) == 1, met_signi, "Mixed")) %>%
  ungroup() %>%
  mutate(Region = region, ExpressionState = exp_signi) %>%
  group_by(Region, ExpressionState, MethylationState) %>%
  summarise(NoGene= length(unique(nearestGeneName)),
            GeneList = gsub("\"", "",paste(list(unique(nearestGeneName))))) %>%
  mutate(GeneList = substr(GeneList,3,nchar(GeneList) -1)) %>%
  mutate(GeneList = gsub(",", ";", GeneList))

write.csv(significant_met_exp, "Signficant_Met_and_Exp.csv", row.names = F)

# Visualization
pdf("Results/Integration_PTSD_Ancestry_Prop_FinalResults.pdf", width = 8)
combined_stat_df_PTSD %>%
  subset(methyl_pval < 0.05) %>%
  mutate(region = factor(region, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub"),
                         labels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")),
         col = ifelse(methyl_pval < 0.05, ifelse(exp_pval < 0.05, ifelse(exp_stat > 0, "DEUp", "DEDown") ,"NotSig"), "NotSig")) %>%
  ggplot(., aes(x = -log10(methyl_pval)*sign(methyl_stat), y = -log10(exp_pval)*sign(exp_stat), fill = col)) + 
  facet_wrap(~region) +
  geom_point(shape = 21, col = "black", show.legend = F) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  scale_fill_manual(values = c("blue", "red", "lightgrey")) + coord_fixed(ratio = 1) + 
  labs(x = expression("Methylation: -Log"[10]*"(p-value)"), 
       y = expression("Expression: -Log"[10]*"(p-value)")) +
  geom_text(data = stat_number_PTSD, aes(x = 5.5*sign_methyl, y = 4.5*sign_exp, label = paste(round(percentage*100, 1), "%")), inherit.aes = F) +
  theme_bw() + theme(text = element_text(size = 12),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank())
dev.off()

pdf("Results/Integration_MDD_Ancestry_Prop_FinalResults.pdf", width = 8)
combined_stat_df_MDD %>%
  subset(methyl_pval < 0.05) %>%
  mutate(region = factor(region, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub"),
                         labels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")),
         col = ifelse(methyl_pval < 0.05, ifelse(exp_pval < 0.05, ifelse(exp_stat > 0, "DEUp", "DEDown") ,"NotSig"), "NotSig")) %>%
  ggplot(., aes(x = -log10(methyl_pval)*sign(methyl_stat), y = -log10(exp_pval)*sign(exp_stat), fill = col)) + 
  facet_wrap(~region) +
  geom_point(shape = 21, col = "black", show.legend = F) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  geom_vline(xintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  scale_fill_manual(values = c("blue", "red", "lightgrey")) + coord_fixed(ratio = 1) + 
  labs(x = expression("Methylation: -Log"[10]*"(p-value)"), 
       y = expression("Expression: -Log"[10]*"(p-value)")) +
  geom_text(data = stat_number_MDD, aes(x = 5.5*sign_methyl, y = 4.5*sign_exp, label = paste(round(percentage*100, 1), "%")), inherit.aes = F) +
  theme_bw() + theme(text = element_text(size = 12),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank())
dev.off()
