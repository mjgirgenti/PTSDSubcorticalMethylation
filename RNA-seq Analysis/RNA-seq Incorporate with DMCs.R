library(dplyr)
library(biomaRt)
library(data.table)
library(ggplot2)

# # Gene Annotation
# ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host = "https://apr2019.archive.ensembl.org") # version 96

gene_annotation <- getBM(attributes=c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'gene_biotype',
                   'start_position', 'end_position', 'strand', 'transcription_start_site'), 
                         filters = 'ensembl_gene_id',
                         values = rownames(fc_geneCount$counts),
                         mart = ensembl)
## CpG -----
# Load differential methylation result
met_res_PTSD <- get(load("../FinalResults/PTSD_DSS_Filtered_NoSmooth_Ancestry_Prop_FinalResults.RData"))
met_res_MDD <- get(load("../FinalResults/MDD_DSS_Filtered_NoSmooth_Ancestry_Prop_FinalResults.RData"))

# Load differential expression result
deg_result <- readRDS("DEG_results.rds")
names(deg_result$Control_PTSD) <- c("Sub", "CA", "DG", "BLA", "MNA", "CNA")
names(deg_result$Control_MDD) <- c("Sub", "CA", "DG", "BLA", "MNA", "CNA")

# Combine methylation and expression results
combine_met_exp <- function(met_res, exp_res) {
  methyl_stat_df <- data.frame(chr_pos = rep(met_res$chr_pos, 6),
                               nearestGeneName = rep(met_res$nearestGeneName, 6),
                               distanceToNearestGene = rep(met_res$distanceToNearestGene, 6),
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
    merge(., geneAnnot, by = "GeneID")
  combined_stat_df <- merge(methyl_stat_df, exp_stat_df, 
                            by.x = c("nearestGeneName", "region"), 
                            by.y = c("GeneSymbol", "region")) 
  return(combined_stat_df)
}

combined_stat_df_PTSD <- combine_met_exp(met_res_PTSD, exp_res = deg_result$Control_PTSD)
combined_stat_df_MDD <- combine_met_exp(met_res_MDD, exp_res = deg_result$Control_MDD)

keep_one_dmc_df_PTSD <- combined_stat_df_PTSD %>%
  subset(distanceToNearestGene < 5000) %>%
  group_by(region, nearestGeneName) %>%
  slice_min(order_by = methyl_pval, n = 1) %>% ungroup()
keep_one_dmc_df_MDD <- combined_stat_df_MDD %>%
  subset(distanceToNearestGene < 5000) %>%
  group_by(region, nearestGeneName) %>%
  slice_min(order_by = methyl_pval, n = 1) %>% ungroup()
combined_stat_df_PTSD %>%
  subset(distanceToNearestGene < 5000) %>% 
  subset(region == "MNA") %>%
  dim()
keep_one_dmc_df_PTSD %>%
  #subset(region == "CNA") %>%
  subset(methyl_pval < 0.005 & exp_pval < 0.05) %>%
  subset(!duplicated(GeneID)) %>%
  #subset(duplicated(GeneID)| duplicated(GeneID, fromLast = T)) %>% 
  mutate(methyl_direction = sign(methyl_stat), exp_direction = sign(exp_stat)) %>% 
  dplyr::select(methyl_direction, exp_direction) %>%
  table()
relocate(methyl_direction, exp_direction) %>% View()
dim()
mutate(methyl_direction = sign(methyl_stat), exp_direction = sign(exp_stat)) %>%
  dplyr::select(methyl_direction, exp_direction) %>%
  table()
deg_PTSD %>%
  merge(geneAnnot, by = "GeneID") %>% 
  subset(GeneSymbol == "GAD2")
# Enrichment
enrichment_table <- list()
for (reg in c("BLA", "CNA", "MNA", "CA", "DG", "Sub")) {
  # Genes near DMCs
  m <- keep_one_dmc_df_PTSD %>% subset(region == reg & methyl_pval < 0.005) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  # DEG
  n <- keep_one_dmc_df_PTSD %>% subset(region == reg & exp_pval < 0.05) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  # Overlapped genes
  x <- m[m %in% n]
  # Considered genes
  k <- keep_one_dmc_df_PTSD %>% subset(region == reg) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  
  M <- matrix(c(length(x), length(m) -length(x), 
                length(n) -length(x), length(k) - length(m) -length(n) + length(x)), nrow = 2)
  chisq_obj <- chisq.test(M)
  fisher_obj <- fisher.test(M)
  
  enrichment_table[[reg]] <- c(chisq_obj$p.value, fisher_obj$p.value, length(x), length(m), length(n), length(k))
}

enrichment_table <- data.frame(enrichment_table, 
                               row.names = c("Chisq_pval", "Fisher_pval", "# overlapped", "# genes near DMC", "# DEGs", "# considered genes")) %>%
  t() %>% as.data.frame() %>% rownames_to_column("region") %>%
  mutate(region = factor(region, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub"),
                         labels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")))

enrichment_table$Chisq_pval_corrected <- p.adjust(enrichment_table$Chisq_pval, method = "bonferroni")

# CpH ----
CpH_result <- list()
for (file in list.files("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/PTSD_Methylation/DSS_Filtered_NoSmooth_Ancestry_Prop_CpH/results/")) {
  reg <- strsplit(file, "_")[[1]][7]
  dis <- strsplit(file, "_")[[1]][8]
  CH_type <- strsplit(file, "_")[[1]][10]
  load(paste0("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/PTSD_Methylation/DSS_Filtered_NoSmooth_Ancestry_Prop_CpH/results/", file))
  CpH_result <- rbind(CpH_result, cbind(dml.results[,1:5], region = reg, Dx = dis, CHType = CH_type))
}

load("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/CpH_annotation/CHG_gene_annotation.RData")
load("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/CpH_annotation/CHH_gene_annotation.RData")

met_CpH_res <- CpH_result %>%
  mutate(seqnames = paste0("chr", chr)) %>%
  merge(rbind(CHG_results, CHH_results), by.x = c("seqnames", "pos"), by.y = c("seqnames", "start")) %>%
  dplyr::select(chr, pos, stat, pvals, fdrs, region.x, CHType, name, description, distance)

keep_one_CpH <- rbindlist(deg_result$Control_PTSD, idcol = "region.x") %>%
  merge(geneAnnot, by = "GeneID") %>%
  merge(., met_CpH_res, by.x = c("region.x", "GeneSymbol"), by.y = c("region.x", "name")) %>%
  subset(description %in% c("inside intron", "inside exon") | distance < 5000) %>%
  group_by(region.x, GeneSymbol) %>%
  slice_min(order_by = pvals, n = 1) %>% 
  ungroup()

# Enrichment
CpH_enrichment_table <- list()
for (reg in c("BLA", "CNA", "MNA", "CA", "DG", "Sub")) {
  # for (type in c("CHG.Rdata", "CHH.Rdata")) {
  
  # Genes near DMCs
  m <- keep_one_CpH %>% subset(region.x == reg & pvals < 0.005) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  # DEG
  n <- keep_one_CpH %>% subset(region.x == reg & pvalue < 0.05) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  # Overlapped genes
  x <- m[m %in% n]
  # Considered genes
  k <- keep_one_CpH %>% subset(region.x == reg) %>% 
    dplyr::select(GeneID) %>% unlist(use.names = F) %>% unique()
  
  M <- matrix(c(length(x), length(m) -length(x), 
                length(n) -length(x), length(k) - length(m) -length(n) + length(x)), nrow = 2)
  chisq_obj <- chisq.test(M)
  fisher_obj <- fisher.test(M)
  
  CpH_enrichment_table[[reg]] <- c(chisq_obj$p.value, fisher_obj$p.value, length(x), length(m), length(n), length(k))
#}
}

CH_enrichment <- data.frame(CpH_enrichment_table,
           row.names = c("Chisq_pval", "Fisher_pval", "# overlapped", "# genes near DMC", "# DEGs", "# considered genes")) %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column("region") %>%
  mutate(Chisq_pval_corrected = p.adjust(Chisq_pval, method = "bonferroni"),
         region = factor(region, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub"),
                         labels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")))

CN_enrichment_table <- enrichment_table %>%
  mutate(region = factor(region, levels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")),
         type = "CG") %>%
  rbind(CH_enrichment %>% mutate(type = "CH"))

pdf("DEG_enrichment_CpN.pdf", height = 10)
keep_one_CpH %>%
  mutate(methyl_stat = stat.y,
         methyl_pval = pvals,
         exp_stat = stat.x,
         exp_pval = pvalue, 
         region = ifelse(Region == "CeA", "CNA", 
                         ifelse(Region == "MeA", "MNA", 
                                ifelse(Region == "S", "Sub", Region))),
         nearestGeneName = GeneSymbol) %>%
  dplyr::select(region, methyl_stat, methyl_pval, exp_stat, exp_pval, nearestGeneName) %>%
  mutate(type = "CH") %>%
  rbind(., keep_one_dmc_df_PTSD %>%
          dplyr::select(region, methyl_stat, methyl_pval, exp_stat, exp_pval, nearestGeneName) %>%
          mutate(type = "CG")) %>%
  mutate(region = factor(region, levels = c("BLA", "CNA", "MNA", "CA", "DG", "Sub"),
                         labels = c("BLA", "CeA", "MeA", "CA", "DG", "Sub")),
         col = ifelse(methyl_pval < 0.005, 
                      ifelse(exp_pval < 0.05, "Overlap", "DMC"),
                      ifelse(exp_pval < 0.05, "DEG", "NotSig")),
         alpha_val = ifelse(col == "Overlap", 1, 0.75)) %>%
  ggplot(., aes(x = -log10(methyl_pval)*sign(methyl_stat), y = -log10(exp_pval)*sign(exp_stat), fill = col)) + 
  ggrastr::rasterize(geom_point(shape = 21, aes(col = as.character(alpha_val), alpha = alpha_val, stroke = 0.25), show.legend = F), dpi = 300) +
  scale_fill_manual(values = c("blue", "red", "snow", "purple")) + 
  scale_color_manual(values = c("grey60", "grey15")) + 
  coord_fixed(ratio = 1) +
  facet_grid(region~type) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)), linetype = "dashed") +
  geom_vline(xintercept = c(-log10(0.005), log10(0.005)), linetype = "dashed") +
  #ggrepel::geom_text_repel(aes(label = nearestGeneName), show.legend = F) +
  geom_text(data = CN_enrichment_table, aes(x = 7.5, y = -5.5, label = paste("Enrichment p", 
                                                                           ifelse(Chisq_pval_corrected < 0.0001, "< 0.0001",
                                                                                  paste("=", signif(Chisq_pval_corrected, 3))))), inherit.aes = F, hjust = 1) +
  labs(x = expression("Methylation: -Log"[10]*"(p-value)"), 
       y = expression("Expression: -Log"[10]*"(p-value)")) +
  theme_bw() + theme(text = element_text(size = 12),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.x = element_blank())
dev.off()

## Export the table of genes
keep_one_CpH %>%
  mutate(methyl_stat = stat.y,
         methyl_pval = pvals,
         exp_stat = stat.x,
         exp_pval = pvalue, 
         region = ifelse(Region == "CeA", "CNA", 
                         ifelse(Region == "MeA", "MNA", 
                                ifelse(Region == "S", "Sub", Region))),
         nearestGeneName = GeneSymbol, 
         chr_pos = paste0(chr, "_", pos)) %>%
  dplyr::select(chr_pos, region, methyl_stat, methyl_pval, exp_stat, exp_pval, nearestGeneName) %>%
  mutate(type = "CH") %>%
  rbind(., keep_one_dmc_df_PTSD %>%
          dplyr::select(chr_pos, region, methyl_stat, methyl_pval, exp_stat, exp_pval, nearestGeneName) %>%
          mutate(type = "CG")) %>%
  subset(methyl_pval < 0.005 & exp_pval < 0.05) %>%
  write.csv("Signicant DMCs near DEGs.csv")
