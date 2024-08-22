
args = (commandArgs(TRUE))
print(args)

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(glue)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(purrr)
library(plyr)
library(rcartocolor)
library(viridis)
library(stringr)
library(ggrepel)
library(Repitools)
library(pcaMethods)
library(pheatmap)
library(Cairo)
library(ComplexHeatmap)
library(tidyverse)
library(writexl)
library(GenomicRanges)

######################################################################################################
################### -------------------------- Prepare Data ---------------------- ###################
######################################################################################################

data_dir = "data_dir"

# import DSS results files ----
results.list = list.files(path = glue({data_dir}, "dataDSS/NewAnnotation/"),
                          pattern = "NewAnnotation_DSS_", full.names = TRUE)
fn = list.files(path = glue({data_dir}, "dataDSS/NewAnnotation/"),
                pattern = "NewAnnotation_DSS_", full.names = FALSE)
fn = unlist(tstrsplit(fn, "_", keep = 3))
names(results.list) = fn

dss.results.list = list()
for (r in names(results.list)) {
  load(results.list[[r]])
  dss.results.list[[r]] = annoGR2DF(data_annotated)
}

# import celldmc results ----

CellDMC_res = read.csv(paaste0(data_dir, "CellDMC_res.csv"))

annotated.df = dss.results.list[[1]]
annotated.df = annotated.df[,1:11]

CellDMC_res = CellDMC_res[,c("chr","position","Estimate","t","p","celltype","region")]
CellDMC_res$direction = ifelse(CellDMC_res$Estimate > 0, "Hyper", "Hypo")
CellDMC_res$cpgID = paste0(CellDMC_res$chr, "_",CellDMC_res$position)
CellDMC_res$category = paste0(CellDMC_res$celltype,"-",CellDMC_res$direction)
CellDMC_res = merge(CellDMC_res, annotated.df[,c("cpgID","relativeLocation", "CpG", "Chromatin")], by = "cpgID")
CellDMC_res = CellDMC_res[!is.na(CellDMC_res$Chromatin),]
CellDMC_res$relativeLocation = droplevels(CellDMC_res$relativeLocation)

CellDMC_res$relativeLocation = factor(CellDMC_res$relativeLocation)
CellDMC_res$CpG = factor(CellDMC_res$CpG)
CellDMC_res$Chromatin = factor(CellDMC_res$Chromatin)

levels(CellDMC_res$relativeLocation) = unique(CellDMC_res$relativeLocation)
levels(CellDMC_res$CpG) = unique(CellDMC_res$CpG)
levels(CellDMC_res$Chromatin) = unique(CellDMC_res$Chromatin)

target.df = CellDMC_res %>%
  filter(category == args & p <= 0.05)

sig_joint_chr_pos = target.df$cpgID

##########################################################################################
## PTSD GWAS
##########################################################################################

ptsd_gwas_hit_1 = data.frame(snp = c("rs7519147", "rs2777888", "rs4697248", "rs7688962",
                                     "rs10235664", "rs67529088", "rs2532252", "rs2123392"),
                             chr = c(1, 3, 4, 4, 7, 7, 17, 18),
                             pos = c(73994416, 49898000, 21931195, 88281182,
                                     2086814, 104907066, 44257783, 53214865),
                             LD_start = NA,
                             LD_end = NA)

ptsd_gwas_2021 = read.table(paste0(data_dir, "PTSD_GWAS_Hits.txt"))
ptsd_gwas_hit_2 = data.frame(snp = ptsd_gwas_2021$V1,
                             chr = ptsd_gwas_2021$V2,
                             pos = ptsd_gwas_2021$V7,
                             LD_start = NA,
                             LD_end = NA)

ptsd_gwas_ebi = read.csv(paste0(data_dir, "efotraits_EFO_0001358-associations-2022-12-28.csv"))
ptsd_gwas_hit_3 = data.frame(study = ptsd_gwas_ebi$Study.accession,
                             snp = unlist(lapply(strsplit(ptsd_gwas_ebi$Variant.and.risk.allele, "-"), function (x) x[[1]])),
                             chr = as.numeric(unlist(lapply(strsplit(ptsd_gwas_ebi$Location, ":"), function (x) x[[1]]))),
                             pos = as.numeric(unlist(lapply(strsplit(ptsd_gwas_ebi$Location, ":"), function (x) x[[2]]))),
                             LD_start = NA,
                             LD_end = NA)

strsplit(ptsd_gwas_ebi$Mapped.gene, ",") %>% unlist %>% unique %>% sort -> selected.genes

selected.genes %>% length()
ptsd_gwas_hit_3$study %>% unique() %>% length()
ptsd_gwas_hit_3 %>% dim()

# ptsd_gwas_hit = ptsd_gwas_hit_3 
ptsd_gwas_hit = rbind(ptsd_gwas_hit_1, ptsd_gwas_hit_2, ptsd_gwas_hit_3[,-1])

bp_range = 500000
ptsd_gwas_hit$LD_start = ptsd_gwas_hit$pos - bp_range
ptsd_gwas_hit$LD_end = ptsd_gwas_hit$pos + bp_range

ptsd_gwas_hit$CHR = paste0("chr", ptsd_gwas_hit$chr)
ptsd_gwas_hit

celldmc.gwas <- GRanges(
  seqnames = annotated.df$chr,
  ranges = IRanges(start = annotated.df$start, end = annotated.df$end),
)

mcols(celldmc.gwas) <- DataFrame(
  nearestGeneName = annotated.df$nearestGeneName,
  distanceToNearestGene = annotated.df$distanceToNearestGene,
  relativeLocation = annotated.df$relativeLocation
)

res_list = list()
for (i in 1:dim(ptsd_gwas_hit)[1]) {
  res_list[[i]] = celldmc.gwas[(celldmc.gwas@seqnames == ptsd_gwas_hit[i, "CHR"]) &
                                 (celldmc.gwas@ranges@start >= ptsd_gwas_hit[i, "LD_start"]) &
                                 ((celldmc.gwas@ranges@start + 1) <= ptsd_gwas_hit[i, "LD_end"])]
}
lapply(res_list, function(dat) paste(dat@seqnames, dat@ranges@start, sep = "_")) %>% unlist %>% unique() -> cpgs.selected

res_all = celldmc.gwas
res_all$chr_pos = paste(res_all@seqnames, res_all@ranges@start, sep = "_")

res.in = res_all[res_all$chr_pos %in% cpgs.selected,]
res.out = res_all[!(res_all$chr_pos %in% cpgs.selected),]

######################################################################################################
################### ---------------------- Fisher's Exact Test ------------------- ###################
######################################################################################################

sum(sig_joint_chr_pos %in% res.in$chr_pos) / length(res.in)
sum(sig_joint_chr_pos %in% res.out$chr_pos) / length(res.out)

100*(sum(res.in$chr_pos %in% sig_joint_chr_pos) / sum(!res.in$chr_pos %in% sig_joint_chr_pos))
100*(sum(res.out$chr_pos %in% sig_joint_chr_pos) / sum(!res.out$chr_pos %in% sig_joint_chr_pos))

x <- matrix(c(sum(sig_joint_chr_pos %in% res.in$chr_pos),
              length(res.in),
              sum(sig_joint_chr_pos %in% res.out$chr_pos),
              length(res.out)), byrow = TRUE, 2, 2)

fisher.test(x)