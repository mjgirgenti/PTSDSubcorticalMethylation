args = (commandArgs(TRUE))
CHR = as.numeric(args[1])
pheno = args[2]
region = args[3]

cat("region =", region, ", disorder =", pheno, ", chr.num =", CHR)

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(quadprog)
library(locfdr)
library(EpiDISH)
library(data.table)
library(plyr)
library(dplyr)
library(purrr)
library(GenomicRanges)
library(Repitools)

data_dir = "data_dir"
load(paste0(data_dir, "metadata_06052020.RData"))
load(paste0(data_dir, "data.depth.normalized.snpFiltered.invariantFiltered.gr.RData"))
load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.gr.RData"))
load(paste0(data_dir, "data.annotation.snpFiltered.gr.RData"))
load(paste0(data_dir, "data.genotype.ancestry.proportion.RData"))

rna.metadata = read.csv(paste0(data_dir, "MetadataForRNASeq.csv"))
AMG = read.csv(paste0(data_dir, "CIBERSORTx_Job2_Results (AMG_1000per).csv"))
HPC = read.csv(paste0(data_dir, "CIBERSORTx_Job4_Results (HPC_1000per).csv"))
celltype = rbind(AMG, HPC)

tmp.df = merge(celltype,rna.metadata[,c("LibName","Brain","Region")], by.x = "Mixture", by.y = "LibName")
tmp.df$Region[tmp.df$Region == "S"] = "Sub"
tmp.df$Region[tmp.df$Region == "CeA"] = "CNA"
tmp.df$Region[tmp.df$Region == "MeA"] = "MNA"

sample.id.df = merge(tmp.df, metadata.df[,c("ID", "SampleID", "RegionSci")], by.x = c("Brain","Region"), by.y = c("ID","RegionSci"))
rownames(sample.id.df) = sample.id.df$SampleID
CTS_Proportion = sample.id.df[,colnames(sample.id.df)[4:9]]

region_disorder <- paste(region, pheno, sep = "_")
region_disorder_chr = paste(region, pheno, paste0("chr", CHR), sep = "_")

######################################################################################################
################### ------------------------ PREPARE DATA ------------------------ ###################
######################################################################################################

# -------------------------------------------------------------------------------------------------- #
#                                         methylation data                                           #
# -------------------------------------------------------------------------------------------------- #

data.ratio.df <- annoGR2DF(data.ratio.gr)
data.depth.df <- annoGR2DF(data.depth.gr)
data.ratio.df <- data.ratio.df[,match(names(data.depth.df), names(data.ratio.df))]

# samples
metadata.df <- metadata.df[!is.na(metadata.df$Smoking),]
sampleSets <- list(BLA_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "BLA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   CNA_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "CNA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   CA_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "CA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   MNA_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "MNA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   DG_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "DG" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   Sub_MDD=unlist(metadata.df[which(metadata.df$RegionSci == "Sub" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="MDD")), "SampleID"]),
                   BLA_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "BLA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]),
                   CNA_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "CNA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]),
                   CA_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "CA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]),
                   MNA_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "MNA" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]),
                   DG_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "DG" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]),
                   Sub_PTSD=unlist(metadata.df[which(metadata.df$RegionSci == "Sub" & (metadata.df$DxGroup1=="Control" | metadata.df$DxGroup1=="PTSD")), "SampleID"]))

samples <- sampleSets[[region_disorder]]
names(samples) <- samples

# -------------------------------------------------------------------------------------------------- #
#                             Split data into data frame for each sample                             #
# -------------------------------------------------------------------------------------------------- #

# Extract data for each chr.num
data.ratio.df <- data.ratio.df[data.ratio.df$chr == CHR,]
data.depth.df <- data.depth.df[data.depth.df$chr == CHR,]

beta.df <- data.frame(chr=data.ratio.df$chr,
                      pos=data.ratio.df$start)
frames <- list()
for (sample in names(samples)) {
  print(sample)
  if (!is.null(data.ratio.df[[sample]])) {
    tab <- data.frame(chr=data.ratio.df$chr,
                      pos=data.ratio.df$start,
                      N=data.depth.df[[sample]],
                      X=data.depth.df[[sample]]*data.ratio.df[[sample]])
    frames[[sample]] <- tab
    beta.df <- cbind(beta.df, sample=data.ratio.df[[sample]])
    setnames(beta.df, "sample", sample)
  }
}

rownames(beta.df) = paste("chr", beta.df$chr, "_", beta.df$pos, sep="")
beta.df = beta.df[,-c(1:2)]
beta.df = na.omit(beta.df)

# kick out any samples not found in actual data but in metadata
samples <- samples[samples %in% names(frames)]

# prep metadata

metadata.df = cbind(metadata.df, data.genotype.ancestry.proportion[match(metadata.df$ID, rownames(data.genotype.ancestry.proportion)),])
metadata.analysis <- metadata.df[metadata.df$SampleID %in% samples,]
metadata.analysis <- metadata.analysis[, c("SampleID","DxGroup1", "Age", "Sex", "PMI", "Source", "Smoking", "AFR", "EUR", "SAS")]

beta.df = beta.df[,intersect(colnames(beta.df) ,rownames(CTS_Proportion))]
frac.m = CTS_Proportion[intersect(colnames(beta.df) ,rownames(CTS_Proportion)),]
rownames(metadata.analysis) = metadata.analysis$SampleID
metadata.analysis = metadata.analysis[intersect(rownames(metadata.analysis),colnames(beta.df)),]

######################################################################################################
################### -------------------------  CellDMC ------------------------- ###################
######################################################################################################

results = CellDMC(as.matrix(beta.df), 
                  metadata.analysis$DxGroup1, 
                  as.matrix(frac.m),
                  cov.mod = model.matrix(~ Age + Sex + PMI  + Smoking + AFR + EUR, data = metadata.analysis))

# write results to Rdata file ----

region_disorder_chr = paste(region, "_", pheno, "_chr", CHR, sep = "")
res_dir = "res_dir"

dmct.outputs = results$dmct
Micro.outputs = results$coe$Micro
Oligo.outputs = results$coe$Oligo
OPC.outputs = results$coe$OPC
Astro.outputs = results$coe$Astro
Inhib.outputs = results$coe$Inhib
Excit.outputs = results$coe$Excit

save(Micro.outputs, file = paste0(res_dir, "Micro_Results_", region_disorder_chr, ".Rdata"))
save(Oligo.outputs, file = paste0(res_dir, "Oligo_Results_", region_disorder_chr, ".Rdata"))
save(OPC.outputs, file = paste0(res_dir, "OPC_Results_", region_disorder_chr, ".Rdata"))
save(Astro.outputs, file = paste0(res_dir, "Astro_Results_", region_disorder_chr, ".Rdata"))
save(Inhib.outputs, file = paste0(res_dir, "Inhib_Results_", region_disorder_chr, ".Rdata"))
save(Excit.outputs, file = paste0(res_dir, "Excit_Results_", region_disorder_chr, ".Rdata"))
