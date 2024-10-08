args = (commandArgs(TRUE))
region = args[[1]]
disorder = args[[2]]
chr.num = args[[3]]
type = args[[4]]

cat("region =", region, ", disorder =", disorder, ", chr.num =", chr.num, ", type =", type)

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(data.table)
library(DSS)
library(plyr)
library(dplyr)
library(purrr)
library(GenomicRanges)
library(Repitools)

# -------------------------------------------------------------------------------------------------- #
#                                                User Inputs                                         #
# -------------------------------------------------------------------------------------------------- #

data_dir = "data_dir"
load(paste0(data_dir, "metadata_06052020.RData"))
load(paste0(data_dir, "data.annotation.snpFiltered.gr.RData"))
load(paste0(data_dir, "MethylationCellTypeProportions.RData"))
load(paste0(data_dir, "data.genotype.ancestry.proportion.RData"))

if (type == "CHG") {
  load(paste0(data_dir, "data.depth.normalized.snpFiltered.invariantFiltered.overlapped.CHG.gr.RData"))
  load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHG.gr.RData"))
} else if (type == "CHH") {
  load(paste0(data_dir, "data.depth.normalized.snpFiltered.invariantFiltered.overlapped.CHH.gr.RData"))
  load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHH.gr.RData"))
}

region_disorder <- paste(region, disorder, sep = "_")
region_disorder_chr_type = paste(region, disorder, paste0("chr", chr.num), type, sep = "_")

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
data.ratio.df <- data.ratio.df[data.ratio.df$chr == chr.num,]
data.depth.df <- data.depth.df[data.depth.df$chr == chr.num,]

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

# kick out any samples not found in actual data but in metadata
samples <- samples[samples %in% names(frames)]

# making bsseq object
bsdata <- makeBSseqData(frames, samples)

# prep metadata
metadata.df$NeuPos = CTS_Proportion[match(metadata.df$SampleID, rownames(CTS_Proportion)),]$`Neu+`
metadata.df = cbind(metadata.df, data.genotype.ancestry.proportion[match(metadata.df$ID, rownames(data.genotype.ancestry.proportion)),])
metadata.analysis <- metadata.df[metadata.df$SampleID %in% samples,]
metadata.analysis <- metadata.analysis[, c("DxGroup1", "Age", "Sex", "PMI", "Source", "Smoking", "NeuPos", "AFR", "EUR", "SAS")]

######################################################################################################
################### -------------------------  FIND DMRS ------------------------- ###################
######################################################################################################

# fit glm for differentially methylated loci (dml)
dml.fit <- DMLfit.multiFactor(bsdata, metadata.analysis,
                              formula = as.formula("~DxGroup1 + Age + Sex + PMI + Source + Smoking + NeuPos + AFR + EUR"),
                              smoothing = FALSE)
# extract dml results
coef_level <- unlist(tstrsplit(region_disorder, "_", keep = 2))
coef <- paste0("DxGroup1", coef_level)

dml.test <- DMLtest.multiFactor(dml.fit, coef=coef)
dml.results <- merge(dml.test, beta.df, by=c("chr", "pos"), all.x=TRUE)
setorder(dml.results, pvals)

res_dir = "res_dir"
save(dml.results, file = paste0(res_dir, "DSS_Filtered_NoSmooth_Ancestry_Prop_Results_", region_disorder_chr_type, ".Rdata"))

