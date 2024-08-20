
##########################################################################################
## Load Packages
##########################################################################################

# Load Packages
library(glue)
library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(gridExtra)
library(purrr)
library(sjstats)
library(plyr)
library(rcartocolor)
library(viridis)
library(stringr)
library(ggrepel)
library(Repitools)
library(pcaMethods)

data_dir = "data_dir"
work_dir = "work_dir"

load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHG.gr.RData"))
data.ratio.CHG.df <- annoGR2DF(data.ratio.CHG.gr)

load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHH.gr.RData"))
data.ratio.CHH.df <- annoGR2DF(data.ratio.CHH.gr)

load(glue({data_dir}, "metadata_06052020.RData"))

# Percent methylation data frame

rownames(data.ratio.CHG.df) <- paste(data.ratio.CHG.df$chr, data.ratio.CHG.df$start, sep = "_")
data.ratio.CHG.df <- data.ratio.CHG.df[,-which(names(data.ratio.CHG.df) %in% c("chr", "start", "end", "width", "strand"))]

rownames(data.ratio.CHH.df) <- paste(data.ratio.CHH.df$chr, data.ratio.CHH.df$start, sep = "_")
data.ratio.CHH.df <- data.ratio.CHH.df[,-which(names(data.ratio.CHH.df) %in% c("chr", "start", "end", "width", "strand"))]

data.ratio.df = rbind(data.ratio.CHG.df[,intersect(colnames(data.ratio.CHG.df), colnames(data.ratio.CHH.df))],
                      data.ratio.CHH.df[,intersect(colnames(data.ratio.CHG.df), colnames(data.ratio.CHH.df))])

print("transposing matrix")
# Need to transpose for PCA
t.data.mat <- t(data.ratio.df)
t.data.df <- as.data.frame(t.data.mat)

# Sample sets
all.samples <- unlist(metadata.df[, "SampleID"])
control.samples <- unlist(metadata.df[metadata.df$DxGroup1 == "Control", "SampleID"])
ptsd.samples <- unlist(metadata.df[metadata.df$DxGroup1 == "PTSD", "SampleID"])
mdd.samples <- unlist(metadata.df[metadata.df$DxGroup1 == "MDD", "SampleID"])

group = "All"
sample.list <- list(All = all.samples,
                    Control = control.samples,
                    PTSD = ptsd.samples,
                    MDD = mdd.samples)
samples <- sample.list[[group]]

# Create datasets
dataList <- list()
dataList[[group]] <- t.data.mat[rownames(t.data.mat) %in% samples,]
data <- dataList[[group]]

# Run PCA
data = apply(data, 2, as.numeric)
results.pca <- pca(data, nPcs = 10, scale = "uv", center = TRUE)

# Export  results
save(results.pca, file = glue({work_dir}, "PCA_", "CpH", ".results.RData"))