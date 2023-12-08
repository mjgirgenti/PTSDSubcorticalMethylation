
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

work_dir = "work_dir"

data_dir = "data_dir"
load(glue({data_dir}, "data.ratio.normalized.snpFiltered.invariantFiltered.gr.RData"))
load(glue({data_dir}, "metadata_06052020.RData"))

# Percent methylation data frame
data.ratio.df <- annoGR2DF(data.ratio.gr)
rownames(data.ratio.df) <- paste(data.ratio.df$chr, data.ratio.df$start, sep = "_")
data.ratio.df <- data.ratio.df[,-which(names(data.ratio.df) %in% c("chr", "start", "end", "width"))]

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
results.pca <- pca(data, nPcs = 10,  scale = "uv", center = TRUE)

# Export  results
save(results.pca, file = glue({work_dir}, "PCA_", group, ".results.RData"))
