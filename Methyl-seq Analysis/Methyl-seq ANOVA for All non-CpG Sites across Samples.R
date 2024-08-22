
args = (commandArgs(TRUE))
index = as.numeric(args[[1]])

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(data.table)
library(dplyr)
library(readxl)
library(writexl)
library(Repitools)

site.index = seq(1,354983, by = 3500)
first.site = site.index[index]

if (index == 102) {
  last.site = 354983
} else {
  last.site = site.index[index+1] - 1
}

data_dir = "data_dir"
load(paste0(data_dir,"metadata_06052020.RData"))
load(paste0(data_dir,"data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHG.gr.RData"))
load(paste0(data_dir,"data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHH.gr.RData"))

######################################################################################################
################### -------------------------- Prepare Data ---------------------- ###################
######################################################################################################

data.ratio.CHG.df = annoGR2DF(data.ratio.CHG.gr)
data.ratio.CHH.df = annoGR2DF(data.ratio.CHH.gr)

data.ratio.df = rbind(data.ratio.CHG.df[,intersect(colnames(data.ratio.CHG.df), colnames(data.ratio.CHH.df))],
                      data.ratio.CHH.df[,intersect(colnames(data.ratio.CHG.df), colnames(data.ratio.CHH.df))])
data.ratio.df$chr = paste0("chr",data.ratio.df$chr)
data.ratio.df$cpgID = paste(data.ratio.df$chr, data.ratio.df$start, sep = "_")

## create test dataframe----
anova.df = data.ratio.df
rownames(anova.df) = anova.df$cpgID
anova.df = anova.df[,7:ncol(anova.df)-1]

anova.df.trans = as.data.frame(t(anova.df))
colnames(anova.df.trans) = rownames(anova.df)
anova.df.trans$SampleID = rownames(anova.df.trans)

anova.df.trans = merge(anova.df.trans, metadata.df[,c("SampleID","RegionSci")], by = "SampleID")
colnames(anova.df.trans)[length(colnames(anova.df.trans))] = "region"

## generating the final results ----
anova.df.trans = anova.df.trans[,-1]

# Precompute column names and number of columns
col_names <- colnames(anova.df.trans)

# Create the results data frame
results <- data.frame(
  site = col_names[first.site:last.site],
  aov.pval = numeric(length(col_names[first.site:last.site]))
)

######################################################################################################
################### ----------------------------- ANOVA -------------------------- ###################
######################################################################################################

# Define a function to calculate p-value and variance for a given column
calculate_stats <- function(column_name) {
  test <- aov(get(column_name) ~ region, data = anova.df.trans)
  test_result <- summary(test)
  aov_pval <- test_result[[1]][["Pr(>F)"]][1]
  return(aov_pval)
}

# Apply the function to each column and store results
stats <- sapply(col_names[first.site:last.site], calculate_stats)

# Assign results to the data frame
results$aov.pval <- stats

res_dir = "res_dir"
save(results, file = paste0(res_dir,"CpH_Region_ANOVA_",index,".Rdata"))
