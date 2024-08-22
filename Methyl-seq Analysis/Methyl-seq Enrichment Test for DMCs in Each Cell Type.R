
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

CellDMC_res = read.csv(paste0(data_dir, "CellDMC_res.csv"))

######################################################################################################
################### -------------------------- Prepare Data ---------------------- ###################
######################################################################################################

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

totalSig = nrow(CellDMC_res[CellDMC_res$p<=0.0001,])
totalNonSig = nrow(CellDMC_res[CellDMC_res$p>0.0001,])

######################################################################################################
################### ---------------------- Fisher's Exact Test ------------------- ###################
######################################################################################################

## Function to carry out the test ----

calculateChiSquareTest = function(category) {
  category.df = CellDMC_res[CellDMC_res$category == category,]
  
  results.dt = as.data.table(category.df)
  
  categories = c("relativeLocation", "CpG", "Chromatin")
  names(categories) = categories
  
  sig.dt = results.dt[results.dt$p <= 0.0001]
  
  nonsig.dt = results.dt[results.dt$p > 0.0001]
  
  SigCat.list = unlist(lapply(categories, 
                              function(cat) split(sig.dt, by = cat)), recursive = FALSE)
  nonSigCat.list = unlist(lapply(categories, 
                                 function(cat) split(nonsig.dt, by = cat)), recursive = FALSE)
  
  Sig.Num.ByCat = lapply(SigCat.list, function(dt) nrow(dt))
  nonSigNum.ByCat = lapply(nonSigCat.list, function(dt) nrow(dt))
  
  Sig.Num.ByCatOut = lapply(Sig.Num.ByCat, function(n) totalSig - n)
  nonSig.Num.ByCatOut = lapply(nonSigNum.ByCat, function(n) totalNonSig - n)
  
  pval = c()
  for (i in names(Sig.Num.ByCat)) {
    contingency.table = matrix(c(Sig.Num.ByCat[[i]][1],
                                 Sig.Num.ByCatOut[[i]][1],
                                 nonSigNum.ByCat[[i]][1],
                                 nonSig.Num.ByCatOut[[i]][1]), byrow = TRUE, 2, 2)
    test = fisher.test(contingency.table)
    pval = c(pval,test$p.value)
  }
  
  result.df = data.frame(feature = names(Sig.Num.ByCat), pval = pval)
  return(result.df)
}

results = calculateChiSquareTest(args)

res.dir = "res_dir"
save(results, file = paste0(res.dir, "Fishers_Exact_Test_",args,".Rdata"))
