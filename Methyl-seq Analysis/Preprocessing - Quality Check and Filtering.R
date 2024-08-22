
args = (commandArgs(TRUE))
type = args[1]

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(GenomicRanges)
library(data.table)

data_dir = "data_dir"
res_dir = "res_dir"

######################################################################################################
################### ------------------- Quality Check and Filtering -------------- ###################
######################################################################################################

if (type == "CpG") {
  # CpG ----
  load(paste0(data_dir, "CpG results/data.ratio.padded.CpG.1.gr.RData"))
  load(paste0(data_dir, "CpG results/data.depth.padded.CpG.1.gr.RData"))
  tmp.ratio.CpG.gr <- data.ratio.CpG.gr 
  tmp.depth.CpG.gr <- data.depth.CpG.gr
  
  for( i in seq(101,901,by=100) ) {
    print(i)
    
    load(paste0(data_dir, "CpG results/data.ratio.padded.CpG.",i,".gr.RData"))
    load(paste0(data_dir, "CpG results/data.depth.padded.CpG.",i,".gr.RData"))
    tmp.ratio.CpG.gr <- merge(tmp.ratio.CpG.gr, data.ratio.CpG.gr, all=T)
    tmp.depth.CpG.gr <- merge(tmp.depth.CpG.gr, data.depth.CpG.gr, all=T)
  }
  
  data.ratio.CpG.gr <- tmp.ratio.CpG.gr
  data.depth.CpG.gr <- tmp.depth.CpG.gr
  remove(i)
  remove(tmp.ratio.CpG.gr)
  remove(tmp.depth.CpG.gr)
  
  names(mcols(data.ratio.CpG.gr)) <- gsub("^X", "L", names(mcols(data.ratio.CpG.gr)))
  names(mcols(data.depth.CpG.gr)) <- gsub("^X", "L", names(mcols(data.depth.CpG.gr)))
  
  ### Remove the 83 year old man
  data.ratio.CpG.gr <- data.ratio.CpG.gr[,-which(names(mcols(data.ratio.CpG.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  data.depth.CpG.gr <- data.depth.CpG.gr[,-which(names(mcols(data.depth.CpG.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  
  ### Fix a couple of the sample names
  colnames(mcols(data.ratio.CpG.gr))[which(colnames(mcols(data.ratio.CpG.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.ratio.CpG.gr))[which(colnames(mcols(data.ratio.CpG.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.ratio.CpG.gr))[which(colnames(mcols(data.ratio.CpG.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.ratio.CpG.gr))[which(colnames(mcols(data.ratio.CpG.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.ratio.CpG.gr))[which(colnames(mcols(data.ratio.CpG.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  colnames(mcols(data.depth.CpG.gr))[which(colnames(mcols(data.depth.CpG.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.depth.CpG.gr))[which(colnames(mcols(data.depth.CpG.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.depth.CpG.gr))[which(colnames(mcols(data.depth.CpG.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.depth.CpG.gr))[which(colnames(mcols(data.depth.CpG.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.depth.CpG.gr))[which(colnames(mcols(data.depth.CpG.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  ### Identify bad samples: 
  ### 2 indicate columns in apply function
  depths <- t(apply(mcols(data.depth.CpG.gr), 2, function(x) { c(mean(x, na.rm=T), sd(x, na.rm=T)) }))
  
  ### 1 indicate rows in apply function
  ### calculate the proportion that of the samples not have methylation depth less than 10 or greater than 125 or is missing
  propPresent <- apply(mcols(data.depth.CpG.gr), 1, function(x) { 
    1 - (length(which(is.na(x) | x < 10 | x > 125 )) / length(x))
  })
  
  ### remove the sites that have proportion less than 90%
  data.ratio.CpG.gr <- data.ratio.CpG.gr[-which(propPresent < 0.9),]
  data.depth.CpG.gr <- data.depth.CpG.gr[-which(propPresent < 0.9),]
  
  ### for each sample, calculate the number of site that have methylation depth less than 10 or greater than 125 or is missing
  numBadSites <- apply(mcols(data.depth.CpG.gr), 2, function(x) { 
    length(which(is.na(x) | x < 10 | x > 125))
  })
  
  print(max(numBadSites))
  print(min(numBadSites))
  print(-which(numBadSites > 150000))
  print(length(numBadSites))
  print(lenth(names(mcols(data.depth.CpG.gr))))
  
  ### remove the samples that have bad sites greater than 150000
  data.ratio.CpG.gr <- data.ratio.CpG.gr[,-which(numBadSites > 150000)]
  data.depth.CpG.gr <- data.depth.CpG.gr[,-which(numBadSites > 150000)]
  
  # data.ratio.CpG.gr <- data.ratio.CpG.gr[,numBadSites <= 150000]
  # data.depth.CpG.gr <- data.depth.CpG.gr[,numBadSites <= 150000]
  
  print(names(mcols(data.depth.CpG.gr)))
  
  ### change the sites that have methylation site less than 10 or greater than 125 to NA
  data.mat <- as.matrix(mcols(data.depth.CpG.gr))
  toReplace <- which(data.mat < 10 | data.mat > 125)
  
  data.mat[toReplace] <- NA
  
  for(i in 1:ncol(data.mat) ) {
    mcols(data.depth.CpG.gr)[,i] <- data.mat[,i]
  }
  
  data.mat <- as.matrix(mcols(data.ratio.CpG.gr))
  data.mat[toReplace] <- NA
  for(i in 1:ncol(data.mat) ) {
    mcols(data.ratio.CpG.gr)[,i] <- data.mat[,i]
  }
  
  save(data.depth.CpG.gr, file=paste0(res_dir, "CpG results/data.depth.filtered.CpG.gr.RData"))
  save(data.ratio.CpG.gr, file=paste0(res_dir, "CpG results/data.ratio.filtered.CpG.gr.RData"))
  
  ### scaling the median coverage of the sequencing depth
  medianCoverage <- apply(mcols(data.depth.CpG.gr), 2, function(x) { median(x, na.rm=T) })
  
  scaleFactors <- max(medianCoverage) / medianCoverage 
  
  for(i in 1:ncol(mcols(data.depth.CpG.gr)) ) {
    mcols(data.depth.CpG.gr)[,i] <- round(mcols(data.depth.CpG.gr)[,i] * scaleFactors[i])
    mcols(data.ratio.CpG.gr)[,i] <- round(mcols(data.depth.CpG.gr)[,i] * mcols(data.ratio.CpG.gr)[,i]) / mcols(data.depth.CpG.gr)[,i]
  }
  
  save(data.depth.CpG.gr, file=paste0(res_dir, "CpG results/data.depth.normalized.CpG.gr.RData"))
  save(data.ratio.CpG.gr, file=paste0(res_dir, "CpG results/data.ratio.normalized.CpG.gr.RData"))
  
  # remove the common SNPs identified by dbSNPs
  dbSNP.df <- read.table("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/Methylation_Raw_Data/Preprocessing/00-common_all.vcf", sep="\t", header=F, stringsAsFactors=F, colClasses = c("character", "integer", "character", rep("NULL", 5)))
  
  dbSNP.gr <- GRanges(seqnames=dbSNP.df$V1, ranges=IRanges(start=dbSNP.df$V2, width=rep(1, length(dbSNP.df$V2))), rsID=dbSNP.df$V3)
  remove(dbSNP.df)
  
  overlapList <- as.matrix(findOverlaps(data.depth.CpG.gr, dbSNP.gr))
  
  sitesToRemove <- unique(overlapList[,1])
  data.depth.CpG.gr <- data.depth.CpG.gr[-sitesToRemove,]
  data.ratio.CpG.gr <- data.ratio.CpG.gr[-sitesToRemove,]
  # 105477 / 2090236
  
  save(data.depth.CpG.gr, file=paste0(res_dir, "CpG results/data.depth.normalized.snpFiltered.CpG.gr.RData"))
  save(data.ratio.CpG.gr, file=paste0(res_dir, "CpG results/data.ratio.normalized.snpFiltered.CpG.gr.RData"))
  
  # remove the invariant methylation sites
  ann.tr.df <- as.matrix(mcols(data.ratio.CpG.gr))
  decileInfo.mat <- t(apply(ann.tr.df, 1, function(x) { quantile(x, probs=c(.1,.9), na.rm=TRUE) }))
  
  meth.stats.dt <- data.table(
    Mean=unlist(apply(ann.tr.df, 1, function(x) mean(x, na.rm=TRUE))),
    MAD=unlist(apply(ann.tr.df, 1, function(x) mad(x, na.rm=TRUE))),
    decile10=decileInfo.mat[,1],
    decile90=decileInfo.mat[,2])
  meth.stats.dt$decile90_10_dif <- meth.stats.dt$decile90-meth.stats.dt$decile10
  
  data.depth.CpG.gr <- data.depth.CpG.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  data.ratio.CpG.gr <- data.ratio.CpG.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  
  save(data.depth.CpG.gr, file=paste0(res_dir, "CpG results/data.depth.normalized.snpFiltered.invariantFiltered.CpG.gr.RData"))
  save(data.ratio.CpG.gr, file=paste0(res_dir, "CpG results/data.ratio.normalized.snpFiltered.invariantFiltered.CpG.gr.RData"))
  
} else if (type == "CHG") {
  # CHG ----
  load(paste0(data_dir, "CHG results/data.ratio.padded.CHG.1.gr.RData"))
  load(paste0(data_dir, "CHG results/data.depth.padded.CHG.1.gr.RData"))
  tmp.ratio.CHG.gr <- data.ratio.CHG.gr 
  tmp.depth.CHG.gr <- data.depth.CHG.gr
  
  for( i in seq(101,901,by=100) ) {
    print(i)
    
    load(paste0(data_dir, "CHG results/data.ratio.padded.CHG.",i,".gr.RData"))
    load(paste0(data_dir, "CHG results/data.depth.padded.CHG.",i,".gr.RData"))
    
    tmp.ratio.CHG.gr <- merge(tmp.ratio.CHG.gr, data.ratio.CHG.gr, all=T)
    tmp.depth.CHG.gr <- merge(tmp.depth.CHG.gr, data.depth.CHG.gr, all=T)
  }
  
  data.ratio.CHG.gr <- tmp.ratio.CHG.gr
  data.depth.CHG.gr <- tmp.depth.CHG.gr
  remove(i)
  remove(tmp.ratio.CHG.gr)
  remove(tmp.depth.CHG.gr)
  
  names(mcols(data.ratio.CHG.gr)) <- gsub("^X", "L", names(mcols(data.ratio.CHG.gr)))
  names(mcols(data.depth.CHG.gr)) <- gsub("^X", "L", names(mcols(data.depth.CHG.gr)))
  
  ### Remove the 83 year old man
  data.ratio.CHG.gr <- data.ratio.CHG.gr[,-which(names(mcols(data.ratio.CHG.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  data.depth.CHG.gr <- data.depth.CHG.gr[,-which(names(mcols(data.depth.CHG.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  
  ### Fix a couple of the sample names
  colnames(mcols(data.ratio.CHG.gr))[which(colnames(mcols(data.ratio.CHG.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.ratio.CHG.gr))[which(colnames(mcols(data.ratio.CHG.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.ratio.CHG.gr))[which(colnames(mcols(data.ratio.CHG.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.ratio.CHG.gr))[which(colnames(mcols(data.ratio.CHG.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.ratio.CHG.gr))[which(colnames(mcols(data.ratio.CHG.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  colnames(mcols(data.depth.CHG.gr))[which(colnames(mcols(data.depth.CHG.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.depth.CHG.gr))[which(colnames(mcols(data.depth.CHG.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.depth.CHG.gr))[which(colnames(mcols(data.depth.CHG.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.depth.CHG.gr))[which(colnames(mcols(data.depth.CHG.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.depth.CHG.gr))[which(colnames(mcols(data.depth.CHG.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  ### Identify bad samples:
  depths <- t(apply(mcols(data.depth.CHG.gr), 2, function(x) { c(mean(x, na.rm=T), sd(x, na.rm=T)) }))
  
  
  propPresent <- apply(mcols(data.depth.CHG.gr), 1, function(x) { 
    1 - (length(which(is.na(x) | x < 10 | x > 125 )) / length(x))
  })
  
  
  data.ratio.CHG.gr <- data.ratio.CHG.gr[-which(propPresent < 0.9),]
  data.depth.CHG.gr <- data.depth.CHG.gr[-which(propPresent < 0.9),]
  
  
  numBadSites <- apply(mcols(data.depth.CHG.gr), 2, function(x) { 
    length(which(is.na(x) | x < 10 | x > 125))
  })
  
  
  data.ratio.CHG.gr <- data.ratio.CHG.gr[,-which(numBadSites > 150000)]
  data.depth.CHG.gr <- data.depth.CHG.gr[,-which(numBadSites > 150000)]
  
  data.mat <- as.matrix(mcols(data.depth.CHG.gr))
  toReplace <- which(data.mat < 10 | data.mat > 125)
  
  data.mat[toReplace] <- NA
  
  for(i in 1:ncol(data.mat) ) {
    mcols(data.depth.CHG.gr)[,i] <- data.mat[,i]
  }
  
  data.mat <- as.matrix(mcols(data.ratio.CHG.gr))
  data.mat[toReplace] <- NA
  for(i in 1:ncol(data.mat) ) {
    mcols(data.ratio.CHG.gr)[,i] <- data.mat[,i]
  }
  
  save(data.depth.CHG.gr, file=paste0(res_dir, "CHG results/data.depth.filtered.CHG.gr.RData"))
  save(data.ratio.CHG.gr, file=paste0(res_dir, "CHG results/data.ratio.filtered.CHG.gr.RData"))
  
  
  medianCoverage <- apply(mcols(data.depth.CHG.gr), 2, function(x) { median(x, na.rm=T) })
  
  scaleFactors <- max(medianCoverage) / medianCoverage 
  
  for(i in 1:ncol(mcols(data.depth.CHG.gr)) ) {
    mcols(data.depth.CHG.gr)[,i] <- round(mcols(data.depth.CHG.gr)[,i] * scaleFactors[i])
    mcols(data.ratio.CHG.gr)[,i] <- round(mcols(data.depth.CHG.gr)[,i] * mcols(data.ratio.CHG.gr)[,i]) / mcols(data.depth.CHG.gr)[,i]
  }
  
  save(data.depth.CHG.gr, file=paste0(res_dir, "CHG results/data.depth.normalized.CHG.gr.RData"))
  save(data.ratio.CHG.gr, file=paste0(res_dir, "CHG results/data.ratio.normalized.CHG.gr.RData"))
  
  
  dbSNP.df <- read.table("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/Methylation_Raw_Data/Preprocessing/00-common_all.vcf", sep="\t", header=F, stringsAsFactors=F, colClasses = c("character", "integer", "character", rep("NULL", 5)))
  
  dbSNP.gr <- GRanges(seqnames=dbSNP.df$V1, ranges=IRanges(start=dbSNP.df$V2, width=rep(1, length(dbSNP.df$V2))), rsID=dbSNP.df$V3)
  remove(dbSNP.df)
  
  overlapList <- as.matrix(findOverlaps(data.depth.CHG.gr, dbSNP.gr))
  
  sitesToRemove <- unique(overlapList[,1])
  data.depth.CHG.gr <- data.depth.CHG.gr[-sitesToRemove,]
  data.ratio.CHG.gr <- data.ratio.CHG.gr[-sitesToRemove,]
  # 105477 / 2090236
  
  save(data.depth.CHG.gr, file=paste0(res_dir, "CHG results/data.depth.normalized.snpFiltered.CHG.gr.RData"))
  save(data.ratio.CHG.gr, file=paste0(res_dir, "CHG results/data.ratio.normalized.snpFiltered.CHG.gr.RData"))
  
  ann.tr.df <- as.matrix(mcols(data.ratio.CHG.gr))
  decileInfo.mat <- t(apply(ann.tr.df, 1, function(x) { quantile(x, probs=c(.1,.9), na.rm=TRUE) }))
  
  meth.stats.dt <- data.table(
    Mean=unlist(apply(ann.tr.df, 1, function(x) mean(x, na.rm=TRUE))),
    MAD=unlist(apply(ann.tr.df, 1, function(x) mad(x, na.rm=TRUE))),
    decile10=decileInfo.mat[,1],
    decile90=decileInfo.mat[,2])
  meth.stats.dt$decile90_10_dif <- meth.stats.dt$decile90-meth.stats.dt$decile10
  
  data.depth.CHG.gr <- data.depth.CHG.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  data.ratio.CHG.gr <- data.ratio.CHG.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  
  save(data.depth.CHG.gr, file=paste0(res_dir, "CHG results/data.depth.normalized.snpFiltered.invariantFiltered.CHG.gr.RData"))
  save(data.ratio.CHG.gr, file=paste0(res_dir, "CHG results/data.ratio.normalized.snpFiltered.invariantFiltered.CHG.gr.RData"))
  
} else if (type == "CHH") {
  # CHH ----
  load(paste0(data_dir, "CHH results/data.ratio.padded.CHH.1.gr.RData"))
  load(paste0(data_dir, "CHH results/data.depth.padded.CHH.1.gr.RData"))
  tmp.ratio.CHH.gr <- data.ratio.CHH.gr 
  tmp.depth.CHH.gr <- data.depth.CHH.gr
  
  for( i in seq(101,901,by=100) ) {
    print(i)
    
    load(paste0(data_dir, "CHH results/data.ratio.padded.CHH.",i,".gr.RData"))
    load(paste0(data_dir, "CHH results/data.depth.padded.CHH.",i,".gr.RData"))
    
    tmp.ratio.CHH.gr <- merge(tmp.ratio.CHH.gr, data.ratio.CHH.gr, all=T)
    tmp.depth.CHH.gr <- merge(tmp.depth.CHH.gr, data.depth.CHH.gr, all=T)
  }
  
  data.ratio.CHH.gr <- tmp.ratio.CHH.gr
  data.depth.CHH.gr <- tmp.depth.CHH.gr
  remove(i)
  remove(tmp.ratio.CHH.gr)
  remove(tmp.depth.CHH.gr)
  
  names(mcols(data.ratio.CHH.gr)) <- gsub("^X", "L", names(mcols(data.ratio.CHH.gr)))
  names(mcols(data.depth.CHH.gr)) <- gsub("^X", "L", names(mcols(data.depth.CHH.gr)))
  
  ### Remove the 83 year old man
  data.ratio.CHH.gr <- data.ratio.CHH.gr[,-which(names(mcols(data.ratio.CHH.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  data.depth.CHH.gr <- data.depth.CHH.gr[,-which(names(mcols(data.depth.CHH.gr)) %in% c("L5149_amygA", "L5149_amygB", "L5149_amygC", "L5149_hipA", "L5149_hipC"))]
  
  ### Fix a couple of the sample names
  colnames(mcols(data.ratio.CHH.gr))[which(colnames(mcols(data.ratio.CHH.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.ratio.CHH.gr))[which(colnames(mcols(data.ratio.CHH.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.ratio.CHH.gr))[which(colnames(mcols(data.ratio.CHH.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.ratio.CHH.gr))[which(colnames(mcols(data.ratio.CHH.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.ratio.CHH.gr))[which(colnames(mcols(data.ratio.CHH.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  colnames(mcols(data.depth.CHH.gr))[which(colnames(mcols(data.depth.CHH.gr)) == "L5626amygA")] <- "L5626_amygA"
  colnames(mcols(data.depth.CHH.gr))[which(colnames(mcols(data.depth.CHH.gr)) == "L5626amygB")] <- "L5626_amygB"
  colnames(mcols(data.depth.CHH.gr))[which(colnames(mcols(data.depth.CHH.gr)) == "L5626amygC")] <- "L5626_amygC"
  colnames(mcols(data.depth.CHH.gr))[which(colnames(mcols(data.depth.CHH.gr)) == "L5236ant_hip_A")] <- "L5236_hipA"
  colnames(mcols(data.depth.CHH.gr))[which(colnames(mcols(data.depth.CHH.gr)) == "L5236ant_hip_B")] <- "L5236_hipB"
  
  ### Identify bad samples:
  depths <- t(apply(mcols(data.depth.CHH.gr), 2, function(x) { c(mean(x, na.rm=T), sd(x, na.rm=T)) }))
  
  
  propPresent <- apply(mcols(data.depth.CHH.gr), 1, function(x) { 
    1 - (length(which(is.na(x) | x < 10 | x > 125 )) / length(x))
  })
  
  
  data.ratio.CHH.gr <- data.ratio.CHH.gr[-which(propPresent < 0.9),]
  data.depth.CHH.gr <- data.depth.CHH.gr[-which(propPresent < 0.9),]
  
  
  numBadSites <- apply(mcols(data.depth.CHH.gr), 2, function(x) { 
    length(which(is.na(x) | x < 10 | x > 125))
  })
  
  
  data.ratio.CHH.gr <- data.ratio.CHH.gr[,-which(numBadSites > 150000)]
  data.depth.CHH.gr <- data.depth.CHH.gr[,-which(numBadSites > 150000)]
  
  data.mat <- as.matrix(mcols(data.depth.CHH.gr))
  toReplace <- which(data.mat < 10 | data.mat > 125)
  
  data.mat[toReplace] <- NA
  
  for(i in 1:ncol(data.mat) ) {
    mcols(data.depth.CHH.gr)[,i] <- data.mat[,i]
  }
  
  data.mat <- as.matrix(mcols(data.ratio.CHH.gr))
  data.mat[toReplace] <- NA
  for(i in 1:ncol(data.mat) ) {
    mcols(data.ratio.CHH.gr)[,i] <- data.mat[,i]
  }
  
  save(data.depth.CHH.gr, file=paste0(res_dir, "CHH results/data.depth.filtered.CHH.gr.RData"))
  save(data.ratio.CHH.gr, file=paste0(res_dir, "CHH results/data.ratio.filtered.CHH.gr.RData"))
  
  
  medianCoverage <- apply(mcols(data.depth.CHH.gr), 2, function(x) { median(x, na.rm=T) })
  
  scaleFactors <- max(medianCoverage) / medianCoverage 
  
  for(i in 1:ncol(mcols(data.depth.CHH.gr)) ) {
    mcols(data.depth.CHH.gr)[,i] <- round(mcols(data.depth.CHH.gr)[,i] * scaleFactors[i])
    mcols(data.ratio.CHH.gr)[,i] <- round(mcols(data.depth.CHH.gr)[,i] * mcols(data.ratio.CHH.gr)[,i]) / mcols(data.depth.CHH.gr)[,i]
  }
  
  save(data.depth.CHH.gr, file=paste0(res_dir, "CHH results/data.depth.normalized.CHH.gr.RData"))
  save(data.ratio.CHH.gr, file=paste0(res_dir, "CHH results/data.ratio.normalized.CHH.gr.RData"))
  
  
  dbSNP.df <- read.table("/gpfs/gibbs/pi/girgenti/yjliu/Methylation_Analysis/Methylation_Raw_Data/Preprocessing/00-common_all.vcf", sep="\t", header=F, stringsAsFactors=F, colClasses = c("character", "integer", "character", rep("NULL", 5)))
  
  dbSNP.gr <- GRanges(seqnames=dbSNP.df$V1, ranges=IRanges(start=dbSNP.df$V2, width=rep(1, length(dbSNP.df$V2))), rsID=dbSNP.df$V3)
  remove(dbSNP.df)
  
  overlapList <- as.matrix(findOverlaps(data.depth.CHH.gr, dbSNP.gr))
  
  sitesToRemove <- unique(overlapList[,1])
  data.depth.CHH.gr <- data.depth.CHH.gr[-sitesToRemove,]
  data.ratio.CHH.gr <- data.ratio.CHH.gr[-sitesToRemove,]
  # 105477 / 2090236
  
  save(data.depth.CHH.gr, file=paste0(res_dir, "CHH results/data.depth.normalized.snpFiltered.CHH.gr.RData"))
  save(data.ratio.CHH.gr, file=paste0(res_dir, "CHH results/data.ratio.normalized.snpFiltered.CHH.gr.RData"))
  
  ann.tr.df <- as.matrix(mcols(data.ratio.CHH.gr))
  decileInfo.mat <- t(apply(ann.tr.df, 1, function(x) { quantile(x, probs=c(.1,.9), na.rm=TRUE) }))
  
  meth.stats.dt <- data.table(
    Mean=unlist(apply(ann.tr.df, 1, function(x) mean(x, na.rm=TRUE))),
    MAD=unlist(apply(ann.tr.df, 1, function(x) mad(x, na.rm=TRUE))),
    decile10=decileInfo.mat[,1],
    decile90=decileInfo.mat[,2])
  meth.stats.dt$decile90_10_dif <- meth.stats.dt$decile90-meth.stats.dt$decile10
  
  data.depth.CHH.gr <- data.depth.CHH.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  data.ratio.CHH.gr <- data.ratio.CHH.gr[-which(meth.stats.dt$decile90_10_dif < 0.05),]
  
  save(data.depth.CHH.gr, file=paste0(res_dir, "CHH results/data.depth.normalized.snpFiltered.invariantFiltered.CHH.gr.RData"))
  save(data.ratio.CHH.gr, file=paste0(res_dir, "CHH results/data.ratio.normalized.snpFiltered.invariantFiltered.CHH.gr.RData"))
  
}
