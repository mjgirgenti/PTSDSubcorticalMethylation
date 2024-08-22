
args = (commandArgs(TRUE))

######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(GenomicRanges)
library(data.table)
library(dplyr)

# -------------------------------------------------------------------------------------------------- #
#                                                User Inputs                                         #
# -------------------------------------------------------------------------------------------------- #

data_dir = "data_dir"

captureRegions.df <- read.table(paste0(data_dir, "130912_HG38_CpGiant_4M_EPI.bed"), sep="\t", stringsAsFactors=F, quote="")
captureRegions.gr <- GRanges(seqnames=gsub("chr", "", captureRegions.df$V1), ranges=IRanges(start=sapply(captureRegions.df$V2, function(x) { max(x-50,0) }), end=captureRegions.df$V3+50))
captureRegions.gr <- reduce(captureRegions.gr)

startFile <- as.numeric(args[1])

fileList <- list.files(path=data_dir, pattern="*.tsv")
fileList.path = paste(data_dir, fileList, sep="")
tmp.df <- as.data.frame(fread(fileList.path[startFile]))

######################################################################################################
################### ------------------------ PREPARE DATA ------------------------ ###################
######################################################################################################

tmp.df.CpG = tmp.df[tmp.df$context == "CG", c(1,2,3,7,8)]
tmp.df.CHG = tmp.df[tmp.df$context == "CHG", c(1,2,3,7,8)]
tmp.df.CHH = tmp.df[tmp.df$context == "CHH", c(1,2,3,7,8)]

nm <- gsub("_+", "_", gsub(".tsv", "", fileList[startFile]))

tmp.df.CpG.pos = tmp.df.CpG[tmp.df.CpG$strand == "+",]
tmp.df.CpG.neg = tmp.df.CpG[tmp.df.CpG$strand == "-",]

tmp.df.CpG.neg$pos1 = tmp.df.CpG.neg$pos - 1
tmp.df.CpG.combined = merge(tmp.df.CpG.pos, tmp.df.CpG.neg, by.x = c("chr", "pos"), by.y = c("chr", "pos1"), all = TRUE )

tmp.df.CpG = tmp.df.CpG.combined %>%
  mutate(
    pos = ifelse(is.na(strand.x), pos.y, pos),
    strand = ifelse(is.na(strand.x), "-", ifelse(is.na(strand.y), "+", "*")),
    C_count.x = ifelse(is.na(C_count.x), 0, C_count.x),
    C_count.y = ifelse(is.na(C_count.y), 0, C_count.y),
    CT_count.x = ifelse(is.na(CT_count.x), 0, CT_count.x),
    CT_count.y = ifelse(is.na(CT_count.y), 0, CT_count.y),
    C_count = C_count.x + C_count.y,
    CT_count = CT_count.x + CT_count.y
  ) %>%
  dplyr::select(chr, pos, strand, C_count, CT_count)

# -------------------------------------------------------------------------------------------------- #
#                                             CpG sites                                              #
# -------------------------------------------------------------------------------------------------- #

data.ratio.CpG.gr <- GRanges(seqnames=tmp.df.CpG$chr, 
                             ranges=IRanges(start=tmp.df.CpG$pos, width=rep(2, nrow(tmp.df.CpG))), 
                             ratio=tmp.df.CpG$C_count/tmp.df.CpG$CT_count,
                             strand = Rle(tmp.df.CpG$strand))

data.depth.CpG.gr <- GRanges(seqnames=tmp.df.CpG$chr, 
                             ranges=IRanges(start=tmp.df.CpG$pos, width=rep(2, nrow(tmp.df.CpG))), 
                             ratio=tmp.df.CpG$CT_count,
                             strand = Rle(tmp.df.CpG$strand))

names(mcols(data.ratio.CpG.gr))[1] <- nm
names(mcols(data.depth.CpG.gr))[1] <- nm
rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, data.ratio.CpG.gr))[,2]
data.ratio.CpG.gr <- data.ratio.CpG.gr[rowsToKeep,]
data.depth.CpG.gr <- data.depth.CpG.gr[rowsToKeep,]

# -------------------------------------------------------------------------------------------------- #
#                                             CHG sites                                              #
# -------------------------------------------------------------------------------------------------- #

data.ratio.CHG.gr <- GRanges(seqnames=tmp.df.CHG$chr, 
                             ranges=IRanges(start=tmp.df.CHG$pos, width=rep(2, nrow(tmp.df.CHG))), 
                             ratio=tmp.df.CHG$C_count/tmp.df.CHG$CT_count,
                             strand = Rle(tmp.df.CHG$strand))

data.depth.CHG.gr <- GRanges(seqnames=tmp.df.CHG$chr, 
                             ranges=IRanges(start=tmp.df.CHG$pos, width=rep(2, nrow(tmp.df.CHG))), 
                             ratio=tmp.df.CHG$CT_count,
                             strand = Rle(tmp.df.CHG$strand))

names(mcols(data.ratio.CHG.gr))[1] <- nm
names(mcols(data.depth.CHG.gr))[1] <- nm
rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, data.ratio.CHG.gr))[,2]
data.ratio.CHG.gr <- data.ratio.CHG.gr[rowsToKeep,]
data.depth.CHG.gr <- data.depth.CHG.gr[rowsToKeep,]

# -------------------------------------------------------------------------------------------------- #
#                                             CHH sites                                              #
# -------------------------------------------------------------------------------------------------- #

data.ratio.CHH.gr <- GRanges(seqnames=tmp.df.CHH$chr, 
                             ranges=IRanges(start=tmp.df.CHH$pos, width=rep(2, nrow(tmp.df.CHH))), 
                             ratio=tmp.df.CHH$C_count/tmp.df.CHH$CT_count,
                             strand = Rle(tmp.df.CHH$strand))

data.depth.CHH.gr <- GRanges(seqnames=tmp.df.CHH$chr, 
                             ranges=IRanges(start=tmp.df.CHH$pos, width=rep(2, nrow(tmp.df.CHH))), 
                             ratio=tmp.df.CHH$CT_count,
                             strand = Rle(tmp.df.CHH$strand))

names(mcols(data.ratio.CHH.gr))[1] <- nm
names(mcols(data.depth.CHH.gr))[1] <- nm
rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, data.ratio.CHH.gr))[,2]
data.ratio.CHH.gr <- data.ratio.CHH.gr[rowsToKeep,]
data.depth.CHH.gr <- data.depth.CHH.gr[rowsToKeep,]

# -------------------------------------------------------------------------------------------------- #
#                                       Process all samples                                          #
# -------------------------------------------------------------------------------------------------- #

endFile <- min(startFile + 99, length(fileList))

for( i in (startFile + 1):endFile ) {
  tmp.df <- as.data.frame(fread(fileList.path[i]))
  
  tmp.df.CpG = tmp.df[tmp.df$context == "CG", c(1,2,3,7,8)]
  tmp.df.CHG = tmp.df[tmp.df$context == "CHG", c(1,2,3,7,8)]
  tmp.df.CHH = tmp.df[tmp.df$context == "CHH", c(1,2,3,7,8)]
  
  nm <- gsub("_+", "_", gsub(".tsv", "", fileList[i]))
  
  # CpG sites
  
  ## Combine methylation reads on both strands
  tmp.df.CpG.pos = tmp.df.CpG[tmp.df.CpG$strand == "+",]
  tmp.df.CpG.neg = tmp.df.CpG[tmp.df.CpG$strand == "-",]
  
  tmp.df.CpG.neg$pos1 = tmp.df.CpG.neg$pos - 1
  tmp.df.CpG.combined = merge(tmp.df.CpG.pos, tmp.df.CpG.neg, by.x = c("chr", "pos"), by.y = c("chr", "pos1"), all = TRUE )
  
  tmp.df.CpG = tmp.df.CpG.combined %>%
    mutate(
      pos = ifelse(is.na(strand.x), pos.y, pos),
      strand = ifelse(is.na(strand.x), "-", ifelse(is.na(strand.y), "+", "*")),
      C_count.x = ifelse(is.na(C_count.x), 0, C_count.x),
      C_count.y = ifelse(is.na(C_count.y), 0, C_count.y),
      CT_count.x = ifelse(is.na(CT_count.x), 0, CT_count.x),
      CT_count.y = ifelse(is.na(CT_count.y), 0, CT_count.y),
      C_count = C_count.x + C_count.y,
      CT_count = CT_count.x + CT_count.y
    ) %>%
    dplyr::select(chr, pos, strand, C_count, CT_count)
  
  tmp.ratio.CpG.gr <- GRanges(seqnames=tmp.df.CpG$chr, 
                              ranges=IRanges(start=tmp.df.CpG$pos, width=rep(2, nrow(tmp.df.CpG))), 
                              ratio=tmp.df.CpG$C_count/tmp.df.CpG$CT_count,
                              strand = Rle(tmp.df.CpG$strand))
  
  tmp.depth.CpG.gr <- GRanges(seqnames=tmp.df.CpG$chr, 
                              ranges=IRanges(start=tmp.df.CpG$pos, width=rep(2, nrow(tmp.df.CpG))), 
                              ratio=tmp.df.CpG$CT_count,
                              strand = Rle(tmp.df.CpG$strand))
  
  names(mcols(tmp.ratio.CpG.gr))[1] <- nm
  names(mcols(tmp.depth.CpG.gr))[1] <- nm
  rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, tmp.ratio.CpG.gr))[,2]
  tmp.ratio.CpG.gr <- tmp.ratio.CpG.gr[rowsToKeep,]
  tmp.depth.CpG.gr <- tmp.depth.CpG.gr[rowsToKeep,]
  data.ratio.CpG.gr <- merge(data.ratio.CpG.gr, tmp.ratio.CpG.gr, all=T)
  data.depth.CpG.gr <- merge(data.depth.CpG.gr, tmp.depth.CpG.gr, all=T)
  
  # CHG sites
  tmp.ratio.CHG.gr <- GRanges(seqnames=tmp.df.CHG$chr, 
                              ranges=IRanges(start=tmp.df.CHG$pos, width=rep(2, nrow(tmp.df.CHG))), 
                              ratio=tmp.df.CHG$C_count/tmp.df.CHG$CT_count,
                              strand = Rle(tmp.df.CHG$strand))
  
  tmp.depth.CHG.gr <- GRanges(seqnames=tmp.df.CHG$chr, 
                              ranges=IRanges(start=tmp.df.CHG$pos, width=rep(2, nrow(tmp.df.CHG))), 
                              ratio=tmp.df.CHG$CT_count,
                              strand = Rle(tmp.df.CHG$strand))
  
  names(mcols(tmp.ratio.CHG.gr))[1] <- nm
  names(mcols(tmp.depth.CHG.gr))[1] <- nm
  rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, tmp.ratio.CHG.gr))[,2]
  tmp.ratio.CHG.gr <- tmp.ratio.CHG.gr[rowsToKeep,]
  tmp.depth.CHG.gr <- tmp.depth.CHG.gr[rowsToKeep,]
  data.ratio.CHG.gr <- merge(data.ratio.CHG.gr, tmp.ratio.CHG.gr, all=T)
  data.depth.CHG.gr <- merge(data.depth.CHG.gr, tmp.depth.CHG.gr, all=T)
  
  # CHH sites
  tmp.ratio.CHH.gr <- GRanges(seqnames=tmp.df.CHH$chr, 
                              ranges=IRanges(start=tmp.df.CHH$pos, width=rep(2, nrow(tmp.df.CHH))), 
                              ratio=tmp.df.CHH$C_count/tmp.df.CHH$CT_count,
                              strand = Rle(tmp.df.CHH$strand))
  
  tmp.depth.CHH.gr <- GRanges(seqnames=tmp.df.CHH$chr, 
                              ranges=IRanges(start=tmp.df.CHH$pos, width=rep(2, nrow(tmp.df.CHH))), 
                              ratio=tmp.df.CHH$CT_count,
                              strand = Rle(tmp.df.CHH$strand))
  
  names(mcols(tmp.ratio.CHH.gr))[1] <- nm
  names(mcols(tmp.depth.CHH.gr))[1] <- nm
  rowsToKeep <- as.matrix(findOverlaps(captureRegions.gr, tmp.ratio.CHH.gr))[,2]
  tmp.ratio.CHH.gr <- tmp.ratio.CHH.gr[rowsToKeep,]
  tmp.depth.CHH.gr <- tmp.depth.CHH.gr[rowsToKeep,]
  data.ratio.CHH.gr <- merge(data.ratio.CHH.gr, tmp.ratio.CHH.gr, all=T)
  data.depth.CHH.gr <- merge(data.depth.CHH.gr, tmp.depth.CHH.gr, all=T)
}

# CpG results
res_dir = "res_dir"
save(data.ratio.CpG.gr, file=paste0(res_dir, "data.ratio.padded.CpG.",startFile,".gr.RData"))
save(data.depth.CpG.gr, file=paste0(res_dir, "data.depth.padded.CpG.",startFile,".gr.RData"))

# CHG results
save(data.ratio.CHG.gr, file=paste0(res_dir, "data.ratio.padded.CHG.",startFile,".gr.RData"))
save(data.depth.CHG.gr, file=paste0(res_dir, "data.depth.padded.CHG.",startFile,".gr.RData"))

# CHH results
save(data.ratio.CHH.gr, file=paste0(res_dir, "data.ratio.padded.CHH.",startFile,".gr.RData"))
save(data.depth.CHH.gr, file=paste0(res_dir, "data.depth.padded.CHH.",startFile,".gr.RData"))


