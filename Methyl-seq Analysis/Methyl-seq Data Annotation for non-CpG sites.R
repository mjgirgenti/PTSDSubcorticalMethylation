
######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

data_dir = "data_dir"

load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHG.gr.RData"))
load(paste0(data_dir, "data.ratio.normalized.snpFiltered.invariantFiltered.overlapped.CHH.gr.RData"))

suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(plyranges))
suppressMessages(library(bumphunter))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GenomicRanges))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.refGene))

######################################################################################################
################### ------------------------- CpH Annotation --------------------- ###################
######################################################################################################

CHG.data = as.data.frame(data.ratio.CHG.gr)
CHH.data = as.data.frame(data.ratio.CHH.gr)

CHG.data$seqnames = paste0("chr",CHG.data$seqnames)
CHH.data$seqnames = paste0("chr",CHH.data$seqnames)

# Reference: UCSC v38
reference = annotateTranscripts(txdb = TxDb.Hsapiens.UCSC.hg38.refGene,
                                annotationPackage = "org.Hs.eg.db")
# Match CHG sites to genes
anno_res_CHG = matchGenes(CHG.data, reference)
CHG_results = cbind(CHG.data[,1:2],anno_res_CHG)

# Match CHH sites to genes
anno_res_CHH = matchGenes(CHH.data, reference)
CHH_results = cbind(CHH.data[,1:2],anno_res_CHH)

CpH_results = rbind(CHG_results,CHH_results)

res_dir = "res_dir"

save(CHG_results, file=paste0(res_dir, "CHG_gene_annotation.RData"))
save(CHH_results, file=paste0(res_dir, "CHH_gene_annotation.RData"))
save(CpH_results, file=paste0(res_dir, "CpH_gene_annotation.RData"))

