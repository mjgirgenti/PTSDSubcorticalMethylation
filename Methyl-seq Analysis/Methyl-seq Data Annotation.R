
##########################################################################################
## Annotate CpG Features
##########################################################################################

args = (commandArgs(TRUE))
region = args[[1]]
cat("region =", region)

data_folder = "data_folder"
results_folder = "results_folder"

suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(annotatr))
suppressMessages(library(plyranges))
suppressMessages(library(bumphunter))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GenomicRanges))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.refGene))

load(glue({data_folder}, "../AllSamples_", {region}, "_AnnotatedDSS_Results.Rdata"))
data = dml.annotated.gr
data_annotated = data[,-c(1:ncol(mcols(data)))]
seqlevels(data_annotated) = paste0("chr", seqlevels(data_annotated))
data_annotated = data_annotated

cpg_annot = build_annotations(genome = 'hg38', annotations = 'hg38_cpgs')
cpg_annotated = annotate_regions(
  regions = data_annotated,
  annotations = cpg_annot,
  ignore.strand = TRUE,
  quiet = FALSE)

cpg_ref_table = data.frame(old_name = c("hg38_cpg_islands", "hg38_cpg_shores", 
                                        "hg38_cpg_shelves", "hg38_cpg_inter"),
                           new_name = c("Island", "Shore", "Shelf", "Open Sea"))

data_annotated$CpG = cpg_ref_table$new_name[match(cpg_annotated$annot$type, cpg_ref_table$old_name)]


save(data_annotated, file = glue({results_folder}, "AllSamples_DSS_", {region}, "_NewAnnotation_CpG.Rdata"))


##########################################################################################
## Annotate Genic Features
##########################################################################################

args = (commandArgs(TRUE))
region = args[[1]]
cat("region =", region)

data_folder = "data_folder"
results_folder = "results_folder"

suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(annotatr))
suppressMessages(library(plyranges))
suppressMessages(library(bumphunter))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GenomicRanges))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.refGene))

load(glue({data_folder}, "dataDSS/AllSamples_", {region}, "_AnnotatedDSS_Results.Rdata"))
data = dml.annotated.gr
data_annotated = data[,-c(1:ncol(mcols(data)))]
seqlevels(data_annotated) = paste0("chr", seqlevels(data_annotated))
data_annotated = data_annotated

genic_annot = annotateTranscripts(txdb = TxDb.Hsapiens.UCSC.hg38.refGene,
                                  annotationPackage = "org.Hs.eg.db")
genic_annotated = matchGenes(data_annotated, genic_annot)

data_annotated$cpgID = paste(data_annotated@seqnames, data_annotated@ranges@start, sep = "_")
data_annotated$nearestGeneID = genic_annotated$Geneid
data_annotated$nearestGeneName = genic_annotated$name
data_annotated$distanceToNearestGene = genic_annotated$distance
data_annotated$relativeLocation = genic_annotated$description

save(data_annotated, file = glue({results_folder}, "AllSamples_DSS_", {region}, "_NewAnnotation_Genic.Rdata"))

##########################################################################################
## Annotate Chromatin Features
##########################################################################################

args = (commandArgs(TRUE))
region = args[[1]]
cat("region =", region)

data_folder = "data_folder"
results_folder = "results_folder"

suppressMessages(library(glue))
suppressMessages(library(dplyr))
suppressMessages(library(annotatr))
suppressMessages(library(plyranges))
suppressMessages(library(bumphunter))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(GenomicRanges))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.refGene))

load(glue({data_folder}, "dataDSS/AllSamples_", {region}, "_AnnotatedDSS_Results.Rdata"))
data = dml.annotated.gr
data_annotated = data[,-c(1:ncol(mcols(data)))]
seqlevels(data_annotated) = paste0("chr", seqlevels(data_annotated))
data_annotated = data_annotated

chromatin_bedfile = as.data.frame(read.table(glue({data_folder}, "dataReference/hg38lift_genome_100_segments.bed"),
                                             header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = ""))
colnames(chromatin_bedfile) = c("chr", "start", "end", "chromatin_old")
chromatin_annot = makeGRangesFromDataFrame(chromatin_bedfile, keep.extra.columns = TRUE)
chromatin_annot$chromatin_new = gsub("[^a-zA-Z]", "", chromatin_annot$chromatin_old)
chromatin_ref_table = data.frame(old_name = c("GapArtf", "Acet", "EnhWk", "Quies", "TxWk", "TxEx", "EnhA", "ReprPC", 
                                              "HET", "DNase", "PromF", "Tx", "TSS", "znf", "TxEnh", "BivProm"),
                                 new_name = c("Assembly Gaps & Artifacts", "Acetylations", "Weak Enhancer", "Quiescent", 
                                              "Weak Transcription", "Exons & Transcription", "Active Enhancer", 
                                              "Polycomb Repressed", "HET", "DNase", "Flanking Promoter", "Strong Transcription", 
                                              "TSS", "ZNF Genes", "Transcribed & Enhancers", "Bivalent Promoter"))
chromatin_annot$chromatin = chromatin_ref_table$new_name[match(chromatin_annot$chromatin_new, chromatin_ref_table$old_name)]

chromatin_annotated = annotate_regions(
  regions = data_annotated,
  annotations = chromatin_annot,
  ignore.strand = TRUE,
  quiet = FALSE)

chromatin_annotated$Chromatin = chromatin_annotated$annot$chromatin
chromatin_annotated$CpG_ID = paste(chromatin_annotated@seqnames, 
                                   chromatin_annotated@ranges@start, sep = "_")
chromatin_annotated = chromatin_annotated[!duplicated(chromatin_annotated$CpG_ID)]

overlaps = findOverlaps(chromatin_annotated, data_annotated)
mcols(data_annotated)$Chromatin <- NA
mcols(data_annotated)[subjectHits(overlaps), "Chromatin"] = mcols(chromatin_annotated)[queryHits(overlaps), "Chromatin"]

save(data_annotated, file = glue({results_folder}, "AllSamples_DSS_", {region}, "_NewAnnotation_Chromatin.Rdata"))



