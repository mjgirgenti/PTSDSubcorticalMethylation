
############################################################
## Obtain Reference Matrix with Single Cell Data
############################################################

library(plyr)
library(MASS)
library(org.Hs.eg.db)
library(annotate)
library(EpiSCORE)
library(ggplot2)
library(ggsci)
library(patchwork)
library(data.table)
library(EpiDISH)
library(parallel)
library(pbapply)
library(matrixStats)
library(enrichR)
library(pbmcapply)
library(rentrez)
library(stringr)
library(DescTools)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(latex2exp)
library(Seurat)
library(Rtsne)
library(ggbeeswarm)
library(umap)
source("../CellTypeProportion_HelperFunctions_1.R")

## Load AD Nancy Data
sc_obj <- readRDS("../seurat_obj_cell_type_labelled.rds")
sc_counts <- GetAssayData(object = sc_obj, slot = "counts", assay = "RNA")
table(sc_obj$cell_type)
# Remove data with cell type = NA
sc_cts <- sc_obj$cell_type
rm_ct_index <- which(is.na(sc_cts))
sc_counts <- sc_counts[,-rm_ct_index]
sc_cts <- sc_cts[-rm_ct_index]
sc_cts <- as.vector(sc_cts)
# Merge Ex and In to Neu+
neuron_cells_index <- which(sc_cts == "Ex" | sc_cts == "In")
sc_cts[neuron_cells_index] <- "Neu+"
# Merge Ast, End, Mic and Oli to Neu-
non_neuron_cells_index <- which(sc_cts == "Ast" | sc_cts == "End" | 
                                  sc_cts == "Mic" | sc_cts == "Oli")
sc_cts[non_neuron_cells_index] <- "Neu-"

ct.idx <- as.integer(factor(sc_cts))
names_ct.v <- levels(factor(sc_cts))

gc(verbose = F)
# Step 1, normalizing the read count matrix 
message("normalizing the count matrix")
maximum_read_count <- max(colSums(sc_counts))
sc_counts <- pbapply(sc_counts, 2, function(cell){
  log(cell*(maximum_read_count/sum(cell)) + 1, base = 2)
})

gc(verbose = F)
# Step 2, constructing the mRNA reference matrix
message("constructing the mRNA reference matrix")
MSS = c(1, 1)
mRNA_ref.m <- ConstExpRef(exp.m = sc_counts, celltype.idx = ct.idx, 
                          namesCellT.v = names_ct.v, markspecTH.v = MSS)

gc(verbose = F)
# Step 3, imputing the DNAm reference matrix
message("imputing the DNAm matrix")
RNA_ref.m <- mRNA_ref.m$ref
DNAm_ref.m <- ImputeDNAmRef(refexp.m = RNA_ref.m$med, db = "RMAP", geneID = "SYMBOL")
NeuProp_ref = list(mRNA_ref = RNA_ref.m, DNAm_ref = DNAm_ref.m)
save(NeuProp_ref, file = "NeuProp_ref.RData")


############################################################
## Prepare the Methylation Data
############################################################

library(glue)
library(dplyr)
library(limma)
library(jamba)
library(readxl)
library(plyranges)
library(bumphunter)
library(org.Hs.eg.db)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.refGene)

main_dir = "DIRECTORY"
load(paste0(main_dir, "data/data.ratio.normalized.snpFiltered.invariantFiltered.gr.RData"))
beta.all = data.frame(data.ratio.gr@elementMetadata)
rownames(beta.all) = paste0("chr", data.ratio.gr@seqnames, "_", data.ratio.gr@ranges@start)
colnames(beta.all) = colnames(data.ratio.gr@elementMetadata)
colnames(beta.all)[77] = "L5236_hipA"
colnames(beta.all)[188] = "L5626_amygA"

load(paste0("/gpfs/gibbs/pi/girgenti/hl732/methylation/dataAnnotated/PTSD_Combined/PTSD_Combined_BLA_Annotated.Robj"))
cpg.meta.all = data.frame(data_annotated@elementMetadata)[, c("cpgID", "nearestGeneID", "nearestGeneName", 
                                                              "distanceToNearestGene", "relativeLocation")]
cpg.meta.all = cbind(data_annotated@seqnames, cpg.meta.all)
rownames(cpg.meta.all) = cpg.meta.all$cpgID
colnames(cpg.meta.all)[1] = "chr"

meta.clinical.all = read_excel(paste0(main_dir, "data/MetadataFinal041720.xlsx"))
meta.clinical.all %>%
  dplyr::select(SampleID, Age, Sex, Race, Source, PMI, DxGroup1) -> meta.clinical.all
meta.clinical.all = data.frame(meta.clinical.all)
rownames(meta.clinical.all) = meta.clinical.all$SampleID

col.common = sort(intersect(colnames(beta.all), meta.clinical.all$SampleID))
beta.all = beta.all[, col.common]
beta.all = beta.all[rowSums(is.na(beta.all)) == 0, ]
row.common = sort(intersect(rownames(beta.all), cpg.meta.all$cpgID))
beta = beta.all[row.common, ]

cpg.meta = cpg.meta.all[row.common, ]
meta.clinical = meta.clinical.all[col.common, ]
colnames(meta.clinical) = c("sampleID", "age", "sex", "race", "brainBank", "PMI", "diagnosis")

data = list(beta = beta, cpg.meta = cpg.meta, meta.clinical = meta.clinical)

save(data, file = "../data.PTSD.RData")

############################################################
## Get Cell Type Proportion using Deconvolution
############################################################

# Load Packages
library(glue)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(plyranges)
library(GenomicRanges)

##############################################################################
## Functions
##############################################################################

source("../CellTypeProportion_HelperFunctions_2.R")

load("../NeuProp_ref.RData")
datMethylSig = NeuProp_ref$DNAm_ref
datMethylSig %>% dim()
head(datMethylSig)
datMethylSig = datMethylSig[complete.cases(datMethylSig), 1:2]
head(datMethylSig)

load("../data.PTSD.RData")
beta = data$beta
beta[1:5, 1:5]
meta = data$cpg.meta
head(meta)
idx = which(meta$relativeLocation == "promoter")

X = beta[idx,]
X$annot = meta$nearestGeneID[idx]
XX = aggregate(.~annot, X, mean)
rownames(XX) = XX$annot
datMethylBulk_AvgByGene = XX[,-1]

selectGeneID = intersect(rownames(datMethylSig), rownames(datMethylBulk_AvgByGene))
bulk = datMethylBulk_AvgByGene[match(selectGeneID, rownames(datMethylBulk_AvgByGene)),]
signature = datMethylSig[match(selectGeneID, rownames(datMethylSig)),]
sig = signature

bulk %>% dim()
sig %>% dim()

# Estimate CTS Proportion 
CTS_Proportion = constraint_ols(sig = sig, bulk = bulk)
dim(CTS_Proportion)
head(CTS_Proportion)
colnames(CTS_Proportion) = colnames(sig)
rownames(CTS_Proportion) = colnames(bulk)
CTS_Proportion = as.data.frame(CTS_Proportion)
save(CTS_Proportion, file = "../MethylationCellTypeProportions.RData")


