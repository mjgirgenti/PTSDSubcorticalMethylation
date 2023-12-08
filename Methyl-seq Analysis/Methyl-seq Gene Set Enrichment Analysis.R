
args = (commandArgs(TRUE))
analysis = args[[1]]
region = args[[2]]
cat("analysis = ", analysis, ", region =", region)

library(dplyr)
library(methylGSA)

load("../protein.coding.genes.RData")
if (analysis == "All") {
  load("../PTSD_DSS_Filtered_NoSmooth_Ancestry_Prop_FinalResults.RData")
} else if (analysis == "Female") {
  load("../PTSD_DSS_Filtered_NoSmooth_Ancestry_Prop_SexSpecific_Female_FinalResults.RData")
} else if (analysis == "Male") {
  load("../PTSD_DSS_Filtered_NoSmooth_Ancestry_Prop_SexSpecific_Male_FinalResults.RData")
}
res_select = res_anno_merged[res_anno_merged$nearestGeneName %in% protein.coding.genes,]

CpG2Gene = res_select[, c("chr_pos", "nearestGeneName")]
colnames(CpG2Gene) = c("CpG", "Gene")
FullAnnot = prepareAnnot(CpG2Gene) 

cpg.pval = res_select[, paste0(region, "_pval")]
names(cpg.pval) = res_select$chr_pos

res = methylglm(cpg.pval = cpg.pval, FullAnnot = FullAnnot, parallel = TRUE)
save(res, file = paste0("../RES_", analysis, "_", region, ".Rdata"))

