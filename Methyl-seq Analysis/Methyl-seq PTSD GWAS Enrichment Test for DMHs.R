
######################################################################################################
################### ------------------ DEPENDENCIES AND USER INPUTS -------------- ###################
######################################################################################################

library(glue)
library(xlsx)
library(dplyr)
library(readr)
library(readxl)
library(plyranges)
library(GenomicRanges)
library(Repitools)

data_dir = "data_dir"

## import CHG data ----
region = c("BLA", "CNA", "MNA", "CA", "DG", "Sub")

chr.num = c(seq(1, 22, 1))

for (i in 1:length(region)) {
  
  PTSD.CHG.dataframe.list = list()
  
  for (j in 1:length(chr.num)) {
    load(paste0(data_dir, 
                "DSS_Filtered_NoSmooth_Ancestry_Prop_Results_",
                region[i],"_PTSD_", "chr", chr.num[j], "_CHG.Rdata"))
    
    assign(paste("PTSD_CHG_", region[i], "_", "chr", chr.num[j], sep=""), dml.results)
    
    this.dataframe = get(paste("PTSD_CHG_", region[i],"_","chr", chr.num[j], sep=""))
    this.dataframe$region = region[i]
    this.dataframe$type = "CHG"
    
    PTSD.CHG.dataframe.list = append(PTSD.CHG.dataframe.list, 
                                     list(this.dataframe))
  }
  
  region.dataframe = do.call(rbind,PTSD.CHG.dataframe.list)
  assign(paste("PTSD_CHG_", region[i], "_full", sep=""), region.dataframe)
  
}

## import CHH data ----

for (i in 1:length(region)) {
  
  PTSD.CHH.dataframe.list = list()
  
  for (j in 1:length(chr.num)) {
    load(paste0(data_dir, 
                "DSS_Filtered_NoSmooth_Ancestry_Prop_Results_",
                region[i],"_PTSD_", "chr", chr.num[j], "_CHH.Rdata"))
    
    assign(paste("PTSD_CHH_", region[i], "_", "chr", chr.num[j], sep=""), dml.results)
    
    this.dataframe = get(paste("PTSD_CHH_", region[i],"_","chr", chr.num[j], sep=""))
    this.dataframe$region = region[i]
    this.dataframe$type = "CHH"
    
    PTSD.CHH.dataframe.list = append(PTSD.CHH.dataframe.list, 
                                     list(this.dataframe))
  }
  
  region.dataframe = do.call(rbind,PTSD.CHH.dataframe.list)
  assign(paste("PTSD_CHH_", region[i], "_full", sep=""), region.dataframe)
  
}

# import DSS results files ----

load(paste0(data_dir,"AllSamples_DSS_CHG_NewAnnotation_Chromatin.Rdata"))
CHG_chromatin = annoGR2DF(data_annotated)
load(paste0(data_dir,"AllSamples_DSS_CHH_NewAnnotation_Chromatin.Rdata"))
CHH_chromatin = annoGR2DF(data_annotated)

load(paste0(data_dir,"AllSamples_DSS_CHG_NewAnnotation_Genic.Rdata"))
CHG_genic = annoGR2DF(data_annotated)
load(paste0(data_dir,"AllSamples_DSS_CHH_NewAnnotation_Genic.Rdata"))
CHH_genic = annoGR2DF(data_annotated)

load(paste0(data_dir,"AllSamples_DSS_CHG_NewAnnotation.Rdata"))
CHG_annotation = annoGR2DF(data_annotated)
load(paste0(data_dir,"AllSamples_DSS_CHH_NewAnnotation.Rdata"))
CHH_annotation = annoGR2DF(data_annotated)

######################################################################################################
################### -------------------------- Prepare Data ---------------------- ###################
######################################################################################################

CHG_all_annotation = merge(merge(CHG_genic,CHG_annotation,by = c("chr","start","end","width","strand")),
                           CHG_chromatin,by = c("chr","start","end","width","strand"))
CHH_all_annotation = merge(merge(CHH_genic,CHH_annotation,by = c("chr","start","end","width","strand")),
                           CHH_chromatin,by = c("chr","start","end","width","strand"))

CpH_all_annotation = rbind(CHG_all_annotation, CHH_all_annotation)

### Calculate FDR and Bonferroni results ----
PTSD_BLA_full = rbind(PTSD_CHG_BLA_full[,1:4], PTSD_CHH_BLA_full[,1:4])
colnames(PTSD_BLA_full) = c("chr", "start", "stat", "pval")
PTSD_BLA_full$chr = paste0("chr",PTSD_BLA_full$chr)
PTSD_BLA_full$pbon = p.adjust(PTSD_BLA_full$pval, method = "bonferroni")
PTSD_BLA_full$pfdr = p.adjust(PTSD_BLA_full$pval, method = "fdr")
CpH_BLA_all_annotation = merge(CpH_all_annotation, PTSD_BLA_full, by = c("chr", "start"))

PTSD_CNA_full = rbind(PTSD_CHG_CNA_full[,1:4], PTSD_CHH_CNA_full[,1:4])
colnames(PTSD_CNA_full) = c("chr", "start", "stat", "pval")
PTSD_CNA_full$chr = paste0("chr",PTSD_CNA_full$chr)
PTSD_CNA_full$pbon = p.adjust(PTSD_CNA_full$pval, method = "bonferroni")
PTSD_CNA_full$pfdr = p.adjust(PTSD_CNA_full$pval, method = "fdr")
CpH_CNA_all_annotation = merge(CpH_all_annotation, PTSD_CNA_full, by = c("chr", "start"))

PTSD_MNA_full = rbind(PTSD_CHG_MNA_full[,1:4], PTSD_CHH_MNA_full[,1:4])
colnames(PTSD_MNA_full) = c("chr", "start", "stat", "pval")
PTSD_MNA_full$chr = paste0("chr",PTSD_MNA_full$chr)
PTSD_MNA_full$pbon = p.adjust(PTSD_MNA_full$pval, method = "bonferroni")
PTSD_MNA_full$pfdr = p.adjust(PTSD_MNA_full$pval, method = "fdr")
CpH_MNA_all_annotation = merge(CpH_all_annotation, PTSD_MNA_full, by = c("chr", "start"))

PTSD_CA_full = rbind(PTSD_CHG_CA_full[,1:4], PTSD_CHH_CA_full[,1:4])
colnames(PTSD_CA_full) = c("chr", "start", "stat", "pval")
PTSD_CA_full$chr = paste0("chr",PTSD_CA_full$chr)
PTSD_CA_full$pbon = p.adjust(PTSD_CA_full$pval, method = "bonferroni")
PTSD_CA_full$pfdr = p.adjust(PTSD_CA_full$pval, method = "fdr")
CpH_CA_all_annotation = merge(CpH_all_annotation, PTSD_CA_full, by = c("chr", "start"))

PTSD_DG_full = rbind(PTSD_CHG_DG_full[,1:4], PTSD_CHH_DG_full[,1:4])
colnames(PTSD_DG_full) = c("chr", "start", "stat", "pval")
PTSD_DG_full$chr = paste0("chr",PTSD_DG_full$chr)
PTSD_DG_full$pbon = p.adjust(PTSD_DG_full$pval, method = "bonferroni")
PTSD_DG_full$pfdr = p.adjust(PTSD_DG_full$pval, method = "fdr")
CpH_DG_all_annotation = merge(CpH_all_annotation, PTSD_DG_full, by = c("chr", "start"))

PTSD_Sub_full = rbind(PTSD_CHG_Sub_full[,1:4], PTSD_CHH_Sub_full[,1:4])
colnames(PTSD_Sub_full) = c("chr", "start", "stat", "pval")
PTSD_Sub_full$chr = paste0("chr",PTSD_Sub_full$chr)
PTSD_Sub_full$pbon = p.adjust(PTSD_Sub_full$pval, method = "bonferroni")
PTSD_Sub_full$pfdr = p.adjust(PTSD_Sub_full$pval, method = "fdr")
CpH_Sub_all_annotation = merge(CpH_all_annotation, PTSD_Sub_full, by = c("chr", "start"))

### Create the combined dataframe ----
CpH_BLA_all_annotation = CpH_BLA_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_BLA_all_annotation) = c("seqnames", "start", "end", "width", "strand", "BLA",
                                     "nearestGeneName", "distanceToNearestGene", "relativeLocation")

CpH_CNA_all_annotation = CpH_CNA_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_CNA_all_annotation) = c("seqnames", "start", "end", "width", "strand", "CNA",
                                     "nearestGeneName", "distanceToNearestGene", "relativeLocation")

CpH_MNA_all_annotation = CpH_MNA_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_MNA_all_annotation) = c("seqnames", "start", "end", "width", "strand", "MNA",
                                     "nearestGeneName", "distanceToNearestGene", "relativeLocation")

CpH_CA_all_annotation = CpH_CA_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_CA_all_annotation) = c("seqnames", "start", "end", "width", "strand", "CA",
                                    "nearestGeneName", "distanceToNearestGene", "relativeLocation")

CpH_DG_all_annotation = CpH_DG_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_DG_all_annotation) = c("seqnames", "start", "end", "width", "strand", "DG",
                                    "nearestGeneName", "distanceToNearestGene", "relativeLocation")

CpH_Sub_all_annotation = CpH_Sub_all_annotation %>%
  select(chr, start, end, width, strand, pval,
         nearestGeneName, distanceToNearestGene, relativeLocation)
colnames(CpH_Sub_all_annotation) = c("seqnames", "start", "end", "width", "strand", "Sub",
                                     "nearestGeneName", "distanceToNearestGene", "relativeLocation")

PTSD.combined.df = merge(CpH_BLA_all_annotation, CpH_CNA_all_annotation, by = c("seqnames", "start", "end", "width", "strand",
                                                                                "nearestGeneName", "distanceToNearestGene", "relativeLocation"))
PTSD.combined.df = merge(PTSD.combined.df, CpH_MNA_all_annotation, by = c("seqnames", "start", "end", "width", "strand",
                                                                          "nearestGeneName", "distanceToNearestGene", "relativeLocation"))
PTSD.combined.df = merge(PTSD.combined.df, CpH_CA_all_annotation, by = c("seqnames", "start", "end", "width", "strand",
                                                                         "nearestGeneName", "distanceToNearestGene", "relativeLocation"))
PTSD.combined.df = merge(PTSD.combined.df, CpH_DG_all_annotation, by = c("seqnames", "start", "end", "width", "strand",
                                                                         "nearestGeneName", "distanceToNearestGene", "relativeLocation"))
PTSD.combined.df = merge(PTSD.combined.df, CpH_Sub_all_annotation, by = c("seqnames", "start", "end", "width", "strand",
                                                                          "nearestGeneName", "distanceToNearestGene", "relativeLocation"))

PTSD.combined <- GRanges(
  seqnames = PTSD.combined.df$seqnames,
  ranges = IRanges(start = PTSD.combined.df$start, end = PTSD.combined.df$end),
  strand = PTSD.combined.df$strand
)

mcols(PTSD.combined) <- DataFrame(
  BLA = PTSD.combined.df$BLA,
  CNA = PTSD.combined.df$CNA,
  MNA = PTSD.combined.df$MNA,
  CA = PTSD.combined.df$CA,
  DG = PTSD.combined.df$DG,
  Sub = PTSD.combined.df$Sub,
  nearestGeneName = PTSD.combined.df$nearestGeneName,
  distanceToNearestGene = PTSD.combined.df$distanceToNearestGene,
  relativeLocation = PTSD.combined.df$relativeLocation
)

ptsd_gwas_hit_1 = data.frame(snp = c("rs7519147", "rs2777888", "rs4697248", "rs7688962",
                                     "rs10235664", "rs67529088", "rs2532252", "rs2123392"),
                             chr = c(1, 3, 4, 4, 7, 7, 17, 18),
                             pos = c(73994416, 49898000, 21931195, 88281182,
                                     2086814, 104907066, 44257783, 53214865),
                             LD_start = NA,
                             LD_end = NA)

ptsd_gwas_2021 = read.table(paste0(data_dir,"PTSD_GWAS_Hits.txt"))
ptsd_gwas_hit_2 = data.frame(snp = ptsd_gwas_2021$V1,
                             chr = ptsd_gwas_2021$V2,
                             pos = ptsd_gwas_2021$V7,
                             LD_start = NA,
                             LD_end = NA)

ptsd_gwas_ebi = read.csv(paste0(data_dir,"efotraits_EFO_0001358-associations-2022-12-28.csv"))
ptsd_gwas_hit_3 = data.frame(study = ptsd_gwas_ebi$Study.accession,
                             snp = unlist(lapply(strsplit(ptsd_gwas_ebi$Variant.and.risk.allele, "-"), function (x) x[[1]])),
                             chr = as.numeric(unlist(lapply(strsplit(ptsd_gwas_ebi$Location, ":"), function (x) x[[1]]))),
                             pos = as.numeric(unlist(lapply(strsplit(ptsd_gwas_ebi$Location, ":"), function (x) x[[2]]))),
                             LD_start = NA,
                             LD_end = NA)

strsplit(ptsd_gwas_ebi$Mapped.gene, ",") %>% unlist %>% unique %>% sort -> selected.genes

selected.genes %>% length()
ptsd_gwas_hit_3$study %>% unique() %>% length()
ptsd_gwas_hit_3 %>% dim()

ptsd_gwas_hit = rbind(ptsd_gwas_hit_1, ptsd_gwas_hit_2, ptsd_gwas_hit_3[,-1])

bp_range = 500000
ptsd_gwas_hit$LD_start = ptsd_gwas_hit$pos - bp_range
ptsd_gwas_hit$LD_end = ptsd_gwas_hit$pos + bp_range

ptsd_gwas_hit$CHR = paste0("chr", ptsd_gwas_hit$chr)
ptsd_gwas_hit

res_list = list()
for (i in 1:dim(ptsd_gwas_hit)[1]) {
  res_list[[i]] = PTSD.combined[(PTSD.combined@seqnames == ptsd_gwas_hit[i, "CHR"]) &
                                  (PTSD.combined@ranges@start >= ptsd_gwas_hit[i, "LD_start"]) &
                                  ((PTSD.combined@ranges@start + 1) <= ptsd_gwas_hit[i, "LD_end"])]
}
lapply(res_list, function(dat) paste(dat@seqnames, dat@ranges@start, sep = "_")) %>% unlist %>% unique() -> cpgs.selected

res_all = PTSD.combined
res_all$chr_pos = paste(res_all@seqnames, res_all@ranges@start, sep = "_")

res.in = res_all[res_all$chr_pos %in% cpgs.selected,]
res.out = res_all[!(res_all$chr_pos %in% cpgs.selected),]

######################################################################################################
################### ---------------------- Fisher's Exact Test ------------------- ###################
######################################################################################################

p.cutoff = 0.05

sig.idx = c(which(PTSD.combined.df$BLA < p.cutoff),
            which(PTSD.combined.df$CNA < p.cutoff),
            which(PTSD.combined.df$MNA < p.cutoff),
            which(PTSD.combined.df$CA < p.cutoff),
            which(PTSD.combined.df$DG < p.cutoff),
            which(PTSD.combined.df$Sub < p.cutoff))
sig.idx %>% unique() %>% sort() -> sig.idx.sorted
sig_joint_chr_pos = res_all$chr_pos[sig.idx.sorted]

sum(sig_joint_chr_pos %in% res.in$chr_pos) / length(res.in)
sum(sig_joint_chr_pos %in% res.out$chr_pos) / length(res.out)

100*(sum(res.in$chr_pos %in% sig_joint_chr_pos) / sum(!res.in$chr_pos %in% sig_joint_chr_pos))
100*(sum(res.out$chr_pos %in% sig_joint_chr_pos) / sum(!res.out$chr_pos %in% sig_joint_chr_pos))

x <- matrix(c(sum(sig_joint_chr_pos %in% res.in$chr_pos),
              length(res.in),
              sum(sig_joint_chr_pos %in% res.out$chr_pos),
              length(res.out)), byrow = TRUE, 2, 2)

fisher.test(x)
