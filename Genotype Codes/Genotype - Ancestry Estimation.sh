
###  step 1: merge studies cohort with reference panel (e.g., 1000G) and perform pruning to reduce total number of SNPs 
plink --bfile Girgenti_GS_GDA_082323_QC --bmerge 1000g_ovp.bed 1000g_ovp.bim 1000g_ovp.fam --make-bed --out Girgenti_1000g_AMRremoved
plink --bfile Girgenti_1000g_merged --make-bed --remove remove_AMR.txt --out Girgenti_1000g_AMRremoved
plink --bfile Girgenti_1000g_AMRremoved --indep-pairwise 5000 1000 0.1 --out pruned5000
plink --bfile Girgenti_1000g_AMRremoved --extract pruned5000.prune.in --make-bed --out Girgenti_1000g_AMRremoved_13w


###  step 2: perform PCA to the merged and pruned bfiles
plink --bfile Girgenti_1000g_AMRremoved_13w --out pca_Girgenti_1000g_AMRremoved_13 --pca 10


###  step 3: According to the results of PCA, choose (AFR, EUR, SAS) as reference ancestry groups and perform admixture 
plink --bfile Girgenti_1000g_AMRremoved_13w --make-bed --out Girgenti_1000g_AMREASremoved_13w --remove remove_EAS.txt
admixture Girgenti_1000g_AMREASremoved_13w.bed 3