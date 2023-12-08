#transfer to bfile
plink --file Girgenti_GS_GDA_082323 --make-bed --out Girgenti_GS_GDA_082323

#step 0: remove samples not in meta data
plink --bfile ../data/Girgenti_GS_GDA_082323 --remove ../data/remove.txt --make-bed --out ../output/step0/Girgenti_GS_GDA_082323_remove

#step 1: delete missing rate > 0.2
plink --bfile ../output/step0/Girgenti_GS_GDA_082323_remove --missing
plink --bfile ../output/step0/Girgenti_GS_GDA_082323_remove --geno 0.2 --make-bed --out ../output/step1/Girgenti_GS_GDA_082323_g02
plink --bfile ../output/step1/Girgenti_GS_GDA_082323_g02 --mind 0.2 --make-bed --out ../output/step1/Girgenti_GS_GDA_082323_g02_i02

#step 2: filter MAF
# Select autosomal SNPs only (i.e., from chromosomes 1 to 22).
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' ../output/step1/Girgenti_GS_GDA_082323_g02_i02.bim > snp_1_22.txt
plink --bfile ../output/step1/Girgenti_GS_GDA_082323_g02_i02 --extract snp_1_22.txt --make-bed --out ../output/step3/Girgenti_GS_GDA_082323_chr122
# Remove SNPs with a low MAF frequency.
plink --bfile ../output/step3/Girgenti_GS_GDA_082323_chr122 --maf 0.01 --make-bed --out ../output/step3/Girgenti_GS_GDA_082323_maf005

#step 3: # Delete SNPs which are not in Hardy-Weinberg equilibrium (HWE).
plink --bfile ../output/step3/Girgenti_GS_GDA_082323_maf005 --hwe 1e-6 --make-bed --out ../output/step4/Girgenti_GS_GDA_082323_hwe1
plink --bfile ../output/step4/Girgenti_GS_GDA_082323_hwe1 --hwe 1e-10 --hwe-all --make-bed --out ../output/step4/Girgenti_GS_GDA_082323_hwe2

#step 4: check datasets for cryptic relatedness
plink --bfile ../output/step4/Girgenti_GS_GDA_082323_hwe2 --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2
awk '{ if ($8 >0.9) print $0 }' pihat_min0.2.genome>zoom_pihat.genome
plink --bfile ../output/step4/Girgenti_GS_GDA_082323_hwe2 --filter-founders --make-bed --out ../output/step6/Girgenti_GS_GDA_082323_relate
plink --bfile ../output/step6/Girgenti_GS_GDA_082323_relate --extract indepSNP.prune.in --genome --min 0.2 --out pihat_min0.2_in_founders
plink --bfile ../output/step6/Girgenti_GS_GDA_082323_relate --missing
vi 0.2_low_call_rate_pihat.txt
plink --bfile ../output/step6/Girgenti_GS_GDA_082323_relate --remove 0.2_low_call_rate_pihat.txt --make-bed --out ../output/final/Girgenti_GS_GDA_082323_QC



