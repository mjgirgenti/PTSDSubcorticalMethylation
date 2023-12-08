library(Rsubread)

fc <- featureCounts(files = list.files("MapFromRaw", pattern = ".out.bam$", full.names = T), 
                         annot.ext = "/home/tpn5/gibbs/HPC/RefGenome/Homo_sapiens.GRCh38.109.gtf", 
                         isGTFAnnotationFile = T,
                         isPairedEnd = T, 
                         strandSpecific = 2,
                         useMetaFeatures = T,
                         countMultiMappingReads = F,
                         nthreads = 50)
                         
fc$targets <- unlist(lapply(fc$targets, function(x) strsplit(x, "_")[[1]][2]))
colnames(fc$counts) <- fc$targets

write.csv(fc$counts, "GeneExpressionCount.csv", row.names = F)
