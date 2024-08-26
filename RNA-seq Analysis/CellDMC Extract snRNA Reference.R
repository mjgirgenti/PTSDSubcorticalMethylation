library(SingleCellExperiment)
library(dplyr)
library(forcats)
# Tran et al. Neuron data (2021)

# Link provided in: https://github.com/LieberInstitute/10xPilot_snRNAseq-human

# Download from: https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_NAc-n8_tran-etal.rda

# dir.create("public_data/Tran")
# system("wget https://libd-snrnaseq-pilot.s3.us-east-2.amazonaws.com/SCE_NAc-n8_tran-etal.rda -P public_data/Tran")
load(".../SCE_AMY-n5_tran-etal.rda", verbose = TRUE)

# Collapse factors
sce.amy.tran$cellType <- fct_relabel(sce.amy.tran$cellType, ~ gsub("_.*", "",.x))
# Remove lowNT columns
sce.amy.tran <- sce.amy.tran[, !grepl("drop", sce.amy.tran$cellType)]
# Remove cell types with less then 100 cells
selected_celltype <- names(table(sce.amy.tran$cellType))[table(sce.amy.tran$cellType) >= 100]
sce.amy.tran <- sce.amy.tran[, sce.amy.tran$cellType %in% selected_celltype]
sce.amy.tran$cellType <- droplevels(sce.amy.tran$cellType)

# Take up to 1000 cells per cell types with highest detected values
top_cells <- colData(sce.amy.tran) %>% 
  as.data.frame() %>%
  group_by(cellType) %>%
  slice_max(order_by = detected, n = 1000)
subset_sce.amy.tran <- sce.amy.tran[,sce.amy.tran$Barcode %in% top_cells$Barcode]

# Remove genes with 0 counts
subset_sce.amy.tran <- subset_sce.amy.tran[rowSums(counts(subset_sce.amy.tran)) != 0,]
counts(subset_sce.amy.tran) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>% 
  write.table(., file = "sc_AMG.txt", sep = "\t", quote = F, row.names = F, col.names = c("gene", as.character(subset_sce.amy.tran$cellType)))

load(".../SCE_HPC-n3_tran-etal.rda", verbose = TRUE)

# Collapse factors
sce.hpc.tran$cellType <- fct_relabel(sce.hpc.tran$cellType, ~ gsub("_.*", "",.x))
# Remove lowNT columns
sce.hpc.tran <- sce.hpc.tran[, !grepl("drop", sce.hpc.tran$cellType)]
# Remove cell types with less then 100 cells
selected_celltype <- names(table(sce.hpc.tran$cellType))[table(sce.hpc.tran$cellType) >= 100]
sce.hpc.tran <- sce.hpc.tran[, sce.hpc.tran$cellType %in% selected_celltype]
sce.hpc.tran$cellType <- droplevels(sce.hpc.tran$cellType)

# Take up to 1000 cells per cell types with highest detected values
top_cells <- colData(sce.hpc.tran) %>% 
  as.data.frame() %>%
  group_by(cellType) %>%
  slice_max(order_by = detected, n = 1000)
subset_sce.hpc.tran <- sce.hpc.tran[,sce.hpc.tran$Barcode %in% top_cells$Barcode]

# Remove genes with 0 counts
subset_sce.hpc.tran <- subset_sce.hpc.tran[rowSums(counts(subset_sce.hpc.tran)) != 0,]
counts(subset_sce.hpc.tran) %>%
  as.matrix() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "gene") %>% 
  write.table(., file = "sc_HPC.txt", sep = "\t", quote = F, row.names = F, col.names = c("gene", as.character(subset_sce.hpc.tran$cellType)))

