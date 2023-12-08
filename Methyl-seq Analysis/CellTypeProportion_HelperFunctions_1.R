cell_type_fraction <- function(DNAm_gene_level_matrix, DNAmRef, disease = FALSE){
  #Returns the estimated cell fraction matrix, barplot, and boxplot
  require(ggplot2)
  require(tidyr)
  require(patchwork)
  require(ggsci)
  require(EpiSCORE)
  #compute estimated cell type fractions
  est_cell_fractions <- wRPC(DNAm_gene_level_matrix, DNAmRef)$estF
  est_cell_fractions_tb <- as_tibble(est_cell_fractions)
  est_cell_fractions_tb$sample_id <- rownames(est_cell_fractions)
  if(!is.character(disease)){
    est_cell_fractions_tb <- est_cell_fractions_tb %>% 
      gather(colnames(est_cell_fractions_tb)[-ncol(est_cell_fractions_tb)], key = "cell_type", value = "fractions")
    bar_plot <- ggplot(est_cell_fractions_tb, 
                       aes(fill = cell_type, y = fractions, x = sample_id)) + 
      geom_bar(position = "stack", stat = "identity", width = 1) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      xlab("Samples") + ylab("Cell Type Proportion") + 
      labs(fill = "Cell Type") + 
      scale_fill_npg()
    
    box_plot <- ggplot(est_cell_fractions_tb, 
                       aes(y = fractions, x = cell_type, fill = cell_type)) + 
      geom_boxplot() +
      theme_minimal() +
      theme(legend.position = "none") +
      xlab("") + ylab("Cell Type Proportion") + 
      scale_fill_npg()
  }else{
    est_cf <- data.frame(est_cell_fractions)
    phenotype <- vapply(rownames(est_cf), function(x){
      ifelse(grepl(disease, x), "Pathological", "Normal")
    }, character(1))
    est_cf$Phenotype <- phenotype
    est_cf$sample_id <- rownames(est_cell_fractions)
    est_cf_plot <- est_cf %>% 
      gather(colnames(est_cf)[1:(ncol(est_cf) - 2)], key = "cell_type", value = "est_cf")
    box_plot <- ggplot(est_cf_plot, 
                       aes(y = est_cf, x = cell_type, fill = Phenotype)) + 
      geom_boxplot() +
      theme_minimal() +
      xlab("") + ylab("Cell Type Proportion") + 
      scale_fill_npg()
    bar_plot <- ggplot(est_cf_plot, 
                       aes(fill = cell_type, y = est_cf, x = sample_id)) + 
      facet_wrap(~Phenotype, scales = "free") + 
      geom_bar(position = "stack", stat = "identity", width = 1) +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank()) +
      xlab("Samples") + ylab("Cell Type Proportion") + 
      labs(fill = "Cell Type") + 
      scale_fill_npg()
  }
  
  list(est_cell_fractions = est_cell_fractions,
       bar_plot = bar_plot,
       box_plot = box_plot)
}

CellDMC <- function(beta.m, pheno.v, frac.m, 
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    sort = FALSE, mc.cores = 1) {
  ### check input
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v)) 
    stop("Number of columns of beta.m should equal to length of pheno.v!")
  if (ncol(beta.m) != nrow(frac.m)) 
    stop("Number of columns of beta.m should equal to number of rows of frac.m!")
  if (length(colnames(frac.m)) != ncol(frac.m)) 
    stop("Pls assign correct names of cell-type to frac.m")
  message("using our own dmc")
  
  ### check whether input is beta value matrix
  is.beta <- ((min(beta.m) >= 0) & (max(beta.m) <= 1))
  
  ### guess factor input
  if (nlevels(factor(pheno.v)) == 2) {
    message("Binary phenotype detected. Predicted change will be 1 - 0.")
    pheno.v <- factor(pheno.v)
  }
  if (!is.factor(pheno.v) & !is.character(pheno.v)) 
    message("pheno.v is not factor or character. Treating as continuous variables.")
  
  ### Fit model
  design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
  if (Matrix::rankMatrix(design) < ncol(design)) {
    stop("The design matrix is not full ranked.\nThis means that you coundn't make inference for all cell-types in your fraction matrix.
         This is usally casued by fractions of a cell-type of one pheno type are all 0 or some fractions in one pheno type are paralle to that of another cell-type.
         You might use which(colSums(model.matrix(~ frac.m + pheno.v:frac.m)[, -1]) == 0) to find the cell type.")
  }
  
  
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  IntNames.v <- str_c(colnames(frac.m), "Pheno")
  colnames(design)[(1 + ncol(frac.m)):(2*ncol(frac.m))] <- IntNames.v 
  
  ### fit linear model for each CpG
  allCoe.m <- do.call(rbind, mclapply(seq_len(nrow(beta.m)), function(i) {
    beta.v <- beta.m[i, ]
    ### model
    Int.o <- lm(beta.v ~ .-1, data = data.frame(design))
    
    ### get coe
    IntCoe.m <- summary(Int.o)$coe[IntNames.v, ]
    IntCoe.v <- unlist(apply(IntCoe.m, 1, function(x) list(x)))
    
    names(IntCoe.v) <- NULL
    return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  return(allCoe.m)
  ### extract coefficients for each cell-type
  # coe.ld <- lapply(seq_len(ncol(frac.m)), function(j) {
  #   idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
  #   tmp.m <- allCoe.m[, idx]
  #   tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
  #   if (is.beta) { 
  #     tmp.m[which(tmp.m[,1] > 1),1] <- 1
  #     tmp.m[which(tmp.m[,1] < -1),1] <- -1
  #   }  ### if input is a beta values matrix, bound the estimated changes
  #   
  #   colnames(tmp.m) <- c("Estimate", "SE", "t", "p", "adjP")
  #   rownames(tmp.m) <- rownames(beta.m)
  #   return(data.frame(tmp.m))
  # })
  # names(coe.ld) <- colnames(frac.m)
  # 
  # ### get dmct matrix
  # dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
  # dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  # dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  # dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  # dmct.m <- cbind(dmc.v, dmct.m)
  # colnames(dmct.m) <- c("DMC", colnames(frac.m))
  # rownames(dmct.m) <- rownames(beta.m)
  # 
  # if(sort) coe.ld <- lapply(coe.ld, function(x) x[order(x$p),] )
  # 
  # return(list(dmct = dmct.m, coe = coe.ld))
}

CellDMC_rlm <- function(beta.m, pheno.v, frac.m, 
                    adjPMethod = "fdr", adjPThresh = 0.05, cov.mod = NULL, 
                    sort = FALSE, mc.cores = 1) {
  require(stringr)
  require(sfsmisc)
  require(pbmcapply)
  ### check input
  if (sum(is.na(pheno.v)) > 0) stop("No NA allowed in pheno.v!")
  if (ncol(beta.m) != length(pheno.v)) 
    stop("Number of columns of beta.m should equal to length of pheno.v!")
  if (ncol(beta.m) != nrow(frac.m)) 
    stop("Number of columns of beta.m should equal to number of rows of frac.m!")
  if (length(colnames(frac.m)) != ncol(frac.m)) 
    stop("Pls assign correct names of cell-type to frac.m")
  
  ### check whether input is beta value matrix
  is.beta <- ((min(beta.m) >= 0) & (max(beta.m) <= 1))
  
  ### guess factor input
  if (nlevels(factor(pheno.v)) == 2) {
    message("Binary phenotype detected. Predicted change will be 1 - 0.")
    pheno.v <- factor(pheno.v)
  }
  if (!is.factor(pheno.v) & !is.character(pheno.v)) 
    message("pheno.v is not factor or character. Treating as continuous variables.")
  
  ### Fit model
  design <- model.matrix(~ frac.m + pheno.v:frac.m)[, -1]
  if (Matrix::rankMatrix(design) < ncol(design)) {
    stop("The design matrix is not full ranked.\nThis means that you coundn't make inference for all cell-types in your fraction matrix.
         This is usally casued by fractions of a cell-type of one pheno type are all 0 or some fractions in one pheno type are paralle to that of another cell-type.
         You might use which(colSums(model.matrix(~ frac.m + pheno.v:frac.m)[, -1]) == 0) to find the cell type.")
  }
  
  
  if (!is.null(cov.mod)) design <- cbind(design, cov.mod[, -1])
  IntNames.v <- str_c(colnames(frac.m), "Pheno")
  colnames(design)[(1 + ncol(frac.m)):(2*ncol(frac.m))] <- IntNames.v
  
  ### fit robust linear model for each CpG
  print("----------Fitting Robust Linear Regression Models----------")
  allCoe.m <- do.call(rbind, pbmclapply(seq_len(nrow(beta.m)), function(i) {
    beta.v <- beta.m[i, ]
    ### model
    Int.o <- rlm(beta.v ~ .-1, data = data.frame(design), model = F)
    ### get p-values of coe and coe
    p_values <- vapply(rownames(summary(Int.o)$coe), function(name_){
      f.robftest(Int.o, var = name_)$p.value
    }, numeric(1))
    IntCoe.m <- cbind(summary(Int.o)$coe, p_values)[IntNames.v, ]
    IntCoe.v <- unlist(apply(IntCoe.m, 1, function(x) list(x)))
    
    names(IntCoe.v) <- NULL
    return(IntCoe.v)
  }, mc.preschedule = TRUE, mc.cores = mc.cores, mc.allow.recursive = TRUE))
  
  print("----------Extracting Coefficients----------")
  ### extract coefficients for each cell-type
  coe.ld <- lapply(seq_len(ncol(frac.m)), function(j) {
    idx <- ((j - 1)*4 + 1):((j - 1)*4 + 4)
    tmp.m <- allCoe.m[, idx]
    tmp.m <- cbind(tmp.m, p.adjust(tmp.m[, 4], method = adjPMethod))
    if (is.beta) { 
      tmp.m[which(tmp.m[,1] > 1),1] <- 1
      tmp.m[which(tmp.m[,1] < -1),1] <- -1
    }  ### if input is a beta values matrix, bound the estimated changes
    
    colnames(tmp.m) <- c("Estimate", "SE", "t", "p", "adjP")
    rownames(tmp.m) <- rownames(beta.m)
    return(data.frame(tmp.m))
  })
  names(coe.ld) <- colnames(frac.m)
  
  print("----------Extracting DMCTs----------")
  ### get dmct matrix
  dmct.m <- matrix(rep(0, ncol(frac.m)*nrow(beta.m)), ncol = ncol(frac.m))
  dmct.idx <- which(sapply(coe.ld, "[[", "adjP") < adjPThresh)
  dmct.m[dmct.idx] <- sign(sapply(coe.ld, "[[", "Estimate")[dmct.idx])
  dmc.v <- ifelse(rowAlls(dmct.m == 0), 0, 1)
  dmct.m <- cbind(dmc.v, dmct.m)
  colnames(dmct.m) <- c("DMC", colnames(frac.m))
  rownames(dmct.m) <- rownames(beta.m)
  
  if(sort) coe.ld <- lapply(coe.ld, function(x) x[order(x$p),] )
  
  return(list(dmct = dmct.m, coe = coe.ld))
}

.enrichment_prep_df <- function(df, showTerms, orderBy) {
  
  if(is.null(showTerms)) {
    showTerms = nrow(df)
  } else if(!is.numeric(showTerms)) {
    stop(paste0("showTerms '", showTerms, "' is invalid."))
  }
  
  Annotated <- as.numeric(sub("^\\d+/", "", as.character(df$Overlap)))
  Significant <- as.numeric(sub("/\\d+$", "", as.character(df$Overlap)))
  
  # Build data frame
  df <- cbind(df, data.frame(Annotated = Annotated, Significant = Significant,
                             stringsAsFactors = FALSE))
  
  # Order data frame (P.value or Combined.Score)
  if(orderBy == "Combined.Score") {
    idx <- order(df$Combined.Score, decreasing = TRUE)
  } else {
    idx <- order(df$P.value, decreasing = FALSE)
  }
  df <- df[idx,]
  
  # Subset to selected number of terms
  if(showTerms <= nrow(df)) {
    df <- df[1:showTerms,]
  }
  
  return(df)
}

plotEnrich <- function(df, showTerms = 20, numChar = 40, y = "Count", orderBy = "P.value",
                       xlab = NULL, ylab = NULL, title = NULL) {
  if(!is.data.frame(df)) stop("df is malformed - must be a data frame")
  if(!is.numeric(numChar)) {
    stop(paste0("numChar '", numChar, "' is invalid."))
  }
  
  df <- .enrichment_prep_df(df, showTerms, orderBy)
  
  # Create trimmed name (as seen in topGO)
  shortName <- paste(substr(df$Term, 1, numChar),
                     ifelse(nchar(df$Term) > numChar, '...', ''), sep = '')
  df$shortName = shortName
  df$shortName <- factor(df$shortName, levels = rev(unique(df$shortName)))
  df$Ratio <- df$Significant/df$Annotated
  
  # Define fill variable (P.value or Combined.Score)
  if(orderBy == "Combined.Score") {
    fill <- "Combined.Score"
  } else {
    fill <- "P.value"
  }
  
  # Define y variable (Count or Ratio)
  if(y != "Ratio") {
    y <- "Significant"
  }
  
  # Define variable mapping
  map <- aes_string(x = "shortName", y = y, fill = fill)
  
  # Define labels
  if(is.null(xlab)) {
    xlab <- "Enriched terms"
  }
  
  if(is.null(ylab)) {
    if(y == "Ratio") {
      ylab <- "Gene ratio"
    } else {
      ylab <- "Gene count"
    }
  }
  
  if(is.null(title)) {
    title <- "Enrichment analysis by Enrichr"
  }
  
  # Make the ggplot
  p <- ggplot(df, map) + 
    geom_bar(stat = "identity", width = 0.4, color = "black") + 
    coord_flip() + 
    theme_minimal()
  
  if(orderBy == "Combined.Score") {
    p <- p + scale_fill_continuous(low = "blue", high = "red") +
      guides(fill = guide_colorbar(title = "Combined Score", reverse = FALSE))
  } else {
    p <- p + scale_fill_continuous(high = "#4DBBD5B2", low = "#DC0000B2") +
      guides(fill = guide_colorbar(title = "P value", reverse = TRUE))
  }
  
  # Adjust theme components
  p <- p + theme(axis.text.x = element_text(colour = "black", vjust = 1),
                 axis.text.y = element_text(colour = "black", hjust = 1),
                 axis.title = element_text(color = "black", margin = margin(10, 5, 0, 0)),
                 axis.title.y = element_text(angle = 90))
  
  p <- p + xlab(xlab) + ylab(ylab) + ggtitle(title)
  
  return(p)
}

epiSCORE_construct_references <- function(sc_count_matrix, ct.idx, names_ct.v, MSS=NULL){
  #normalizing sc_count_matrix
  message("normalizing the count matrix")
  maximum_read_count <- max(colSums(sc_count_matrix))
  #step I, normalizing the read count matrix and constructing the mRNA reference matrix
  sc_count_matrix <- pbapply(sc_count_matrix, 2, function(cell){
    log(cell*(maximum_read_count/sum(cell)) + 1, base = 2)
  })
  gc(verbose = F)
  message("constructing the mRNA reference matrix")
  #construct reference expression matrix
  mRNA_ref_ad.m <- ConstExpRef(exp.m = sc_count_matrix, celltype.idx = ct.idx, namesCellT.v = names_ct.v, 
                               markspecTH.v = MSS)
  gc(verbose = F)
  RNArefBrain_ROSMAP.m <- mRNA_ref_ad.m$ref
  message("imputing the DNAm matrix")
  #DNAm reference matrix, merging to a final one from 2 databases
  DNAm_ref_ROSMAP_1.m <- ImputeDNAmRef(refexp.m = RNArefBrain_ROSMAP.m$med, db = "SCM2", geneID = "SYMBOL")
  DNAm_ref_ROSMAP_2.m <- ImputeDNAmRef(refexp.m = RNArefBrain_ROSMAP.m$med, db = "RMAP", geneID = "SYMBOL")
  DNAm_ref_sc.m <- ConstMergedDNAmRef(DNAm_ref_ROSMAP_1.m, DNAm_ref_ROSMAP_2.m)
  list(mRNA_ref = mRNA_ref_ad.m, DNAm_ref = DNAm_ref_sc.m)
}
