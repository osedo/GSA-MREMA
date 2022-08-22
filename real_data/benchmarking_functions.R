
# wrap the mrema functions so that only and SummarizedExperiment object
# and a list of gene sets is required

## change DESeq2.STAT to limma.STAT if limma was used in the runDE() step 

mrema_bench_1DF <- function(se,gs, threshold,  sn = NULL, par = NULL){
  
  GS.MIN.SIZE <- 5
  GS.MAX.SIZE <- 500
  postdata <- tibble("Ensembl" = rownames(rowData(se)), "effect" = rowData(se)$FC, "variance" = (rowData(se)$FC/rowData(se)$DESeq2.STAT)^2)
  postdata_sim <- postdata[complete.cases(postdata),]
  postdata_sim <- postdata_sim[which(postdata$Ensembl %in% unlist(gs)),]
  gs <- lapply(gs, function(s) s[s %in% postdata_sim$Ensembl]) 
  lens <- lengths(gs)
  gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
  
  if(is.null(sn) == FALSE) sn <- which(names(gs) == sn)

  res <- mrema(postdata_sim, gs, DF =1, threshold = threshold, overlap = 0.25, set_number = sn, params = par)
  return(res)
}

mrema_bench_2DF <- function(se,gs, threshold,  sn = NULL, par = NULL){
  GS.MIN.SIZE <- 5
  GS.MAX.SIZE <- 500
  postdata <- tibble("Ensembl" = rownames(rowData(se)), "effect" = rowData(se)$FC, "variance" = (rowData(se)$FC/rowData(se)$DESeq2.STAT)^2)
  postdata_sim <- postdata[complete.cases(postdata),]
  postdata_sim <- postdata_sim[which(postdata$Ensembl %in% unlist(gs)),]
  gs <- lapply(gs, function(s) s[s %in% postdata_sim$Ensembl]) 
  lens <- lengths(gs)
  gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
  
  if(is.null(sn) == FALSE) sn <- which(names(gs) == sn)
  res <- mrema(postdata_sim, gs, DF =2, threshold = threshold, overlap = 0.25)
  return(res)
}

mrema_bench_6DF <- function(se,gs, threshold, sn = NULL, par = NULL){
  
  GS.MIN.SIZE <- 5
  GS.MAX.SIZE <- 500
  postdata <- tibble("Ensembl" = rownames(rowData(se)), "effect" = rowData(se)$FC, "variance" = (rowData(se)$FC/rowData(se)$DESeq2.STAT)^2)
  postdata_sim <- postdata[complete.cases(postdata),]
  postdata_sim <- postdata_sim[which(postdata$Ensembl %in% unlist(gs)),]
  gs <- lapply(gs, function(s) s[s %in% postdata_sim$Ensembl]) 
  lens <- lengths(gs)
  gs <- gs[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
  
  if(is.null(sn) == FALSE) sn <- which(names(gs) == sn)
  
  res <- mrema(postdata_sim, gs, DF =6, threshold = threshold, overlap = 0.25, set_number = sn, params = par)
  return(res)
}


gsea_abs_bench <- function(se,gs){
  res <- sbea_abs(se, method = "gsea", perm = 10, gs = gs)
  results <- res$res.tbl
  results <- tibble::as_data_frame(results)
  results <- results[order(results$GENE.SET),]
  return(results)
}

