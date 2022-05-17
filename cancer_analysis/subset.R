#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args <- unlist(lapply(strsplit(args, ","), as.numeric))

source("packages.R")
source("mrema.R")
load("~/data/BLCA_SummarizedExperiment.RData")
load("kegg.gs.RData")

womrema_n_results <- c()
mrema_n_results <- c()
ora_n_results <- c()
gsea_n_results <- c()
camera_n_results <- c()

n_cases <- args[1]

for(i in 1:args[2]){
  subset_samples <- sample(unique(colData(SummarizedExperiment)$BLOCK), n_cases, replace = FALSE)
  subset_colData <- colData(SummarizedExperiment)[which(colData(SummarizedExperiment)$BLOCK %in% subset_samples),]
  group <- subset_colData$GROUP
  group[which(group != "Solid Tissue Normal")] <- 1
  group[which(group == "Solid Tissue Normal")] <- 0
  subset_colData$GROUP <- group
  se <- assay(SummarizedExperiment)
  se_subset <- se[,which(colData(SummarizedExperiment)$BLOCK %in% subset_samples)]
  
  SummarizedExperiment_subset <- SummarizedExperiment(assays= SimpleList(counts=se_subset), colData=subset_colData, rowData= DataFrame(row.names = substr(rownames(se_subset),1,15)))
  
  
  
  ## run DESeq2 
  dea_raw_count<- DESeqDataSetFromMatrix(countData = se_subset, colData = subset_colData, design = ~ BLOCK + GROUP)
  dea_raw_count$GROUP <- relevel(as.factor(dea_raw_count$GROUP), ref = "0")
  dds <- DESeq(dea_raw_count)
  res <- results(dds)
  res <- res[complete.cases(res),]
  postdata <- tibble("Ensembl" = substr(rownames(res),1,15), "effect" = res[,2], "variance" = res[,3]^2)
  postdata <- postdata[complete.cases(postdata),]
  postdata <- postdata[which(postdata$Ensembl %in% unlist(kegg.gs)),]
  kegg.gs.subset <- lapply(kegg.gs, function(s) s[s %in% postdata$Ensembl])
  
  
  
  womrema_1.5 <- mrema(postdata, kegg.gs.subset, DF = 1, threshold = 1.5, set_number = 121)
  womrema_n_results <- rbind(womrema_n_results, womrema_1.5)
  
  mrema_1.5 <- mrema(postdata, kegg.gs.subset, DF = 7, threshold = 1.5, set_number = 121)
  mrema_n_results <- rbind(mrema_n_results, mrema_1.5)
  
  subset_de <- deAna(SummarizedExperiment_subset, grp = "GROUP", blk = "BLOCK", de.method = "DESeq2", padj.method = "BH")
  
  blca_ora <- sbea(subset_de, kegg.gs, method = "ora", perm = 0, padj.method = "BH")
  ORA_subset <- blca_ora$res.tbl[which(blca_ora$res.tbl$GENE.SET == "hsa04022_cGMP-PKG_signaling_pathway"),]
  ora_n_results <- rbind(ora_n_results, ORA_subset)
  
  blca_gsea <- sbea(subset_de, kegg.gs, method = "gsea", perm = 1000, padj.method = "BH")
  GSEA_subset <- blca_gsea$res.tbl[which(blca_gsea$res.tbl$GENE.SET == "hsa04020_Calcium_signaling_pathway"),]
  gsea_n_results <- rbind(gsea_n_results, GSEA_subset)
  
  blca_camera <- sbea(subset_de, kegg.gs, method = "camera", perm = 1000, padj.method = "BH")
  CAMERA_subset <- blca_camera$res.tbl[which(blca_camera$res.tbl$GENE.SET == "hsa04020_Calcium_signaling_pathway"),]
  camera_n_results <- rbind(camera_n_results, CAMERA_subset)
  rm(.Random.seed, envir=globalenv())
  
}


write.csv(womrema_n_results, file = paste0("BLCA_Subset_N_samples_", args[1] ,"_N_subs_",args[2],"_womrema.csv") , quote = FALSE, row.names = FALSE)
write.csv(gsea_n_results, file = paste0("BLCA_Subset_N_samples_", args[1] ,"_N_subs_",args[2],"_gsea.csv"),  quote = FALSE, row.names = FALSE)
write.csv(camera_n_results, file = paste0("BLCA_Subset_N_samples_", args[1] ,"_N_subs_",args[2],"_camera.csv") , quote = FALSE, row.names = FALSE)
write.csv(mrema_n_results, file = paste0("BLCA_Subset_N_samples_", args[1] ,"_N_subs_",args[2],"_mrema.csv") , quote = FALSE, row.names = FALSE)
write.csv(ora_n_results, file = paste0("BLCA_Subset_N_samples_", args[1] ,"_N_subs_",args[2],"_ora_cGMP.csv") , quote = FALSE, row.names = FALSE)


