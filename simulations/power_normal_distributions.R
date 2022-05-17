#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args <- unlist(lapply(strsplit(args, ","), as.numeric))

source("../mrema.R")
source("../packages.R")
source("../sbea_gsea_abs.R")
source("sim_distribution.R")
j <- read.csv("BRCA_sim_param.csv")

n_samples <- args[1]
sd_in <- args[2]
sd_out <- args[3]
n_sims <- args[4]

T1.1_1DF_res <- c() 
T1.3_1DF_res <- c()
T1.5_1DF_res <- c()

res_gsea <- c()
res_ora <- c()
res_safe <- c()

res_ora_1.1 <- c()
res_ora_1.3 <- c()
res_ora_1.5 <- c()

for(i in 1:n_sims){
  
  
  fc <- 2^rhalfnorm(100, sd2theta(sd_in))
  simulation <- simulation_dist(s = 100, p = 100, N = n_samples, beta = 0.01, gamma = 1, theta = 1, foldchange = fc, inset_fc = sd_out, upreg = 0.5)
  
  # filter out genes not in list of gene sets and lowly expressed genes
  igenes <- intersect(rownames(assay(simulation$sum_exp_raw_count)), unique(unlist(simulation$raw.gs)))
  if(!length(igenes)) stop("Expression dataset (se)", " and ", "gene sets (gs) have no gene IDs in common")
  se <- assay(simulation$sum_exp_raw_count)[igenes,]
  group <- colData(simulation$sum_exp_raw_count)$GROUP
  keep <- filterByExpr(se, group = group)
  se <- se[keep,]
  gene_sets <- lapply(simulation$raw.gs, function(s) s[s %in% rownames(se)]) 
  
  ### use Deseq2 to get group effect estimate and it's standard error
  dea_raw_count<- DESeqDataSetFromMatrix(countData = se, colData = colData(simulation$sum_exp_raw_count), design = ~ GROUP)
  dds <- DESeq(dea_raw_count)
  res <- results(dds)
  postdata <- tibble("Ensembl" = rownames(res), "effect" = res[,2], "variance" = res[,3]^2)
  postdata_sim <- postdata[complete.cases(postdata),]
  
  
  T1.5_1DF <- mrema(postdata_sim, gene_sets, DF =1, threshold = 1.5, set_number = 100, CI_interval = FALSE)
  T1.5_1DF_res <- rbind(T1.5_1DF, T1.5_1DF_res)
  T1.3_1DF <- mrema(postdata_sim, gene_sets, DF =1, threshold = 1.3, set_number = 100, CI_interval = FALSE)
  T1.3_1DF_res <- rbind(T1.3_1DF, T1.3_1DF_res)
  T1.1_1DF <- mrema(postdata_sim, gene_sets, DF =1, threshold = 1.1, set_number = 100, CI_interval = FALSE)
  T1.1_1DF_res <- rbind(T1.1_1DF, T1.1_1DF_res)

  
  
   de_sim <- deAna(simulation$sum_exp_raw_count, de.method = "DESeq2", padj.method = "BH")
   
   safe_res <- sbea(de_sim, method = "safe", perm = 1000, gs = simulation$raw.gs)
   res_safe <- rbind(res_safe, safe_res$res.tbl[which(safe_res$res.tbl$GENE.SET == "gs100"),])
   ora_res <- sbea(de_sim, method = "ora", perm = 0, gs = simulation$raw.gs, padj.method = "BH")
   res_ora <- rbind(res_ora, ora_res$res.tbl[which(ora_res$res.tbl$GENE.SET == "gs100"),])
   gsea_res <- sbea_abs(de_sim, method = "gsea", perm = 1000, gs = simulation$raw.gs, padj.method = "BH")
   res_gsea <- rbind(res_gsea, gsea_res$res.tbl[which(gsea_res$res.tbl$GENE.SET == "gs100"),])
   rm(.Random.seed, envir=globalenv())
   
   
   res_1.5 <- results(dds, lfcThreshold = log2(1.5))
   rowData_1.5 <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res_1.5$log2FoldChange, "ADJ.PVAL" = res_1.5$padj)
   res_1.3 <- results(dds, lfcThreshold = log2(1.3))
   rowData_1.3 <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res_1.3$log2FoldChange, "ADJ.PVAL" = res_1.3$padj)
   res_1.1 <- results(dds, lfcThreshold = log2(1.1))
   rowData_1.1 <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res_1.1$log2FoldChange, "ADJ.PVAL" = res_1.1$padj)
   
   sumexp_1.1 <- SummarizedExperiment(assays=SimpleList(counts=se), colData=colData(simulation$sum_exp_raw_count), rowData=rowData_1.1)
   sumexp_1.3 <- SummarizedExperiment(assays=SimpleList(counts=se), colData=colData(simulation$sum_exp_raw_count), rowData=rowData_1.3)
   sumexp_1.5 <- SummarizedExperiment(assays=SimpleList(counts=se), colData=colData(simulation$sum_exp_raw_count), rowData=rowData_1.5)
   
   ora_res_1.1 <- sbea(sumexp_1.1, method = "ora", perm = 0, gs = simulation$raw.gs, padj.method = "BH")
   res_ora_1.1 <- rbind(res_ora_1.1, ora_res_1.1$res.tbl[which(ora_res_1.1$res.tbl$GENE.SET == "gs100"),])
   ora_res_1.3 <- sbea(sumexp_1.3, method = "ora", perm = 0, gs = simulation$raw.gs, padj.method = "BH")
   res_ora_1.3 <- rbind(res_ora_1.3, ora_res_1.3$res.tbl[which(ora_res_1.3$res.tbl$GENE.SET == "gs100"),])
   ora_res_1.5 <- sbea(sumexp_1.5, method = "ora", perm = 0, gs = simulation$raw.gs, padj.method = "BH")
   res_ora_1.5 <- rbind(res_ora_1.5, ora_res_1.5$res.tbl[which(ora_res_1.5$res.tbl$GENE.SET == "gs100"),])
   
   print(i)
}




write.csv(T1.5_1DF_res, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_T1.5_1DF.csv"), row.names = FALSE, quote = FALSE)
write.csv(T1.1_1DF_res, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_T1.1_1DF.csv"), row.names = FALSE, quote = FALSE)
write.csv(T1.3_1DF_res, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_T1.3_1DF.csv"), row.names = FALSE, quote = FALSE)

write.csv(res_gsea, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_gsea.csv"), row.names = FALSE, quote = FALSE)
write.csv(res_ora, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_ora.csv"), row.names = FALSE, quote = FALSE )
write.csv(res_safe, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_safe.csv"), row.names = FALSE, quote = FALSE )

write.csv(res_ora_1.1, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_ora_1.1.csv"), row.names = FALSE, quote = FALSE)
write.csv(res_ora_1.3, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_ora_1.3.csv"), row.names = FALSE, quote = FALSE)
write.csv(res_ora_1.5, file = paste0("Power_sims_", args[5], "_N", n_samples, "_SDinset_", sd_in, "_SDoutset_", sd_out, "_ora_1.5.csv"), row.names = FALSE, quote = FALSE)
