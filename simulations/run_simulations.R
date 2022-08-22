#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
args <- unlist(lapply(strsplit(args, ","), as.numeric))

source("../packages.R")
source("simulation_power.R")
source("sbea_gsea_abs.R")
j <- read.csv("BRCA_sim_param.csv")
library(EnrichmentBrowser)
library(GSEABenchmarkeR)
n <- args[1]
outset <- args[2]
inset <- args[3]
n_sims <- args[4]
array <- args[5]

gsa <- vector(mode = "list", len = n_sims)
start <- Sys.time()

for(i in 1:n_sims){
  # simulate distribution in gene set
  lfc_inset <- rnorm(100, 0, sd = inset)
  # simulate distribution in background...background is split into 99 gene sets of 100 genes to allow all methods to run.
  lfc_outset <- sample(rnorm(9900, 0, sd = outset))
  simulation <- simulation_power(inset_LFC = lfc_inset, outset_LFC = lfc_outset, n_samples = n, disp = j$disp_para, mean = j$mean_para, outset_split = 100)

  dea_raw_count<- DESeqDataSetFromMatrix(countData = assay(simulation$sum_exp_raw_count), colData = colData(simulation$sum_exp_raw_count), design = ~ GROUP)
  dds <- DESeq(dea_raw_count)
  res <- results(dds)
  
  
  # run mrema
  postdata <- tibble("Ensembl" = rownames(res), "effect" = res[,2], "variance" = res[,3]^2)
  postdata_sim <- postdata[complete.cases(postdata),]
  
  gsa$mrema_1DF_T1.5[[i]] <- mrema(postdata_sim, simulation$raw.gs, DF =1, threshold = 1.5, overlap = 0.25)
  gsa$mrema_2DF_T1.5[[i]] <- mrema(postdata_sim, simulation$raw.gs, DF= 2, threshold = 1.5, overlap = 0.25)
  gsa$mrema_6DF_T1.5[[i]] <- mrema(postdata_sim, simulation$raw.gs, DF =6, threshold = 1.5, overlap = 0.25)
  
  
  #run established methods
  rowData <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res$log2FoldChange,"LFC.SE" = res$lfcSE, "ADJ.PVAL" = res$padj)
  de_sim <- SummarizedExperiment(assays=SimpleList(counts=assay(simulation$sum_exp_raw_count)), colData=colData(simulation$sum_exp_raw_count), rowData=rowData)

  ora <- sbea(de_sim, gs = simulation$raw.gs, perm = 0, method = "ora", padj.method = "BH")
  gsa$ora[[i]] <- ora$res.tbl
  
  gsea <- sbea(de_sim, gs = simulation$raw.gs, perm = 1000, method = "gsea", padj.method = "BH")
  gsa$gsea[[i]] <- gsea$res.tbl

  gsea <- sbea_abs(de_sim, gs = simulation$raw.gs, perm = 1000, method = "gsea", padj.method = "BH")
  gsa$GSEA_absolute[[i]] <- gsea$res.tbl

  res <- results(dds, lfcThreshold = log2(1.5))
  rowData <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res$log2FoldChange, "ADJ.PVAL" = res$padj)
  de_sim <- SummarizedExperiment(assays=SimpleList(counts=assay(simulation$sum_exp_raw_count)), colData=colData(simulation$sum_exp_raw_count), rowData=rowData)
  ora <- sbea(de_sim, gs = simulation$raw.gs, perm = 0, method = "ora", padj.method = "BH")
  gsa$ora_res_1.5[[i]] <- ora$res.tbl

  res <- results(dds, lfcThreshold = log2(1.3))
  rowData <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res$log2FoldChange, "ADJ.PVAL" = res$padj)
  de_sim <- SummarizedExperiment(assays=SimpleList(counts=assay(simulation$sum_exp_raw_count)), colData=colData(simulation$sum_exp_raw_count), rowData=rowData)
  ora <- sbea(de_sim, gs = simulation$raw.gs, perm = 0, method = "ora", padj.method = "BH")
  gsa$ora_res_1.3[[i]] <- ora$res.tbl

  res <- results(dds, lfcThreshold = log2(1.1))
  rowData <- data.frame(row.names = paste0("gene", 1:10000), "FC" = res$log2FoldChange, "ADJ.PVAL" = res$padj)
  de_sim <- SummarizedExperiment(assays=SimpleList(counts=assay(simulation$sum_exp_raw_count)), colData=colData(simulation$sum_exp_raw_count), rowData=rowData)
  ora <- sbea(de_sim, gs = simulation$raw.gs, perm = 0, method = "ora", padj.method = "BH")
  gsa$ora_res_1.1[[i]] <- ora$res.tbl
  rm(.Random.seed, envir=globalenv())
  print(i)
}

Sys.time() - start


saveRDS(gsa, file = paste0("Normal_dist_sims_N_",n,"_inset_",inset,"_outset_",outset,"_",array,".RDS"))


