args <- commandArgs(TRUE)

cancer <- args[1]

source("../packages.R")
source("../mrema.R")
load(paste0("~/data/",cancer,"_SummarizedExperiment.RData"))

kegg.gs <- getGenesets(org="hsa", db="kegg", gene.id.type = "ENSEMBL")

## run DESeq2 with individual/block as covariate
se <- assay(SummarizedExperiment)
group <- colData(SummarizedExperiment)$GROUP
keep <- filterByExpr(se, group = group)
se <- se[keep,]
dea_raw_count<- DESeqDataSetFromMatrix(countData = se, colData = colData(SummarizedExperiment), design = ~ BLOCK + GROUP)
dea_raw_count$GROUP <- relevel(as.factor(dea_raw_count$GROUP), ref = "Solid Tissue Normal")
dds <- DESeq(dea_raw_count)
res <- results(dds)
res <- res[complete.cases(res),]
postdata <- tibble("Ensembl" = substr(rownames(res),1,15), "effect" = res[,2], "variance" = res[,3]^2, "padj" = res$padj)
postdata <- postdata[complete.cases(postdata),]

postdata_overlap <- postdata[which(postdata$Ensembl %in% unlist(kegg.gs)),]
kegg.gs.overalp <- lapply(kegg.gs, function(s) s[s %in% postdata_overlap$Ensembl]) 

# run mrema
womrema_1.5 <- mrema(postdata_overlap, kegg.gs.overalp, DF = 1, threshold = 1.5)
sig <- womrema_1.5[which(womrema_1.5$`Adj-Pval` < 0.05 & womrema_1.5$mid_dist == TRUE),]
sig <- sig[order(-sig$Upper_CI_nonDE_inset),]

### enrichmentbrowser other methods
se <- assay(SummarizedExperiment)
rownames(se) <- substr(rownames(se), 1, 15)
group <- colData(SummarizedExperiment)$GROUP
group[which(group == "Solid Tissue Normal")] <- 0
group[which(group == "Primary Tumor")] <- 1
colData_2 <- data_frame("GROUP" = group, "BLOCK" = colData(SummarizedExperiment)$BLOCK) 
sum_enrich <- SummarizedExperiment(assays=SimpleList(counts=se), colData=colData_2)

de <- deAna(sum_enrich, grp = "GROUP", blk = "BLOCK", de.method = "DESeq2", padj.method = "BH")
ora <- sbea(de, kegg.gs, method = "ora", perm = 0, padj.method = "BH")
safe <- sbea(de, kegg.gs.overalp, method = "safe", perm = 1000, padj.method = "BH")

## read in an altered version of sbea() which uses the absolute S2N instead of the S2N as the ranking metric
source("../sbea_gsea_abs.R")
gsea <- sbea_abs(de, kegg.gs, method = "gsea", perm = 1000, padj.method = "BH")

gsea_sig <- gsea$res.tbl[which(gsea$res.tbl$PVAL < 0.05),]
ora_sig <- ora$res.tbl[which(ora$res.tbl$PVAL < 0.05),]
safe_sig <- safe$res.tbl[which(safe$res.tbl$PVAL < 0.05),]



library(dplyr)
sig_sets_comparison <- bind_rows(
  as.data.frame(ora_sig[which(ora_sig$ADJ.PVAL < 0.05),c(1,5)]),
  as.data.frame(gsea_sig[which(gsea_sig$ADJ.PVAL < 0.05),c(1,5)]))

method <- c(rep("ORA", length(which(ora_sig$ADJ.PVAL < 0.05))), rep("GSEA", length(which(gsea_sig$ADJ.PVAL < 0.05))))

sig_sets_comparison <- cbind(sig_sets_comparison, method)

sig_sets_comparison$GENE.SET <- substr(sig_sets_comparison$GENE.SET, 10, 100)
sig_sets_comparison$GENE.SET <- gsub("_", " ", sig_sets_comparison$GENE.SET)
sig_sets_comparison$ADJ.PVAL <- signif(sig_sets_comparison$ADJ.PVAL, 2)
sig_sets_comparison <- sig_sets_comparison


sig$GeneSets <- substr(sig$GeneSets, 10, 100)
sig$GeneSets <- gsub("_", " ", sig$GeneSets)
sig$`Adj-Pval` <- signif(sig$`Adj-Pval`, 2)
sig <- sig[order(-sig$Upper_CI_nonDE_inset),]
sig <- sig[,c(1,4,3,6)]

write.csv(sig_sets_comparison, file = paste0(cancer,"_other_tools.csv"), row.names = FALSE, quote = FALSE)
write.csv(sig, file = paste0(cancer,"_1DF_1.5.csv"), row.names = FALSE, quote = FALSE)