load("kegg.RData")
load("de_tcga.RData")
source("benchmarking_functions.R")
source("../mrema.R")
source("../simulations/sbea_gsea_abs.R")
source("../packages.R")
library(EnrichmentBrowser)
library(GSEABenchmarkeR)


gsa <- list()
gsa$mrema_1DF_1.5 <- vector(mode = "list",  length(de_tcga))
gsa$mrema_2DF_1.5 <- vector(mode = "list",  length(de_tcga))
gsa$mrema_6DF_1.5 <- vector(mode = "list",  length(de_tcga))
gsa$ora <-  vector(mode = "list", length(de_tcga))
gsa$gsea_abs <-  vector(mode = "list", length(de_tcga))


for(i in 1:length(de_tcga)){
  
  gsa$mrema_1DF_1.5[[i]]<- DataFrame(mrema_bench_1DF(de_tcga[[i]],threshold = 1.5, kegg.gs))
  gsa$mrema_2DF_1.5[[i]]<- DataFrame(mrema_bench_2DF(de_tcga[[i]],threshold = 1.5, kegg.gs))
  gsa$mrema_6DF_1.5[[i]] <- DataFrame(mrema_bench_6DF(de_tcga[[i]],threshold = 1.5, kegg.gs))
  
  
  gsa$ora[[i]]<- sbea(de_tcga[[i]], gs = kegg.gs, method = "ora", padj.method = "BH", perm = 0)$res.tbl
  gsa$gsea_abs[[i]]<- sbea_abs(de_tcga[[i]], gs = kegg.gs, method = "gsea", padj.method = "BH", perm = 1000)$res.tbl
  
  
}


names(gsa$ora) <- names(de_tcga)
names(gsa$gsea_abs) <- names(de_tcga)
names(gsa$mrema_1DF_1.5) <- names(de_tcga)
names(gsa$mrema_2DF_1.5) <- names(de_tcga)
names(gsa$mrema_6DF_1.5) <- names(de_tcga)


save(gsa, file = "GSA_tcga.RData")

