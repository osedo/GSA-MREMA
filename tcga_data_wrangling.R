
library(msigdbr)
library(tidyverse)
library(data.table)
library(SummarizedExperiment)


generate_count_mat <- function(path, pattern) {
  files = list.files(path, pattern, full.names = TRUE, recursive=TRUE, include.dirs=TRUE)
  mat = as.data.frame(do.call(cbind, lapply(files, function(x) fread(x, stringsAsFactors = FALSE))))
  rownames(mat) = mat[,1]
  mat = as.data.frame(mat[, seq(2, ncol(mat), 2)])
  return(mat)
}


#create gene expression matrix
counts <- generate_count_mat("TCGA-COAD/", "\\.htseq.counts$")


# get sample information
samplesheet <- read.table("gdc_sample_sheet_colon.tsv", header=T, sep="\t")
colnames(counts) <- samplesheet$Sample.ID
counts <- counts[,!duplicated(colnames(counts))]
meta <- subset(samplesheet, select=c(Sample.ID, Sample.Type))
meta <- meta[!duplicated(meta$Sample.ID),]
table(meta$Sample.Type)

### remove all samples that are not primary tumour or solid tissue normal
m <- which(meta$Sample.Type == "Metastatic")
counts <- counts[,-m]
meta <- meta[-m,]
r <- which(meta$Sample.Type == "Recurrent Tumor")
counts <- counts[,-r]
meta <- meta[-r,]

## left with 503 samples 
table(meta$Sample.ID == colnames(counts))


### create group data for analysis
group <- rep(0, length(meta$Sample.Type))
group[which(meta$Sample.Type == "Solid Tissue Normal")] <- 0
group[which(meta$Sample.Type == "Primary Tumor")] <- 1


### download kegg gene sets....we only keep genes that are in our sets and expression matrix
canon_gene_sets = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
row.names(counts) <- gsub("\\..*","",row.names(counts))
canon_gene_sets_overlap <- canon_gene_sets[which(canon_gene_sets$ensembl_gene %in% row.names(counts) == TRUE),]
counts_overlap <- counts[which(row.names(counts) %in% canon_gene_sets_overlap$ensembl_gene == TRUE),]


## create summarizedexperiment object
rowData<- DataFrame(row.names = rownames(counts_overlap))
colData <- DataFrame(GROUP=group, row.names= colnames(counts_overlap))
COAD_SummarizedExperiment <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts_overlap)), colData=colData, rowData=rowData)
COAD_SummarizedExperiment


## create gene set list
gene_sets <- list()
set_names <- unique(canon_gene_sets_overlap$gs_name)
for(i in 1:length(set_names)){
  inset <- canon_gene_sets_overlap[which(canon_gene_sets_overlap$gs_name == set_names[i]),6]
  inset <- as.character(inset$ensembl_gene)
  gene_sets[[i]] <- c(inset)
}
names(gene_sets) <- set_names


## extracting data from the clinical.tsv file and adding to the the coldata
clinical <- read_tsv("clinical.tsv")
case_expr_id <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1", colnames(assay(COAD_SummarizedExperiment)))
## match based on case submitter id, here we get age and gender
new_coldata <- clinical[match(case_expr_id, clinical$case_submitter_id),c(2,4,12)]
rownames(new_coldata) <- colnames(assay(COAD_SummarizedExperiment))
colData_full <- cbind(colData(COAD_SummarizedExperiment), new_coldata)
colnames(colData_full) <- c("GROUP", "Patient", "Age", "Sex")
colData_full$GROUP <- as.factor(colData_full$GROUP)
colData_full$Sex <- as.factor(colData_full$Sex)

## categrize samples based on age and gender
meno <- c()
meno[which(colData_full$Age <= 44 & colData_full$Sex == "female")] <- "pre"
meno[which(colData_full$Age > 44 | colData_full$Sex == "male")] <- "post.none"
colData_full$menopause <- as.factor(meno)
assay <- assay(COAD_SummarizedExperiment)[,complete.cases(colData_full)]
colData_full <- colData_full[complete.cases(colData_full),]
rowData<- DataFrame(row.names = rownames(assay))
COAD_SummarizedExperiment <- SummarizedExperiment(assays=SimpleList(counts=assay), colData=colData_full, rowData=rowData)
COAD_SummarizedExperiment


