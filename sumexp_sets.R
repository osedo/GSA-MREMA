

library(TCGAbiolinks)
library(msigdbr)
library(biomaRt)

# download tumour subtype data 
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/tcga-study-abbreviations

cancer_type = "BRCA"
tumour_subtypes <- TCGAquery_subtype(tumor = cancer_type)
# pick rownumber with the categories/subtypes, changes depeneding on tumour, in BRCA it is 12
table(tumour_subtypes[,12])
row_number <- 12

## project name to download tumour expression data
project_name <- "TCGA-BRCA"

## the subtypes you want to compare - take from the table above
group_0 <- "LumA"
group_1 <- "LumB"

## genes you want to analyse - pick a category and up to six subcategories
## see collections here http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
gs_category <- "C2"
gs_sub_category = c("CP:KEGG")

### purity data - careful not available for all tumour type
### maybe an updated file around?
purity <- read.table("Combined_Purity_Estimates.txt", header = TRUE)
table(purity$Cancer_type)


## Everthing below this point could go in a function? Returning summarized experiment and a gene set


## download tumour gene expression counts 
query <- GDCquery(project = project_name,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = "Primary Tumor")
GDCdownload(query)
data <- GDCprepare(query)

# pull out Group 0 and Group 1 samples
Group0_samples <- tumour_subtypes[which(tumour_subtypes[,row_number] == group_0 ),1]
Group1_samples <- tumour_subtypes[which( tumour_subtypes[,row_number] == group_1 ),1]

# get same name format as subgroup data
colnames(data) <- substr(colnames(data), start = 1, stop = 12)


# create two summarized experiment object luminal a and luminal b
group_0_SEO <- subset(data, select = colnames(data) %in% Group0_samples$patient)
group_1_SEO <- subset(data, select = colnames(data) %in% Group1_samples$patient)

# remove any duplicated samples
dup_0 <- group_0_SEO$patient[which(duplicated(group_0_SEO$patient)== TRUE)]
dup_1 <- group_1_SEO$patient[which(duplicated(group_1_SEO$patient)== TRUE)]
group_0_SEO <- subset(group_0_SEO, select = !(colnames(group_0_SEO) %in% dup_0))
group_1_SEO <- subset(group_1_SEO, select = !(colnames(group_1_SEO) %in% dup_1))


#read in purity data
purity <- purity[which(purity$Cancer_type == cancer_type),c(1,3)]
purity$Sample_ID <- substr(purity$Sample_ID, start = 1, stop = 12)
purity <- purity[complete.cases(purity),]
purity <- purity[!duplicated(purity$Sample_ID),]


# remove samples without purity data
p_data_1 <- group_1_SEO$patient[which(group_1_SEO$patient %in% purity$Sample_ID == TRUE)]
group_1_SEO <- subset(group_1_SEO, select = colnames(group_1_SEO) %in% p_data_1)

p_data_0 <- group_0_SEO$patient[which(group_0_SEO$patient %in% purity$Sample_ID == TRUE)]
group_0_SEO <- subset(group_0_SEO, select = colnames(group_0_SEO) %in% p_data_0)



# order purity by samples in summarizedexperiment object
purity_1 <- purity %>% arrange(factor(Sample_ID, levels = group_1_SEO$patient))
purity_1 <- purity_1[c(1:length(group_1_SEO$patient)),]

purity_0 <- purity %>% arrange(factor(Sample_ID, levels = group_0_SEO$patient))
purity_0 <- purity_0[c(1:length(group_0_SEO$patient)),]



# create group data for summarizedExperiment object
group_0 <- rep(0,length(group_0_SEO$patient))
group_1 <- rep(1,length(group_1_SEO$patient))
group <- c(group_0, group_1)
# create purity object
purity_vector <- c(purity_0$ESTIMATE, purity_1$ESTIMATE)


# make a combined expression matrix
group0_counts <- as.data.frame(assay(group_0_SEO, "HTSeq - Counts"))
group1_counts <- as.data.frame(assay(group_1_SEO, "HTSeq - Counts"))
combined_counts <- cbind(group0_counts, group1_counts)


# download gene sets of interest
# http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
msigdb_gene_sets <- data_frame()

if (length(gs_sub_category) == 0){
  
  msigdb_gene_sets = msigdbr(species = "Homo sapiens", category = gs_category)
} else if (length(gs_sub_category) == 1){
  msigdb_gene_sets = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category)
} else if (length(gs_sub_category) == 2){
  msigdb_gene_sets_1 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[1])
  msigdb_gene_sets_2 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[2])
  msigdb_gene_sets <- rbind(msigdb_gene_sets_1, msigdb_gene_sets_2)
} else if (length(gs_sub_category) == 3){
  msigdb_gene_sets_1 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[1])
  msigdb_gene_sets_2 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[2])
  msigdb_gene_sets_3 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[3])
  msigdb_gene_sets <- rbind(msigdb_gene_sets_1, msigdb_gene_sets_2)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_3)
} else if (length(gs_sub_category) == 4){
  msigdb_gene_sets_1 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[1])
  msigdb_gene_sets_2 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[2])
  msigdb_gene_sets_3 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[3])
  msigdb_gene_sets_4 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[4])
  msigdb_gene_sets <- rbind(msigdb_gene_sets_1, msigdb_gene_sets_2)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_3)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_4)
} else if (length(gs_sub_category) == 5){
  msigdb_gene_sets_1 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[1])
  msigdb_gene_sets_2 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[2])
  msigdb_gene_sets_3 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[3])
  msigdb_gene_sets_4 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[4])
  msigdb_gene_sets_5 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[5])
  msigdb_gene_sets <- rbind(msigdb_gene_sets_1, msigdb_gene_sets_2)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_3)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_4)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_5)
} else if (length(gs_sub_category) == 6){
  msigdb_gene_sets_1 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[1])
  msigdb_gene_sets_2 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[2])
  msigdb_gene_sets_3 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[3])
  msigdb_gene_sets_4 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[4])
  msigdb_gene_sets_5 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[5])
  msigdb_gene_sets_6 = msigdbr(species = "Homo sapiens", category = gs_category, subcategory = gs_sub_category[6])
  msigdb_gene_sets <- rbind(msigdb_gene_sets_1, msigdb_gene_sets_2)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_3)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_4)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_5)
  msigdb_gene_sets <- rbind(msigdb_gene_sets, msigdb_gene_sets_6)
}


## we have to change our rownames (gene names) in the expression matrix to a format that matches the gene sets we downloaded
ensembl.genes <- rownames(combined_counts)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c( "ensembl_gene_id", "hgnc_symbol"),
  values=ensembl.genes,
  mart=mart)


## only keep genes that are in our gene sets
genes_insets <- genes[which(genes$hgnc_symbol %in% msigdb_gene_sets$human_gene_symbol == TRUE),]

## remove genes whose names/ids are duplicated
dup_ens <- which(duplicated(genes_insets$ensembl_gene_id) == TRUE)

if (length(dup_ens) > 0){
  genes_insets <- genes_insets[-dup_ens,]
}

dup_sym <- which(duplicated(genes_insets$hgnc_symbol) == TRUE)
if (length(dup_sym) > 0){
  genes_insets <- genes_insets[-dup_sym,]
}


### only keep the count data for genes that are insets/have names from biomart
combined_counts <- combined_counts[which(rownames(combined_counts) %in% genes_insets$ensembl_gene_id == TRUE),]
## check order is correct
which((rownames(combined_counts) == genes_insets$ensembl_gene_id) == FALSE)
## change rownames to corresponing symbol
rownames(combined_counts) <- genes_insets$hgnc_symbol


## create summarizedexperiment object
colData <- DataFrame(GROUP=group, row.names=colnames(combined_counts), purity = (1 - purity_vector))
rowData <- DataFrame(rownames(combined_counts))
rownames(combined_counts) <- NULL
colnames(combined_counts) <- NULL
sum_exp_raw_count<- SummarizedExperiment(assays=SimpleList(counts=as.matrix(combined_counts)), colData=colData, rowData=rowData)
rownames(sum_exp_raw_count) <- genes_insets$hgnc_symbol



gene_sets <- list()
set_names <- unique(msigdb_gene_sets$gs_name)
for(i in 1:length(set_names)){
  inset <- msigdb_gene_sets[which(msigdb_gene_sets$gs_name == set_names[i]),7]
  inset <- as.character(inset$human_gene_symbol)
  gene_sets[[i]] <- c(inset)
}
names(gene_sets) <- set_names


save(gene_sets, file = "BR_KEGG.Rdata")
save(sum_exp_raw_count, file = "Luminal_Sum_Exp_KEGG.Rdata")
