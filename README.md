## Gene set analysis & cell-type specific gene set analysis
&nbsp;
### 1. Required
1. A SummarizedExperiment object with count data and group and purity information for each sample.
2. A list of gene sets.
3. The functions and packages found in the scripts above. 



### 2. Load Packages, Functions & Data
```R 
source("packages.R")
source("functions.R")
source("simulateddata.R")
source("globalfunctions.R")
```

&nbsp;

```R
load("BR_KEGG.Rdata") # gene_sets
load("Luminal_Sum_Exp_KEGG.Rdata") # sum_exp_raw_count
```

The data consists of a SummarizedExperiment object which has the raw count data for TCGA breast cancer samples. Luminal A tumours are classed as group 0 and luminal B tumours are classed as group 1. The purity given is the proportion of normal cells in a sample.  
The gene sets given are from the KEGG subset of the canonical pathways, taken from the MSigDB Collections.  
Only genes that were in these gene sets were kept in the count data. 

```R
sum_exp_raw_count
```
```R
class: SummarizedExperiment 
dim: 5102 760 
metadata(0):
assays(1): counts
rownames(5102): DPM1 FGR ... OR12D2 OR8K3
rowData names(1): rownames.combined_counts.
colnames(760): TCGA-E9-A245 TCGA-AC-A6NO ... TCGA-BH-A0BD TCGA-AN-A0AM
colData names(2): GROUP purity
```
Note that purity corresponds to the proportion of normal cells in a sample.



&nbsp;

### 3. Normalize Data (Optional)
We follow DESeq2 workflow to normalize data, this is optional and users can use other packages/workflows. It should be noted that users could also use microarray data.

```R
dea_raw_count<- DESeqDataSetFromMatrix(countData = assay(sum_exp_raw_count), colData = colData(sum_exp_raw_count), design = ~ GROUP)
# now we estimate the size factors
dea_raw_count<- estimateSizeFactors(dea_raw_count)
# Extract the normalized counts, this matrix must then go to the gsa
sum_exp_normalized_counts<- counts(dea_raw_count, normalized=TRUE)
sum_exp_normalized_counts<- log2(sum_exp_normalized_counts+0.5)
# create new data frame with normalized counts.
sum_exp_normalized_object <- SummarizedExperiment(assays=SimpleList(counts=sum_exp_normalized_counts), colData=colData(sum_exp_raw_count), rowData=DataFrame(rownames(sum_exp_normalized_counts)))
```
We have created a new SummarizedExperiment object with normalized counts instead of raw counts.

### 4. a) GSA

A standard GSA can be run to find gene sets enriched for differentially expressed genes between luminal A and luminal B subtypes of breast cancer. 
```R
w_o_mrema(sum_exp_normalized_object, gene_sets)
```

`mrema()` and `wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively. 

The p-value is given for the likelihood-ratio test. This should be corrected for multiple testing and is conditional on the mid_dist being TRUE. mid_dist tests whether the weight assigned to the non-DE genes in a gene set is smaller in the gene set than in the genes outside the gene set.
The size of the gene set is also given.

Finally it is tested whether the Bayesian information criterion (BIC) of the gene is lower than the BIC for all.



### 4. b) csGSA (cell type specific GSA)

Gene sets enriched differentially expressed genes in luminal B cancer cells when compared to luminal A cancer cells.

```R
cancer_w_o_mrema(sum_exp_normalized_object, gene_sets)
```
`cancer_mrema()` and `cancer_wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively.
&nbsp;

Genes sets found to be enriched in normal cells in luminal B tumours when compared to normal cells in luminal A tumours.

```R
normal_w_o_mrema(simulation_cancer$sum_exp_raw_count, simulation_cancer$raw.gs)
```

`normal_mrema()` and `normal_wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively.




