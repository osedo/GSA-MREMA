# GSA-MREMA

### Files
1. **packages.R**- Loads the required packages.
2. **functions.R** - Defines functions used for the EM algorithm.
3. **globalfunctions.R** - Functions for Gene Set Anlysis and csGSA.
4. **sumexp_sets.R** - Downloads and orders count data from TCGA and Gene sets from MSigDB. Returns gene set and summarizedExperiment object.
5. **Combined_Purity_Estimates.txt** - Purity Estimates - Data not all available for cancer types.
6. **Luminal_Sum_Exp_KEGG.Rdata** - SummarizedExperiment object - Luminal A and Luminal B breast cancer samples.
7. **BR_KEGG.Rdata** - Gene set list - KEGG subset of canonical pathways.

&nbsp;

&nbsp;


### Load Packages & Functions 
```R 
source("packages.R")
source("functions.R")
source("simulateddata.R")
source("globalfunctions.R")
```


&nbsp;
&nbsp;

### Data Sets

```R
load("BR_KEGG.Rdata") # gene_sets
load("Luminal_Sum_Exp_KEGG.Rdata") # sum_exp_raw_count
```

The data consists of a SummarizedExperiment object which has the raw count data for TCGA breast cancer samples. Luminal A tumours are classed as group 0 and luminal B tumours are classed as group 1. The purity given is the proportion of normal cells in a sample.  
The gene sets given are from the KEGG subset of the canonical pathways, taken from the MSigDB Collections.  
Only genes that overlapped between sets were kept for analysis.  

```R
sum_exp_raw_count
```



&nbsp;

### Normalize Data
We follow DESeq2 workflow to normalize data, this is optional users can use other packages. It should be noted that users could also use microarray data.

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

### RUN GSA

```R
w_o_mrema(sum_exp_normalized_object, gene_sets)
```

`mrema()` and `wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively. 

The p-value is given for the likelihood-ratio test. This should be corrected for multiple testing and is conditional on the mid_dist being TRUE. mid_dist tests whether the weight assigned to the non-DE genes in a gene set is smaller in the gene set than in the genes outside the gene set.
The size of the gene set is also given.

Finally it is tested whether the Bayesian information criterion (BIC) of the gene is lower than the BIC for all.


&nbsp;

&nbsp;

### Run csGSA (cell type specific GSA)

Gene sets enriched in luminal B cancer cells when compared to luminal A cancer cells.

```R
cancer_w_o_mrema(sum_exp_normalized_object, gene_sets)
```
`cancer_mrema()` and `cancer_wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively.
&nbsp;

Genes sets found to be enriched in normal cells in luminal B tumour when compared to normal cells in luminal A tumours.

```R
normal_w_o_mrema(simulation_cancer$sum_exp_raw_count, simulation_cancer$raw.gs)
```

`normal_mrema()` and `normal_wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively.




