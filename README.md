## Gene set analysis & cell-type specific gene set analysis
&nbsp;

A novel approach to GSA based on random effect meta-analysis, which compares the distribution of the effect of a sample group on genes in a gene set with the distribution found in genes outside the gene set. Extending this approach to the cell-type specific GSA problem allows the identification of gene sets for which the group effect differs between cell types in the samples. 

&nbsp;

### Input 

A SummarizedExperiment object is required. Details on creating a SummarizedExperiment object can be found [here](https://www.bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html), in part 3. If you wish to continue with the cell-type specific problem the estimates of the cell propotion of interst should be added to the colData of the SummarizedExperiment object. In this example the proportion of cancer cells in a sample is supplied in the colData as purity.

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

An R list of gene sets with each element corresponding to a gene set. The name of the element should correspond to the gene set while the element itself contains the constituent genes. 

```R
gene_sets[1]
```

```R
$KEGG_ABC_TRANSPORTERS
 [1] "ABCA1"  "ABCA10" "ABCA12" "ABCA13" "ABCA2"  "ABCA3"  "ABCA4"  "ABCA5"  "ABCA6"  "ABCA7"  "ABCA8"  "ABCA9"  "ABCB1"  "ABCB10" "ABCB11" "ABCB4"  "ABCB5"  "ABCB6"  "ABCB7"  "ABCB8"  "ABCB9" 
[22] "ABCC1"  "ABCC10" "ABCC11" "ABCC12" "ABCC2"  "ABCC3"  "ABCC4"  "ABCC5"  "ABCC6"  "ABCC8"  "ABCC9"  "ABCD1"  "ABCD2"  "ABCD3"  "ABCD4"  "ABCG1"  "ABCG2"  "ABCG4"  "ABCG5"  "ABCG8"  "CFTR"  
[43] "TAP1"   "TAP2"  
```




&nbsp;

### Preprocessing
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

### GSA

A standard GSA can be run to find gene sets enriched for differentially expressed genes between luminal A and luminal B subtypes of breast cancer. 
```R
w_o_mrema(sum_exp_normalized_object, gene_sets)
```

`mrema()` and `wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively. 

The p-value is given for the likelihood-ratio test. This should be corrected for multiple testing and is conditional on the mid_dist being TRUE. mid_dist tests whether the weight assigned to the non-DE genes in a gene set is smaller in the gene set than in the genes outside the gene set.
The size of the gene set is also given.

Finally it is tested whether the Bayesian information criterion (BIC) of the gene is lower than the BIC for all.



### csGSA (cell type specific GSA)

Gene sets enriched differentially expressed genes in luminal B cancer cells when compared to luminal A cancer cells.

```R
cell_type_w_o_mrema(sum_exp_normalized_object, gene_sets)
```
`cell_type_mrema()` and `cell_type_wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively.
&nbsp;

Genes sets found to be enriched in normal cells in luminal B tumours when compared to normal cells in luminal A tumours can be found by changing the purity column values to (1 - purity) and running the cell-type specific functions.




