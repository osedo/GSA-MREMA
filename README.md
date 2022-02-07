## Mixture of Random-Effect Meta-Analysis Components
&nbsp;

A novel approach to GSA based on random effect meta-analysis, which compares the distribution of the effect of a sample group on genes in a gene set with the distribution found in genes outside the gene set.

&nbsp;

### Input 

Read in packages, functions and data
```R
source("packages.R")
source("GMM_7DF.R")
load("COAD_SummarizedExperiment.RDs")
```

A SummarizedExperiment object is required. Details on creating a SummarizedExperiment object can be found [here](https://www.bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html), in part 3. The data used in our analysis was downloaded from the [GDC portal](https://portal.gdc.cancer.gov/repository) with the GDC client. There is a useful tutorial on how to do this and create gene expression matrices [here](https://www.biostars.org/p/9500223/). The code used for creating the COAD_SummarizedExperiment object data is shown in ```tcga_data_wrangling.R```.   

```R
COAD_SummarizedExperiment
```
```R
class: SummarizedExperiment 
dim: 5230 501 
metadata(0):
assays(1): counts
rownames(5230): ENSG00000000419 ENSG00000000938 ... ENSG00000280236 ENSG00000280314
rowData names(0):
colnames(501): TCGA-D5-6530-01A TCGA-G4-6320-01A ... TCGA-A6-5657-01A TCGA-AA-3688-01A
colData names(5): GROUP Patient Age Sex menopause
```

An R list of gene sets with each element corresponding to a gene set. 

```R
gene_sets[c(1,2)]
```

```R
$KEGG_ABC_TRANSPORTERS
 [1] "ENSG00000165029" "ENSG00000154263" "ENSG00000144452" "ENSG00000179869" "ENSG00000107331" "ENSG00000167972"
 [7] "ENSG00000198691" "ENSG00000154265" "ENSG00000154262" "ENSG00000064687" "ENSG00000141338" "ENSG00000154258"
[13] "ENSG00000085563" "ENSG00000135776" "ENSG00000073734" "ENSG00000005471" "ENSG00000004846" "ENSG00000115657"
[19] "ENSG00000131269" "ENSG00000197150" "ENSG00000150967" "ENSG00000103222" "ENSG00000124574" "ENSG00000121270"
[25] "ENSG00000140798" "ENSG00000023839" "ENSG00000108846" "ENSG00000125257" "ENSG00000114770" "ENSG00000091262"
[31] "ENSG00000006071" "ENSG00000069431" "ENSG00000101986" "ENSG00000173208" "ENSG00000117528" "ENSG00000119688"
[37] "ENSG00000160179" "ENSG00000118777" "ENSG00000172350" "ENSG00000138075" "ENSG00000143921" "ENSG00000001626"
[43] "ENSG00000168394" "ENSG00000204267"

$KEGG_ACUTE_MYELOID_LEUKEMIA
 [1] "ENSG00000142208" "ENSG00000105221" "ENSG00000117020" "ENSG00000078061" "ENSG00000002330" "ENSG00000157764"
 [7] "ENSG00000133101" "ENSG00000110092" "ENSG00000245848" "ENSG00000213341" "ENSG00000187840" "ENSG00000122025"
[13] "ENSG00000177885" "ENSG00000174775" "ENSG00000104365" "ENSG00000269335" "ENSG00000173801" "ENSG00000157404"
[19] "ENSG00000133703" "ENSG00000138795" "ENSG00000169032" "ENSG00000126934" "ENSG00000100030" "ENSG00000102882"
[25] "ENSG00000198793" "ENSG00000136997" "ENSG00000109320" "ENSG00000213281" "ENSG00000121879" "ENSG00000051382"
[31] "ENSG00000171608" "ENSG00000105851" "ENSG00000145675" "ENSG00000105647" "ENSG00000117461" "ENSG00000141506"
[37] "ENSG00000137193" "ENSG00000102096" "ENSG00000140464" "ENSG00000112033" "ENSG00000132155" "ENSG00000131759"
[43] "ENSG00000173039" "ENSG00000108443" "ENSG00000175634" "ENSG00000159216" "ENSG00000079102" "ENSG00000115904"
[49] "ENSG00000100485" "ENSG00000066336" "ENSG00000168610" "ENSG00000126561" "ENSG00000173757" "ENSG00000081059"
[55] "ENSG00000152284" "ENSG00000148737" "ENSG00000109906" 
```




&nbsp;



### GSA

Looking for significantly enriched gene sets between tumour and normal tissue samples.

```R
COAD_enrich <- GMM_1DF(COAD_SummarizedExperiment, gene_sets)
head(COAD_enrich)
```

```R
GeneSets                                        `p-values` mid_dist  size BIC_value `Adj-Pval`
1 KEGG_ABC_TRANSPORTERS                                1     TRUE        43 FALSE              1
2 KEGG_ACUTE_MYELOID_LEUKEMIA                          0.165 FALSE       57 FALSE              1
3 KEGG_ADHERENS_JUNCTION                               0.241 FALSE       72 FALSE              1
4 KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY                 0.878 TRUE        63 FALSE              1
5 KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM      0.271 FALSE       30 FALSE              1
6 KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION       0.263 FALSE       37 FALSE              1
```

The p-value and adjusted p-values are shown. The criteria for the middle component having a smaller weight in the inset than the outset is shown in the mid_dist column.

&nbsp;

Looking for gene sets with a significantly different underlying LFC distribution between the inset and the outset.

```R
COAD_difference <- GMM_7DF(COAD_SummarizedExperiment, gene_sets)
head(COAD_difference)
```
```R
  GeneSets                                        `p-values` mid_dist  size BIC_value `Adj-Pval`
1 KEGG_ABC_TRANSPORTERS                                0.214 TRUE        43 FALSE          0.924
2 KEGG_ACUTE_MYELOID_LEUKEMIA                          0.826 TRUE        57 FALSE          1    
3 KEGG_ADHERENS_JUNCTION                               0.151 TRUE        72 FALSE          0.787
4 KEGG_ADIPOCYTOKINE_SIGNALING_PATHWAY                 0.940 TRUE        63 FALSE          1    
5 KEGG_ALANINE_ASPARTATE_AND_GLUTAMATE_METABOLISM      0.552 TRUE        30 FALSE          1    
6 KEGG_ALDOSTERONE_REGULATED_SODIUM_REABSORPTION       0.937 TRUE        37 FALSE          1    
```

&nbsp;

### More comparisons

In order include covariates and interaction terms and to carry out more specific comparisons alter the DESeq2 design in the second block of code in the ```GMM_1DF()``` and ```GMM_7DF()``` functions. 

&nbsp;
Example - adding a covariate and interaction term for whether an individual in pre-menopausal or post-menopausal/male.
```R
### use Deseq2 to get group effect estimate and it's standard error
dea_raw_count<- DESeqDataSetFromMatrix(countData = se, colData = colData(sum_exp_raw_count), design = ~ GROUP + menopause + GROUP:menopause)
### this line sets pre-menopausal females as the reference group
dea_raw_count$menopause <- relevel(dea_raw_count$menopause, ref = "pre")
dds <- DESeq(dea_raw_count)
### this gives the effect of cancer on pre-menopausal female samples
res <- results(dds, name = "GROUP_1_vs_0")
### uncomment this line for the effect of cancer on samples that are not pre-menopausal females
#res <- results(dds, list(c("GROUP_1_vs_0", "GROUP1.menopausepost.none"))
## uncomment this line for the difference in cancer effect between our two subsets of samples
#res <- results(dds, name = "GROUP1.menopausepost.none")
```

There are more details on DESeq2 expirmental design [here](https://rstudio-pubs-static.s3.amazonaws.com/329027_593046fb6d7a427da6b2c538caf601e1.html)
