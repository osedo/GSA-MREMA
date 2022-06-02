## Mixture of Random-Effect Meta-Analysis
&nbsp;

A novel approach to GSA based on random effect meta-analysis, which compares the distribution of the effect of a sample group on genes in a gene set with the distribution found in genes outside the gene set.

&nbsp;

### Running MREMA

Looking for significantly enriched gene sets between tumour and normal tissue samples.
```R
source("packages.R")
source("mrema.R")
# using the DF parameter to call 1DF test
womrema_1.5 <- mrema(postdata_overlap, kegg.gs.overalp, DF = 1, threshold = 1.5, CI_interval = FALSE)
head(womrema_1.5)
```
```R
  GeneSets                                          `p-values` mid_dist  size BIC_value `Adj-Pval`
  <chr>                                                  <dbl> <lgl>    <dbl> <lgl>          <dbl>
1 hsa00010_Glycolysis_/_Gluconeogenesis                 0.789  FALSE       53 FALSE         0.883 
2 hsa00020_Citrate_cycle_(TCA_cycle)                    0.479  FALSE       29 FALSE         0.676 
3 hsa00030_Pentose_phosphate_pathway                    0.882  FALSE       25 FALSE         0.931 
4 hsa00040_Pentose_and_glucuronate_interconversions     0.0128 TRUE        23 FALSE         0.0471
5 hsa00051_Fructose_and_mannose_metabolism              0.417  FALSE       32 FALSE         0.624 
6 hsa00052_Galactose_metabolism                         0.251  TRUE        27 FALSE         0.433 
```
The p-value and adjusted p-values are shown. The criteria for the middle component having a smaller weight in the inset than the outset is shown in the mid_dist column.

&nbsp;

Looking for gene sets with a significantly different underlying LFC distribution between the inset and the outset, this allows for a parameters in the mixture of each group of genes to be independent across gene groups.

```R
mrema_1.5 <- mrema(postdata_overlap, kegg.gs.overalp, DF = 6, threshold = 1.5, CI_interval = FALSE)
head(mrema_1.5)
```
```R
  GeneSets                                          `p-values` mid_dist  size BIC_value `Adj-Pval`
1 hsa00010_Glycolysis_/_Gluconeogenesis              0.000297  TRUE        53 FALSE      0.000748 
2 hsa00020_Citrate_cycle_(TCA_cycle)                 0.0465    TRUE        29 FALSE      0.0656   
3 hsa00030_Pentose_phosphate_pathway                 0.0000299 TRUE        25 FALSE      0.0000960
4 hsa00040_Pentose_and_glucuronate_interconversions  0.0000910 TRUE        23 FALSE      0.000271 
5 hsa00051_Fructose_and_mannose_metabolism           0.0131    TRUE        32 FALSE      0.0217   
6 hsa00052_Galactose_metabolism                      0.198     TRUE        27 FALSE      0.237  
```

### Inputs
 Require a list of gene sets and a dataframe with the gene names, logfold change estimate and standard error of the estimate in three columns. Any DE method can be used and effects combined with interation terms to make comparison between subsets.

```R
head(postdata_overlap)
```
```
  Ensembl          effect variance
   <chr>             <dbl>    <dbl>
 1 ENSG00000000419  0.620   0.0151 
 2 ENSG00000000938  0.0837  0.0293 
 3 ENSG00000000971 -1.40    0.0255 
 4 ENSG00000001036  0.136   0.00663
 5 ENSG00000001084  0.0182  0.00636
 6 ENSG00000001167  0.486   0.00540
 7 ENSG00000001561 -0.617   0.0149 
 8 ENSG00000001617  0.986   0.0317 
 9 ENSG00000001626 -0.499   0.0687 
10 ENSG00000001630  0.976   0.0328 
```

```R
kegg.gs.overalp[c(1,2)]
```

```R
$`hsa00010_Glycolysis_/_Gluconeogenesis`
 [1] "ENSG00000006534" "ENSG00000067057" "ENSG00000067225" "ENSG00000072210" "ENSG00000074800" "ENSG00000079739"
 [7] "ENSG00000091140" "ENSG00000100889" "ENSG00000102144" "ENSG00000105220" "ENSG00000106633" "ENSG00000107789"
[13] "ENSG00000108515" "ENSG00000108602" "ENSG00000109107" "ENSG00000111275" "ENSG00000111640" "ENSG00000111669"
[19] "ENSG00000111674" "ENSG00000111716" "ENSG00000117448" "ENSG00000124253" "ENSG00000131069" "ENSG00000131828"
[25] "ENSG00000132746" "ENSG00000134333" "ENSG00000136872" "ENSG00000137124" "ENSG00000141349" "ENSG00000141959"
[31] "ENSG00000143149" "ENSG00000143627" "ENSG00000143891" "ENSG00000149925" "ENSG00000150768" "ENSG00000152556"
[37] "ENSG00000154930" "ENSG00000156510" "ENSG00000156515" "ENSG00000159322" "ENSG00000159399" "ENSG00000160883"
[43] "ENSG00000164904" "ENSG00000165140" "ENSG00000168291" "ENSG00000169299" "ENSG00000171314" "ENSG00000172331"
[49] "ENSG00000172955" "ENSG00000196616" "ENSG00000197894" "ENSG00000198099" "ENSG00000248144"

$`hsa00020_Citrate_cycle_(TCA_cycle)`
 [1] "ENSG00000014641" "ENSG00000062485" "ENSG00000067829" "ENSG00000073578" "ENSG00000091140" "ENSG00000091483"
 [7] "ENSG00000100412" "ENSG00000100889" "ENSG00000101365" "ENSG00000105953" "ENSG00000117118" "ENSG00000119689"
[13] "ENSG00000122729" "ENSG00000124253" "ENSG00000131473" "ENSG00000131828" "ENSG00000136143" "ENSG00000138413"
[19] "ENSG00000143252" "ENSG00000146701" "ENSG00000150768" "ENSG00000163541" "ENSG00000166411" "ENSG00000168291"
[25] "ENSG00000172340" "ENSG00000173599" "ENSG00000182054" "ENSG00000197444" "ENSG00000204370"
```


&nbsp;


### Running DE analysis 

The input for the mrema() function is a dataframe with gene names, the effect size being examined and the standard error of that effect size as columns. Here we use DESeq2 to estimate the effect of cancer on gene expression. Other tools can be used and other effects can be tested. 

```R
load("~/data/COAD__SummarizedExperiment.RData")
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
postdata <- tibble("Ensembl" = substr(rownames(res),1,15), "effect" = res[,2], "variance" = res[,3]^2)
postdata <- postdata[complete.cases(postdata),]

postdata_overlap <- postdata[which(postdata$Ensembl %in% unlist(kegg.gs)),]
kegg.gs.overalp <- lapply(kegg.gs, function(s) s[s %in% postdata_overlap$Ensembl]) 
```
&nbsp;

