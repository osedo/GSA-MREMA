## Mixture of Random-Effect Meta-Analysis
&nbsp;

A novel approach to GSA based on random effect meta-analysis, which compares the distribution of the effect of a sample group on genes in a gene set with the distribution found in genes outside the gene set.

&nbsp;


### Running DE analysis 

The input for the mrema() function is a dataframe with gene names, the effect size being examined and the variance of that effect size as columns. Here we use DESeq2 to estimate the effect of cancer on gene expression. Other tools can be used and other effects can be tested. 

Load scripts and data.
```R
load("real_data/kegg.RData")
load("real_data/de_tcga_subset.RData")
source("mrema.R")
source("packages.R")
```

Run DESeq2 on data
```R
## run DESeq2 with individual/block as covariate
se <- assay(de_tcga$BLCA)
group <- colData(de_tcga$BLCA)$GROUP
keep <- filterByExpr(se, group = group)
se <- se[keep,]

dea_raw_count<- DESeqDataSetFromMatrix(countData = se, colData = colData(de_tcga$BLCA), design = ~ BLOCK + GROUP)
dea_raw_count$GROUP <- relevel(as.factor(dea_raw_count$GROUP), ref = "0")
dds <- DESeq(dea_raw_count)
res <- results(dds)
res <- res[complete.cases(res),]
postdata <- tibble("Ensembl" = substr(rownames(res),1,15), "effect" = res[,2], "variance" = res[,3]^2)
postdata <- postdata[complete.cases(postdata),]

# only keep genes which are in a kegg dataset
postdata_overlap <- postdata[which(postdata$Ensembl %in% unlist(kegg.gs)),]
kegg.gs.overlap <- lapply(kegg.gs, function(s) s[s %in% postdata_overlap$Ensembl])

# remove gene sets that are too small/large
GS.MIN.SIZE <- 5
GS.MAX.SIZE <- 500
lens <- lengths(kegg.gs.overlap)
kegg.gs.overlap <- kegg.gs.overlap[lens >= GS.MIN.SIZE & lens <= GS.MAX.SIZE]
```

### Running MREMA

Looking for significantly enriched gene sets between tumour and normal tissue samples.

As the first input there is the dataframe ``` postdata_overlap ``` with the first column being gene names, the second column the LFC estimate and the third column the variance of the LFC estimate. 

```R
DataFrame(postdata_overlap)
```
```R
DataFrame with 5401 rows and 3 columns
         Ensembl     effect  variance
     <character>  <numeric> <numeric>
1              2  -2.813973 0.0623559
2          53947  -1.086577 0.0796909
3           8086   0.431819 0.0127647
4          65985   0.608141 0.0141467
5          51166   0.162179 0.0680566
...          ...        ...       ...
5397        7643  1.2021809 0.1573218
5398        7644  0.2892590 0.0781725
5399       81931  1.2179119 0.0954465
5400       84133  0.0187835 0.0258423
5401        7791 -1.1828887 0.0578679
```
The second input is a list of gene sets ``` kegg.gs.overlap ```.

```R
kegg.gs.overlap[1:2]
```
```
Show in New Window
$`hsa00010_Glycolysis_/_Gluconeogenesis`
 [1] "10327"  "125"    "126"    "128"    "130589" "1737"   "1738"   "2023"   "2026"   "2027"   "217"    "218"    "219"   
[14] "2203"   "221"    "222"    "223"    "224"    "226"    "230"    "2597"   "2821"   "3098"   "3099"   "3101"   "3939"  
[27] "3945"   "441531" "501"    "5106"   "5160"   "5162"   "5211"   "5213"   "5214"   "5223"   "5224"   "5230"   "5232"  
[40] "5236"   "55276"  "55902"  "669"    "7167"   "83440"  "84532"  "92579"  "9562"  

$`hsa00020_Citrate_cycle_(TCA_cycle)`
 [1] "1431" "1737" "1738" "1743" "2271" "3417" "3418" "3419" "3420" "3421" "4190" "4191" "47"   "48"   "4967" "50"   "5091"
[18] "5106" "5160" "5162" "6389" "6390" "6391" "6392" "8801" "8802" "8803"
```

We can now run MREMA. Here we are looking for gene sets enriched for differentially expressed genes so we set the ``` DF ``` parameter to 1. Other options for the parameter are ``` DF = 2 ``` which tests for gene sets with significantly different proportion of upregulated and downregulated genes, and ``` DF = 6 ``` which tests for gene sets with significantly different LFC distributions. 

The ``` threshold ``` parameter allows the researcher to define what they mean by differential expression by passing the absolute fold-change value. 

```R
mrema_1DF <- mrema(postdata = postdata_overlap, raw.gs = kegg.gs.overlap, DF = 1, threshold = 1.5)
DataFrame(mrema_1DF)
```
```R
DataFrame with 337 rows and 6 columns
                  GENE.SET Prop.DE.Increased Estimated.Difference  NR.GENES        PVAL    ADJ.PVAL
               <character>         <numeric>            <numeric> <numeric>   <numeric>   <numeric>
1   hsa00010_Glycolysis_..                 0           -0.0575534        48  0.59059118   0.7453711
2   hsa00020_Citrate_cyc..                 0           -0.3068993        27  0.00586442   0.0270728
3   hsa00030_Pentose_pho..                 0           -0.2473814        24  0.05117969   0.1390932
4   hsa00040_Pentose_and..                 1            0.2961189        15  0.16317360   0.3197064
5   hsa00051_Fructose_an..                 0           -0.1681622        25  0.27955449   0.4595603
...                    ...               ...                  ...       ...         ...         ...
333 hsa05414_Dilated_car..                 1          0.452302082        68 1.46610e-07 4.94077e-06
334 hsa05415_Diabetic_ca..                 0         -0.130236598       156 7.65592e-03 3.20014e-02
335 hsa05416_Viral_myoca..                 1          0.100874883        52 3.33263e-01 5.10498e-01
336 hsa05417_Lipid_and_a..                 0         -0.000100034       172 9.93958e-01 1.00000e+00
337 hsa05418_Fluid_shear..                 0         -0.007349065       120 8.70802e-01 9.34587e-01
```
The p-value and adjusted p-values are shown. The criteria for a higher proportion of weight in the DE components in the gene set when compared to the background is shown in the Prop.DE.Increased column. In order for a gene set to be consider enriched for DE genes, the p-value should be significant and the proportion of DE genes should increased in comparison to the background. 

&nbsp;


&nbsp;

