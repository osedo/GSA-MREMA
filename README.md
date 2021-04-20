# GSA-MREMA

### Scripts
1. **packages.R**- Loads the required packages. BiocManager version 3.12 needed.
2. **functions.R** - Defines functions used for the EM algorithm and for simulating data.
3. **simulateddata.R** - Functions for the simulation of a SummarizedExperiment object - one or two cell types simulations possible.
4. **globalfunctions.R** - Functions for Gene Set Anlysis and csGSA. 

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

### MREMA-GSA


#### Simulate data

```R
simulation_1 <- simulation(s=100 , p=10, N=20, beta=0.02 ,gamma=.8 , foldchange=3 , upreg=.5)
```
s - number of gene sets.  
p - number of genes in each set.  
N - number of samples.  
beta - proportion of gene sets in data set that have DE genes.  
gamma - proportion of DE genes in a gene set with DE genes.  
foldchange - foldchange of DE genes.  
upreg - proportion of DE genes that are upregulated.  

`simulation_1$sum_exp_raw_count` - SummarizedExperiment object.

`simulation_1$raw.gs` - List of gene sets. If gene sets in DE they will be last gene sets in the list. In this example gene sets 99 and 100.

&nbsp;

#### RUN MREMA-GSA

```R
w_o_mrema(simulation_1$sum_exp_raw_count, simulation_1$raw.gs)
```

`mrema()` and `wm_o_mrema()` can be used here for MREMA and WM-MREMA respectively. 

The p-value is given for the likelihood-ratio test. This should be corrected for multiple testing and is conditional on the mid_dist being TRUE. mid_dist tests whether the weight assigned to the non-DE genes in a gene set is smaller in the gene set than in the genes outside the gene set.
The size of the gene set is also given.

Finally it is tested whether the Bayesian information criterion (BIC) of the gene is lower than the BIC for all.


&nbsp;

&nbsp;

### csGSA (cell type specific GSA)

#### Simulate data
```R
prop = rbeta(20,1,1)
simulation_cancer <- simulation_2Cell_types(s=100 , p=10, N=20, beta=0.02 ,gamma=.8 , foldchange=3 , upreg=.5, prop_norm = prop)
```
Same parameter as above for DE in cell type one. The proportions for cell type two (with no DE genes or sets) are given as an arguement.

Counts for genes in cell type two are drawn from the same distribution as their corresponding gene in cell type one controls, but no DE is seen in cases.

The gene set list returned is based on cell type one, as above if their is DE the last gene sets listed will contain the DE genes.

#### Cancer_WO_MREMA
```R
cancer_w_o_mrema(simulation_cancer$sum_exp_raw_count, simulation_cancer$raw.gs)
```
&nbsp;

Compare with genes sets found to be enriched in normal cells.

```R
normal_w_o_mrema(simulation_cancer$sum_exp_raw_count, simulation_cancer$raw.gs)
```





