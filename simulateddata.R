# simulation to get a gene set and a SummarizedExperiment object
# from this simulation we get a summarized espression matrix and a gene set
#s<- number of gene sets
#p<- number of genes in each set
#N<- number of samples
#Beta<- the proportion of gene sets in the data set that have truly DE genes
#gamma<- the percentage of genes, truly DE in each gene set
#upreg<- the proportion of upregulated gene in de genes
#foldchange<- the fold change




#setwd("/data/adib/metapower1")
source("packages.R")
source("functions.R")
simulation<- function(s,p,N,beta,gamma,foldchange,upreg){
  
# the idea is to generate gene expression sampled from negative binomial distribution, 
# therefore, we require a vector of means and a vector of dispersion parameters, 
# we use Airway data set to obtain those vectors
  
  data(airway)
  
  air.eset <- as(airway, "ExpressionSet")
  annotation(air.eset) <- "hsa"
  air.eset <- air.eset[grep("^ENSG", rownames(air.eset)), ]
  dim(air.eset)
  # there are two groups of samples 
  pData(air.eset)$GROUP <- ifelse(colData(airway)$dex == "trt", 1, 0)
  pData(air.eset)$BLOCK <- colData(airway)$cell
  table(pData(air.eset)$BLOCK)
  
  
  rseqSE<- deAna(exprs(air.eset),grp =  pData(air.eset)$GROUP , de.method = "DESeq2")
  rseqSE_count<- exprs(air.eset)[rownames(rseqSE),]
  deairway<- DESeqDataSetFromMatrix(countData = rseqSE_count, colData = pData(air.eset),
                                    design = ~ dex)
  
  deairway<- DESeq(deairway)
  dis_obj<- estimateDispersions(deairway)
  disp_para<- dispersions(dis_obj)
  # this disperation vector has NA value due to removing genes by DEseq2, we have to remove them 
  # and also removing the same genes from the mean count vector
  ###########
  disp_para<-disp_para[!is.na(disp_para)]
  # now we find the mean count vector, mean count of each gene, and then we romve the filtered genes
  mean_para<-apply(rseqSE_count, 1, mean)
  ###########
  mean_para<- as.vector(mean_para)
  mean_para<-mean_para[mean_para>200]
  ######################
  # now we use sim.counts function from ssizeRNA package to generate genes from negative binomial 
  # distribution with mean_i=mean_para_i and disperation_i=disp_para_i, for each gene i
  
  simulation<-sim.counts(nGenes = s*p, pi0=1-gamma*beta, m=N/2, mu=mean_para, disp = disp_para,
                         fc=foldchange, up=upreg, replace = TRUE) 
  raw_countmat<- simulation$counts
  pre1_group<- simulation$group
  pre2_group<- replace(pre1_group, pre1_group==1,0)
  group<- replace(pre2_group, pre2_group==2,1)
  # now we must convert this raw count into a SummarizedExperiment object
  nrow<- s*p
  ncols <- N
  counts <- raw_countmat
  rowData<- DataFrame(row.names = paste0("gene", 1:(s*p)))
  
  colData <- DataFrame(GROUP=group,
                       row.names=paste0("sample", 1:N))
  
  sum_exp_raw_count<- SummarizedExperiment(assays=SimpleList(counts=counts), colData=colData, rowData=rowData)
  if(gamma==0 & beta==0){
    
    raw.gs<- makeExampleGeneSets(gnames=rownames(assay(sum_exp_raw_count)), n=s, size=p)
  }else if(beta==1){
    vec<- rownames(assay(sum_exp_raw_count))
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(upDEgenes, (gamma*p)/2), sample(downDEgenes,(gamma*p)/2))
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][1:(gamma*p)/2]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][1:(gamma*p)/2]]
    }
    
    names(DEsets) <- paste0("gs", 1:s)
    raw.gs<- DEsets 
  }else if(upreg==1){
    vec<- rownames(assay(sum_exp_raw_count))
    # remember the number of non-DE genes is ps-Bs(gp)
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
    # the number of nonDE sets is (1-B)s and there are p(1-B)s genes in them,
    # so first from all the non DE genes I am going to create (1-B)s sets each of size p
    
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(upDEgenes, (gamma*p)))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
    
    
  }else if(upreg==0){
    vec<- rownames(assay(sum_exp_raw_count))
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
   
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(downDEgenes,(gamma*p)))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
  }else{
    vec<- rownames(assay(sum_exp_raw_count))
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
    
    
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(upDEgenes, (gamma*p)/2), sample(downDEgenes,(gamma*p)/2))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
    
    
    
    
    
  }
  
  
  
  simulateddata<- list("raw.gs"=raw.gs, "sum_exp_raw_count"=sum_exp_raw_count)
  return(simulateddata)
}


#tt<- simulation(10,10,8,0.2,0.8,3,0.5)

#sum_exp_raw_count<- tt$sum_exp_raw_count
#raw.gs<- tt$raw.gs
#save(sum_exp_raw_count, file="sum_exp_raw_count")
#save(raw.gs, file = "raw.gs")

simulation_2Cell_types<- function(s,p,N,beta,gamma,foldchange,upreg, prop_norm){
  
  # the idea is to generate gene expression sampled from negative binomial distribution, 
  # therefore, we require a vector of means and a vector of dispersion parameters, 
  # we use Airway data set to obtain those vectors
  
  data(airway)
  
  air.eset <- as(airway, "ExpressionSet")
  annotation(air.eset) <- "hsa"
  air.eset <- air.eset[grep("^ENSG", rownames(air.eset)), ]
  dim(air.eset)
  # there are two groups of samples 
  pData(air.eset)$GROUP <- ifelse(colData(airway)$dex == "trt", 1, 0)
  pData(air.eset)$BLOCK <- colData(airway)$cell
  table(pData(air.eset)$BLOCK)
  
  
  rseqSE<- deAna(exprs(air.eset),grp =  pData(air.eset)$GROUP , de.method = "DESeq2")
  rseqSE_count<- exprs(air.eset)[rownames(rseqSE),]
  deairway<- DESeqDataSetFromMatrix(countData = rseqSE_count, colData = pData(air.eset),
                                    design = ~ dex)
  
  deairway<- DESeq(deairway)
  dis_obj<- estimateDispersions(deairway)
  disp_para<- dispersions(dis_obj)
  # this disperation vector has NA value due to removing genes by DEseq2, we have to remove them 
  # and also removing the same genes from the mean count vector
  ###########
  disp_para<-disp_para[!is.na(disp_para)]
  
  
  
  # now we find the mean count vector, mean count of each gene, and then we romve the filtered genes
  mean_para<-apply(rseqSE_count, 1, mean)
  ###########
  mean_para<- as.vector(mean_para)
  mean_para<-mean_para[mean_para>200]
  ######################
  # now we use sim.counts function from ssizeRNA package to generate genes from negative binomial 
  # distribution with mean_i=mean_para_i and disperation_i=disp_para_i, for each gene i
  
  simulation<-sim.counts(nGenes = s*p, pi0=1-gamma*beta, m=N/2, mu=mean_para, disp = disp_para,
                         fc=foldchange, up=upreg, replace = TRUE)
  simulation2 <- sim.counts(nGenes = s*p, pi0=1, m=N/2, mu=simulation$lambda0, disp = simulation$phi0,
                            fc=0, up=0, replace = FALSE)
  
  sum_exp_raw_count_cancer<- simulation$counts
  
  sum_exp_raw_count_normal<- simulation2$counts
  
  
  
  cancer.matrix<- sum_exp_raw_count_cancer
  
  normal.matrix<- sum_exp_raw_count_normal
  
  #purity<- rbeta(N,shape1, shape2)
  
  purity.cancer<- prop # it shows the proportion of normal cells
  
  purity.normal<- 1-purity.cancer
  
  
  
  # multiply the expression value by purity
  
  post.cancer.matrix<- t(apply(cancer.matrix, 1, function(x) x*purity.normal))
  
  post.normal.matrix<- t(apply(normal.matrix, 1, function(x) x*purity.cancer))
  
  
  
  # add up the two matrices to obtain the weighted sum
  
  mix.matrix<- post.cancer.matrix + post.normal.matrix
  
  total_counts <- round(mix.matrix)
  total_counts <- round(total_counts)
  raw_countmat<- total_counts
  pre1_group<- simulation$group
  pre2_group<- replace(pre1_group, pre1_group==1,0)
  group<- replace(pre2_group, pre2_group==2,1)
  # now we must convert this raw count into a SummarizedExperiment object
  nrow<- s*p
  ncols <- N
  counts <- raw_countmat
  rowData<- DataFrame(row.names = paste0("gene", 1:(s*p)))
  
  colData <- DataFrame(GROUP=group,
                       row.names=paste0("sample", 1:N), purity = prop)
  
  sum_exp_raw_count<- SummarizedExperiment(assays=SimpleList(counts=counts), colData=colData, rowData=rowData)
  if(gamma==0 & beta==0){
    
    raw.gs<- makeExampleGeneSets(gnames=rownames(assay(sum_exp_raw_count)), n=s, size=p)
  }else if(beta==1){
    vec<- rownames(assay(sum_exp_raw_count))
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(upDEgenes, (gamma*p)/2), sample(downDEgenes,(gamma*p)/2))
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][1:(gamma*p)/2]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][1:(gamma*p)/2]]
    }
    
    names(DEsets) <- paste0("gs", 1:s)
    raw.gs<- DEsets 
  }else if(upreg==1){
    vec<- rownames(assay(sum_exp_raw_count))
    # remember the number of non-DE genes is ps-Bs(gp)
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
    # the number of nonDE sets is (1-B)s and there are p(1-B)s genes in them,
    # so first from all the non DE genes I am going to create (1-B)s sets each of size p
    
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(upDEgenes, (gamma*p)))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
    
    
  }else if(upreg==0){
    vec<- rownames(assay(sum_exp_raw_count))
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
    
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(downDEgenes,(gamma*p)))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
  }else{
    vec<- rownames(assay(sum_exp_raw_count))
    nonDE<- which(simulation$de==0)        
    nonDEgene<- vec[nonDE]
    
    
    nonDEsets<- DisSample(nonDEgene, (1-beta)*s, p)
    names(nonDEsets) <- paste0("gs", seq_len((1-beta)*s))
    rr<- unlist(nonDEsets)
    rrr<- as.vector(rr)
    nonDEgenesInDEsets<- nonDEgene[!nonDEgene %in% rrr]
    
    upDE<- which(simulation$de==1)        
    upDEgenes<- vec[upDE]
    downDE<- which(simulation$de==-1)   
    downDEgenes<- vec[downDE]
    
    
    DEsets<- list()
    fer<- beta*s
    for(i in 1:fer){
      DEsets[[i]]<- c(sample(nonDEgenesInDEsets, (p-p*gamma)), 
                      sample(upDEgenes, (gamma*p)/2), sample(downDEgenes,(gamma*p)/2))
      nonDEgenesInDEsets<- nonDEgenesInDEsets[!nonDEgenesInDEsets %in% DEsets[[i]][1:((1-gamma)*p)]]
      upDEgenes<- upDEgenes[!upDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
      downDEgenes<- downDEgenes[!downDEgenes %in% DEsets[[i]][((1-gamma)*p):p]]
    }
    
    names(DEsets) <- paste0("gs", (((1-beta)*s)+1):s)
    raw.gs<- append(nonDEsets,DEsets)  
    
    
    
    
    
  }
  
  
  
  simulateddata<- list("raw.gs"=raw.gs, "sum_exp_raw_count"=sum_exp_raw_count)
  return(simulateddata)
}

