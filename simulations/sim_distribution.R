

# simulation to get a gene set and a SummarizedExperiment object
# from this simulation we get a summarized espression matrix and a gene set
#s<- number of gene sets
#p<- number of genes in each set
#N<- number of samples
#Beta<- the proportion of gene sets in the data set that are truly enriched for DE genes
#gamma<- the percentage of genes, truly DE in each gene set....set to 1 for our simulations
#theta <- proportion of DE genes in non-Enriched gene sets...set to 1 for our simulations
#upreg<- the proportion of upregulated gene in de genes...set to 0.5 for bi-directional enrichment 
#foldchange<- the fold change, a list, with each item containing values for DE genes in each DE gene set





simulation_dist <- function(s,p,N,beta,gamma, theta, foldchange ,upreg, inset_fc){
  
 # j <- read.csv("BRCA_sim_param.csv")
  mean_para <- j$mean_para
  disp_para <- j$disp_para
  
  raw.gs <- list()

  fc_DE <- 2^rhalfnorm(9900*theta, sd2theta(inset_fc))

  ### Here we simulate all nonDE sets together, pi0 = 1 ..... no DE genes
  
  simulation_nonDE <-sim.counts(nGenes = 9900, pi0=(1-theta), m=N/2, mu=mean_para, disp = disp_para,
                                fc=fc_DE, up=0.5, replace = TRUE)
  
  nonDE_counts <- simulation_nonDE$counts
  
  delta <- simulation_nonDE$delta
  
  #### for each DE set we call sim.counts() allowing us to give each DE set the DE structure given in the foldchange list. 
  #### we then add this to the expression matrix created above. 
  if(beta == 1) total_counts <- 0 else total_counts <- nonDE_counts
  GeneSetDE <- c(rep(1,s*beta))
  
  for(i in 1:length(GeneSetDE)){
    fc <- foldchange
    sim_DE <- sim.counts(nGenes = p, pi0=(1-gamma), m=N/2, mu=mean_para, disp = disp_para, fc=fc, up=upreg, replace = FALSE)
    DE_counts <- sim_DE$counts
    total_counts <- rbind(total_counts, DE_counts)
    delta <- c(delta,sim_DE$delta)
  }
  
  if(beta == 1) total_counts <- total_counts[-1,]
  
  #### sets groups to 0 and 1
  if(beta == 0) pre1_group<- simulation_nonDE$group else pre1_group <- sim_DE$group
  pre2_group<- replace(pre1_group, pre1_group==1,0)
  group<- replace(pre2_group, pre2_group==2,1)
  rm(pre2_group)
  rm(pre1_group)
  #### create summarizedexperiment object
  rowData<- data.frame(row.names = paste0("gene", 1:(s*p)))
  colData <- data.frame(GROUP=group, row.names=paste0("sample", 1:N))
  sum_exp_raw_count<- SummarizedExperiment(assays=SimpleList(counts=total_counts), colData=colData, rowData=rowData)
  
  ### create gene sets 
  #### just split by the number of genes in each set, this returns the correct amount of sets with the DE structure specified
  #### it should be easy enough to change this if we wanted uneven set size? 
 # genes <- rownames(rowData)
#  raw.gs[[1]] <- genes[c(1:(length(genes)-p))]
  #raw.gs[[2]] <- genes[c((length(genes)-p+1):length(genes))]
  
raw.gs <- split(rownames(rowData), ceiling(seq_along(rownames(rowData))/p))
names(raw.gs) <- paste0("gs",(1:s))
  
#  names(raw.gs) <- paste0("gs",(1:2))
  
  ### saving the delta data from DE count
 
 
  simulateddata<- list("raw.gs"=raw.gs, "sum_exp_raw_count"=sum_exp_raw_count, "LFC" = delta)
  
  return(simulateddata)
}
