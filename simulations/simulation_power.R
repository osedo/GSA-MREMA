
simulation_power <- function(inset_LFC, outset_LFC, n_samples, means, disp, outset_split){
  # inset_LFC - the lfc distribution of the enriched gene set
  # outset_LFC - the lfc distribution of the outset....this is broken into different equal size sets to use enrichemnt browser on.
  # n_samples - total number of samples, spit equally into cases and controls
  # means - mean values for genes
  # disp - dispersion parameter for genes
  # outset_split - number of genes in each set of the outset
  
  
  # simulate the outset - number of genes is the number of LFC values given in outset_LFC
  non_enriched <- matrix(0, nrow = length(outset_LFC), ncol = n_samples)
  case_FC <- 2^outset_LFC
  sampled <- sample(c(1:length(means)), length(outset_LFC), replace = TRUE)
  non_enriched_means <- means[sampled]
  non_enriched_disp <- disp[sampled]
  for(i in 1:length(case_FC)){
    controls <- rnegbin(n_samples/2, mu = non_enriched_means[i], theta = 1/non_enriched_disp[i])
    cases <- rnegbin(n_samples/2, mu = non_enriched_means[i]*case_FC[i], theta = 1/non_enriched_disp[i])
    non_enriched[i,] <- c(controls, cases)
  }
  
  
  # simulate the inset - number of genes is the number of LFC values given in inset_LFC
  enriched <- matrix(0, nrow = length(inset_LFC), ncol = n_samples)
  case_FC <- 2^inset_LFC
  sampled <- sample(c(1:length(means)), length(inset_LFC), replace = TRUE)
  enriched_means <- means[sampled]
  enriched_disp <- disp[sampled]
  for(i in 1:length(case_FC)){
    controls <- rnegbin(n_samples/2, mu = enriched_means[i], theta = 1/enriched_disp[i])
    cases <- rnegbin(n_samples/2, mu = enriched_means[i]*case_FC[i], theta = 1/enriched_disp[i])
    enriched[i,] <- c(controls, cases)
  }
  
  #combine into one expression matrix
  counts <- rbind(non_enriched, enriched)
  rownames(counts) <- paste0("gene",1:dim(counts)[1])
  colnames(counts) <- paste0("sample",1:dim(counts)[2])
  sample_data <- data.frame("GROUP" = as.factor(c(rep(0,dim(counts)[2]/2), rep(1,dim(counts)[2]/2))))
  gene_data <- data.frame(row.names = paste0("gene", 1:dim(counts)[1]))
  rownames(sample_data) <- paste0("sample",1:dim(counts)[2])
  
  #make a summarized experiment object to store
  experiment <- SummarizedExperiment(assays=SimpleList(counts=counts), colData=sample_data, rowData=gene_data)
  
  # default to making outset into 99 approx equal sized sets
  gene_set <- split(rownames(gene_data)[1:length(outset_LFC)], ceiling(seq_along(rownames(gene_data)[1:length(outset_LFC)])/outset_split))
  names(gene_set) <- paste0("gs",(1:length(gene_set)))
  gene_set[[(length(gene_set) +1)]] <- rownames(gene_data)[(length(outset_LFC) + 1):(length(outset_LFC) + length(inset_LFC))]
  names(gene_set)[length(gene_set)] <- paste0("gs", length(gene_set))
  list("sum_exp_raw_count" = experiment, "raw.gs" = gene_set)
}