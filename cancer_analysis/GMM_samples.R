GMM_side_by_side <- function(postdata, raw.gs, set_number = NULL, DF = NULL, params = NULL, threshold = NULL){
  
  if(is.null(set_number) == FALSE){
    set_name <- names(raw.gs[set_number]) 
    raw.gs <- raw.gs[which(names(raw.gs) == set_name)]}
  
  ## thrshold max for middle component
  comp1_var_max <- seq(0,1, by = 0.00001)
  comp1_var_max <- comp1_var_max[which(pnorm(log2(threshold), 0, sqrt(comp1_var_max)) < 0.975)[1]]
  
  v<- c()
  nonDE_criterion<- c()
  tol<- c()
  BIC<- c()
  nonDE_criterion_upper_CI <- c()
  for(j in 1:length(raw.gs)){
    tryCatch({
  postdata_set <- postdata[which(postdata$genesEnsembl %in% raw.gs[[j]]),] 
  
  effect <- c(postdata_set$GROUP_effect_1, postdata_set$GROUP_effect_2)
  variance <- c(postdata_set$GROUP_variance_1, postdata_set$GROUP_variance_2)
  
  
  
  
  # fit ggm to all genes without regard for set membership
  all_genes_mixture <- .EM_7FP_fixed(effect, variance, comp1_var_max =  comp1_var_max, threshold = threshold)
  loglike_all_genes <- all_genes_mixture$loglike
  #print(all_genes_mixture$param)
 

  if(DF == 6){
 
    # for each gene set in the list a gmm model is fit for the genes in the set and then genes out of the set

      # same as above we create three groups out of genes in the set to initialize the EM algorithm
      set_specific_post<- postdata_set
      set_specific_post_in<- set_specific_post
      effect_inset <- set_specific_post_in$GROUP_effect_1
      variance_inset<- set_specific_post_in$GROUP_variance_1 
      
      
      # fit the gmm to the genes in the gene set
      inset_mixture <- .EM_7FP_fixed(effect_inset, variance_inset, comp1_var_max =  comp1_var_max, threshold = threshold)
      loglike_Inset_genes <- inset_mixture$loglike
      inset_parameters <- inset_mixture$param
      
      # get all genes in the outset
      set_specific_post_out<- set_specific_post 
      effect_outset <- set_specific_post_out$GROUP_effect_2
      variance_outset<- set_specific_post_out$GROUP_variance_2 
      
      # fit the gmm to genes outside the gene set
      outset_mixture <- .EM_7FP_fixed(effect_outset, variance_outset, comp1_var_max =  comp1_var_max, threshold = threshold)
      loglike_Outset_genes <- outset_mixture$loglike
      outset_parameters <- outset_mixture$param
      
      parameters_h1 <- c(inset_parameters, outset_parameters)
      # compare the two models
      aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 6)
      BIC_all<- 6*log(nrow(postdata))-2*loglike_all_genes
      BIC_set<- 12*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
      v[j]<- aa$p.value[[1]]
      nonDE_criterion[j] <- inset_parameters$alpha[1] < outset_parameters$alpha[1]
      tol[j]<- length(raw.gs[[j]])
      BIC[j]<- BIC_set < BIC_all
      parameters_h1 <- c(inset_parameters, outset_parameters)
      
      
      progress(j, max.value = length(raw.gs)) }
  
      else if(DF == 1){
        
        set <- c(rep(0,length(effect)/2) , rep(1,length(effect)/2))
      
        set_mixture <- .EM_1FP_fixed(effect, variance, set, comp1_var_max, threshold = threshold)
        loglike_set_genes <- set_mixture$loglike
        set_parameters <- set_mixture$param
        ll_trace <- set_mixture$ll.vector
        
        
        
        # compare the two models
        aa<-lr.test(-loglike_all_genes, - loglike_set_genes,alpha = 0.05, 1)
        BIC_all<- 6*log(nrow(postdata))-2*loglike_all_genes
        BIC_set<- 7*log(nrow(postdata))-2*(loglike_set_genes)
        parameters_h1 <- set_parameters
        nonDE_criterion[j]<- set_parameters$alpha[1] < set_parameters$alpha[4]
        v[j]<- aa$p.value[[1]]
        tol[j]<- length(raw.gs[[j]])
        BIC[j]<- BIC_set < BIC_all
        
        # CI for middle weights only for sets with m dist true to save time
        p <- c(set_parameters$alpha[1], set_parameters$alpha[4], set_parameters$alpha[5]/(1-set_parameters$alpha[4]), set_parameters$mu[2], set_parameters$mu[3], set_parameters$var[2], set_parameters$var[3])
        hess <- hessian(func = .ll_set, x = p, lfc = effect, var_gene = variance, threshold = threshold, set = set, method="Richardson", method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
        CI <- 1.96*sqrt(diag(inv(hess)))
        nonDE_criterion_upper_CI[j] <- ((set_parameters$alpha[4] - CI[2]) - (set_parameters$alpha[1] + CI[1]))
        
        if(is.nan(nonDE_criterion_upper_CI[j]) == TRUE){
          print("BOOTSTRAPPING")
          data <- data.frame("effect" = effect, "variance" = variance, "set" = set)
          boot_out <- boot(data = data,threshold = threshold,statistic = .get_se, R = 200)
          nonDE_criterion_upper_CI[j] <- (set_parameters$alpha[4] - 1.96*sd(boot_out$t[,2])) - (set_parameters$alpha[1] + 1.96*sd(boot_out$t[,1]))
        }
        progress(j, max.value = length(raw.gs))
      }
      }, error = function(e){
      cat("Error in ", names(raw.gs[j]),"set\n")
        print(e)
      nonDE_criterion[j]<- NA
      v[j]<- NA
      tol[j]<- NA
      BIC[j]<- NA  
      nonDE_criterion_upper_CI[j] <- NA
      progress(j, max.value = length(raw.gs))
    })
     # print(inset_parameters)
     # print(outset_parameters)
    }
  

  p_adj <- p.adjust(v, method = "BH", n = length(v))
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"= nonDE_criterion, "size"=tol, "BIC_value"=BIC, "Adj-Pval" = p_adj)
  if(DF == 1) detected_gene_sets <- cbind(detected_gene_sets, nonDE_criterion_upper_CI)
  if(is.null(params) == TRUE) return(detected_gene_sets) else return(parameters_h1)
  
}