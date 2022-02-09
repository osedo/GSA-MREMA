mrema <- function(postdata, raw.gs, set_number = NULL, DF = NULL){
  
  effect <- pull(postdata,2)

  variance<- pull(postdata,3)
  if(is.null(set_number) == FALSE){
    set_name <- names(raw.gs[set_number]) 
    raw.gs <- raw.gs[which(names(raw.gs) == set_name)]}

  ## create three groups out of the genes to provide initial parameters for EM algorithm
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  effect_df <- dplyr::data_frame(x = effect, cluster = effect_kmeans_cluster)
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = dplyr::n())
  
  # Set one component with a mean of zero
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # order the remaining two components
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2 <- effect_summary_df[2,]
    m3 <- effect_summary_df[3,]
    effect_summary_df[2,] <- m3
    effect_summary_df[3,] <- m2
  }
  
  # generate the initial mixing weights
  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # fit ggm to all genes without regard for set membership
  all_genes_mixture <- .EM_7FP(effect, variance, effect_summary_df)
  loglike_all_genes <- all_genes_mixture$loglike

  initial_alpha<- c( effect_summary_df[["alpha"]][1],effect_summary_df[["alpha"]][2],
                     effect_summary_df[["alpha"]][3],effect_summary_df[["alpha"]][1],
                     effect_summary_df[["alpha"]][2],effect_summary_df[["alpha"]][3])
  
  v<- c()
  nonDE_criterion<- c()
  tol<- c()
  BIC<- c()
  
  ### if we want to fit the model with 7 degrees of freedom between the null and alternative model 
  if(DF == 7){
  # for each gene set in the list a gmm model is fit for the genes in the set and then genes out of the set
  for(j in 1:length(raw.gs)){
    tryCatch({
    # same as above we create three groups out of genes in the set to initialize the EM algorithm
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(pull(postdata,1) %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    effect_inset <- pull(set_specific_post_in,2)
    variance_inset<- pull(set_specific_post_in,3) 
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
        summarize(mu = mean(x), variance = var(x), std = sd(x), size = dplyr::n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    # fit the gmm to the genes in the gene set
    inset_mixture <- .EM_7FP(effect_inset, variance_inset, effect_summary_df_inset)
    loglike_Inset_genes <- inset_mixture$loglike
    inset_parameters <- inset_mixture$param
    
    # get all genes in the outset
    set_specific_post_out<- set_specific_post %>% filter(set==0)
    effect_outset <- pull(set_specific_post_out,2)
    variance_outset<- pull(set_specific_post_out,3)
    
    # fit the gmm to genes outside the gene set
    outset_mixture <- .EM_7FP(effect_outset, variance_outset, effect_summary_df)
    loglike_Outset_genes <- outset_mixture$loglike
    outset_parameters <- outset_mixture$param
    
    
    # compare the two models
    aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 7)
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 14*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
    v[j]<- aa$p.value[[1]]
    nonDE_criterion[j] <- inset_parameters$alpha[1] < outset_parameters$alpha[1]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
    }, error = function(e){
      cat("Error in ", names(raw.gs[j]),"set\n")
      cat("Check set overlap with expression matrix: Overlap = ",length(which(set_specific_post$set == 1)) , "\n")
      nonDE_criterion[j]<- NA
      v[j]<- NA
      tol[j]<- NA
      BIC[j]<- NA  
    }) 
   }
  } else if(DF == 1) {
    ### if we want to fit the model with 1 degree of freedom between the null and alternative model 
    for(j in 1:length(raw.gs)){
      tryCatch({
      set_specific_post<- postdata
      set_specific_post$set<- ifelse(pull(postdata,1) %in% raw.gs[[j]], 1, 0)
      effect_set <- pull(set_specific_post,2)
      variance_set<- pull(set_specific_post,3)
      set<- set_specific_post$set
      
      
      
      set_mixture <- .EM_1FP(effect_set, variance_set, effect_summary_df, set, initial_alpha)
      
      loglike_set_genes <- set_mixture$loglike
      set_parameters <- set_mixture$param
      
      aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 1)
      
      BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
      BIC_set<- 8*log(nrow(postdata))-2*(loglike_set_genes)
      
      
      nonDE_criterion[j]<- set_parameters$alpha[1] < set_parameters$alpha[4]
      v[j]<- aa$p.value[[1]]
      tol[j]<- length(raw.gs[[j]])
      BIC[j]<- BIC_set < BIC_all 
      }, error = function(e){
        cat("Error in ", names(raw.gs[j]),"set\n")
        cat("Check set overlap with expression matrix: Overlap = ", length(which(set == 1)), "\n")
        nonDE_criterion[j]<- NA
        v[j]<- NA
        tol[j]<- NA
        BIC[j]<- NA  
      })

    }
  }
  
  p_adj <- p.adjust(v, method = "BH", n = length(v))
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"= nonDE_criterion, "size"=tol, "BIC_value"=BIC, "Adj-Pval" = p_adj)
  return(detected_gene_sets) 
  
}



# fitting the ggm with all parameters (apart from nonDE component mean) free 
.EM_7FP <- function(effect, variance, effect_summary_df){
  n<- 1000
  m<- 1e-6
  iter<- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step <- .e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]], effect_summary_df[["alpha"]])
      m.step <- .m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- .e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], m.step[["alpha"]])
      m.step <- .m_step_iter(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]])
      loglik.vector <- c(loglik.vector, e.step[["loglik"]])
      loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
      if(loglik.diff < m) {
        break
      } else {
        cur.loglik <- e.step[["loglik"]]
      }
    }
  }  
  loglike_all_genes<-tail(loglik.vector, n=1)
  parameters <- list("loglike" = loglike_all_genes, "param" = m.step)
  return(parameters)
}


.EM_1FP <- function(effect_set, variance_set, effect_summary_df, set, initial_alpha){
  n<- 1000
  m<- 1e-6
  iter<- 10
for (i in 1:n) {
  if (i == 1) {
    # Initialization
    e.step.set <- .e_step_set_iter(effect_set, variance_set, set, effect_summary_df[["mu"]], effect_summary_df[["variance"]], initial_alpha)
    
    m.step.set <- .m_step_set_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
    cur.loglik.set <- e.step.set[["loglik"]]
    loglik.vector.set <- e.step.set[["loglik"]]
  } else {
    # Repeat E and M steps till convergence
    e.step.set <- .e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]],m.step.set[["alpha"]])
    
    m.step.set <- .m_step_set_iter(effect_set,variance_set,set,m.step.set[["var"]],iter,  e.step.set[["posterior_df"]])
    loglik.vector.set <- c(loglik.vector.set, e.step.set[["loglik"]])
    
    loglik.diff.set <- abs((cur.loglik.set - e.step.set[["loglik"]]))
    
    if(loglik.diff.set < m) {
      break
    } else {
      cur.loglik.set <- e.step.set[["loglik"]]
    }
  }
}  

loglike_set_genes<-tail(loglik.vector.set, n=1)
parameters <- list("loglike" = loglike_set_genes, "param" = m.step.set)
return(parameters)
}


# the e_step and m_step functions in the EM algorithm
.e_step_iter <- function(x,lfc_se, mu_vector, component_var, alpha_vector) {
  # both lfc_se and component_var contribute to total variance 
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_se)) * alpha_vector[1]
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_se)) * alpha_vector[2]
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_se)) * alpha_vector[3]
  
  sum_of_comps <- comp1_prod + comp2_prod + comp3_prod
  comp1_post <- comp1_prod / sum_of_comps
  comp2_post <- comp2_prod / sum_of_comps
  comp3_post <- comp3_prod / sum_of_comps
  
  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)
  
  list("loglik" = sum_of_comps_ln_sum,
       "posterior_df" = cbind(comp1_post, comp2_post, comp3_post),
       "prod_df" = cbind(comp1_prod, comp2_prod, comp3_prod))
}


# the e_step and m_step functions in the EM algorithm for w_o_mrema (weights only)
.e_step_set_iter <- function(x,lfc_se,set, mu_vector, component_var, alpha_vector) {
  # both sd_vector and sigma are variance 
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_se)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_se)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_se)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_se)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_se)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_se)) * alpha_vector[6]*(1-set)
  
  sum_of_comps1 <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps2 <- comp4_prod + comp5_prod + comp6_prod
  sum_of_comps<- sum_of_comps1 + sum_of_comps2
  
  
  comp1_post <- comp1_prod / sum_of_comps1
  comp1_post[is.na(comp1_post)] <- 0
  
  comp2_post <- comp2_prod / sum_of_comps1
  comp2_post[is.na(comp2_post)] <- 0
  
  comp3_post <- comp3_prod / sum_of_comps1
  comp3_post[is.na(comp3_post)] <- 0
  
  comp4_post <- comp4_prod / sum_of_comps2
  comp4_post[is.na(comp4_post)] <- 0
  
  comp5_post <- comp5_prod / sum_of_comps2
  comp5_post[is.na(comp5_post)] <- 0
  
  comp6_post <- comp6_prod / sum_of_comps2
  comp6_post[is.na(comp6_post)] <- 0

  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)
  
  list("loglik" = sum_of_comps_ln_sum,
       "posterior_df" = cbind(comp1_post, comp2_post, comp3_post, comp4_post, comp5_post, comp6_post))
}



.m_step_iter <- function(x, lfc_se, component_var, t, posterior_df) {
  
  
  comp1_n <- sum(posterior_df[, 1]) 
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  
  
  ###########################
  for (i in 1:t) {     
    if (i == 1) { 
      # Initialization
      comp2_var <- component_var[2]
      w2_i <- 1/(comp2_var+lfc_se)
      comp2_mu <- max(sum(posterior_df[,2] * w2_i * x)/sum(posterior_df[,2] * w2_i), log2(1.1))
      comp3_var <- component_var[3]
      w3_i <- 1/(comp3_var+lfc_se)
      comp3_mu <- min(sum(posterior_df[,3] * w3_i * x)/sum(posterior_df[,3] * w3_i), -log2(1.1))
      comp1_var <- min(component_var[1],0.05)
      w1_i <- 1/(comp1_var+lfc_se)
      comp1_mu<- 0
      
      
      
    } else {
      
      comp2_var <- max(0.005,sum(posterior_df[, 2] * ((w2_i)^2) * (((x - comp2_mu)^2)-lfc_se))/(sum(posterior_df[, 2]*w2_i^2)))
      w2_i <- 1/(comp2_var+lfc_se)
      comp2_mu <- max( sum(posterior_df[, 2] * w2_i * x)/sum(posterior_df[, 2] * w2_i), log2(1.1))
      
      comp3_var <- max(0.005,sum(posterior_df[, 3] * ((w3_i)^2) * (((x - comp3_mu)^2)-lfc_se))/(sum(posterior_df[, 3]*w3_i^2)))
      w3_i <- 1/(comp3_var+lfc_se)
      comp3_mu <- min(sum(posterior_df[, 3] * w3_i * x)/sum(posterior_df[, 3] * w3_i), -log2(1.1))
      
      comp1_var <-  min(max(0.005,sum(posterior_df[, 1] * ((w1_i)^2) * (((x - comp1_mu)^2)-lfc_se))/(sum(posterior_df[, 1]*w1_i^2))),0.05)
      w1_i <- 1/(comp1_var+lfc_se)
      comp1_mu<- 0
      
    }
  }  
  

  comp1_alpha <- max(comp1_n / length(x),lower_bound)
  comp2_alpha <- max(comp2_n / length(x),lower_bound)
  comp3_alpha <- max(comp3_n / length(x),lower_bound)
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha))
}


.m_step_set_iter <- function(x, lfc_se, set, component_var, t, posterior_df) {
  
  
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4]) 
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])
  
  comp1_alpha <- max(comp1_n / sum(set),lower_bound)
  comp4_alpha <- max(comp4_n / (length(set)-sum(set)),lower_bound)
  
  c <- (comp2_n + comp5_n)/(comp2_n + comp5_n + comp3_n + comp6_n)
  
  comp2_alpha <- (1 - comp1_alpha)*(c)
  comp3_alpha <- (1 - comp1_alpha)*(1-c)
  
  comp5_alpha <- (1 - comp4_alpha)*(c)
  comp6_alpha <- (1 - comp4_alpha)*(1-c)
  
  
  #### doesn't speed up or improve fit
  #posterior_df[,2] <- (posterior_df[, 1]) *(c)
  #posterior_df[,3] <- (posterior_df[, 2]) *(1-c)
  # posterior_df[,5] <- (posterior_df[, 5]) *(c)
  #  posterior_df[,6] <- (posterior_df[, 6]) *(1-c)
  # the following for loop iterates between mean and variance, the number of iteration is t, 
  comp1_pd <- posterior_df[,1] + posterior_df[,4]
  comp2_pd <- posterior_df[,2] + posterior_df[,5]
  comp3_pd <- posterior_df[,3] + posterior_df[,6]
  ###########################
  for (i in 1:t) {

    if (i == 1) {       
      comp2_var <- component_var[2]
      w2_i <- 1/(comp2_var+lfc_se)
      comp2_mu <- max(sum(comp2_pd * w2_i * x)/sum(comp2_pd * w2_i), log2(1.1))
      comp3_var <- component_var[3]
      w3_i <- 1/(comp3_var+lfc_se)
      comp3_mu <- min(sum(comp3_pd * w3_i * x)/sum(comp3_pd * w3_i), -log2(1.1))
      comp1_var <- min(component_var[1],0.05)
      w1_i <- 1/(comp1_var+lfc_se)
      comp1_mu<- 0
      
      
      
    } else {
      comp2_var <- max(0.005,sum(comp2_pd * ((w2_i)^2) * (((x - comp2_mu)^2)-lfc_se))/(sum(comp2_pd*w2_i^2)))
      w2_i <- 1/(comp2_var+lfc_se)
      comp2_mu <- max( sum(comp2_pd * w2_i * x)/sum(comp2_pd * w2_i), log2(1.1))
      
      comp3_var <- max(0.005,sum(comp3_pd * ((w3_i)^2) * (((x - comp3_mu)^2)-lfc_se))/(sum(comp3_pd*w3_i^2)))
      w3_i <- 1/(comp3_var+lfc_se)
      comp3_mu <- min(sum(comp3_pd * w3_i * x)/sum(comp3_pd * w3_i), -log2(1.1))
      
      comp1_var <-  min(max(0.005,sum(comp1_pd * ((w1_i)^2) * (((x - comp1_mu)^2)-lfc_se))/(sum(comp1_pd*w1_i^2))),0.05)
      w1_i <- 1/(comp1_var+lfc_se)
      comp1_mu<- 0

    }
  }  
  
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha ))
}


lower_bound<- 0.000000005
