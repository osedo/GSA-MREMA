## postdata: a dataframe with gene names, effect size estimate and variance of the estimate as columns, in that order
## raw.gs: a list of gene sets
## DF: degrees of freedom, 1 or 6
## threshold: the soft fold-change threshold chosen for DE
## set_number: run test on one gene set only
## params: set to TRUE to return parameters of distribution, used for visualisation along with plot_distributions() function
## CI_interval; set to false to stop calculation of lower bound of difference in weights.....useful for speeding up simulations

mrema <- function(postdata, raw.gs, DF = NULL, threshold = NULL, set_number = NULL, params = NULL, CI_interval = TRUE){
  
  tictoc::tic("\nTime to run: ")
  postdata <- postdata[complete.cases(postdata),]
  effect <- dplyr::pull(postdata,2)
  variance<- dplyr::pull(postdata,3)
  if(is.null(set_number) == FALSE){
    set_name <- names(raw.gs[set_number]) 
    raw.gs <- raw.gs[which(names(raw.gs) == set_name)]}
  
  
  
  
  ## threshold max for middle component
  comp1_var_max <- seq(0,1, by = 0.00001)
  comp1_var_max <- comp1_var_max[which(pnorm(log2(threshold), 0, sqrt(comp1_var_max)) < 0.975)[1]]
  
  
  # fit ggm to all genes without regard for set membership
  all_genes_mixture <- .EM_7FP_fixed(effect, variance, comp1_var_max =  comp1_var_max, threshold = threshold)
  loglike_all_genes <- all_genes_mixture$loglike
  
  
  v<- rep(1,length(raw.gs))
  nonDE_criterion<- rep(FALSE, length(raw.gs))
  tol<- rep(0, length(raw.gs))
  BIC<- rep(FALSE, length(raw.gs))
  Upper_CI_nonDE_inset <- rep(1,length(raw.gs))
  for(j in 1:length(raw.gs)){
    tryCatch({
      if(DF == 1){
        
        ### run 1DF approach
        set_specific_post<- postdata
        set_specific_post$set<- ifelse(pull(postdata,1) %in% raw.gs[[j]], 1, 0)
        effect_set <- pull(set_specific_post,2)
        variance_set<- pull(set_specific_post,3)
        set<- set_specific_post$set
        set_mixture <- .EM_1FP_fixed(effect_set, variance_set, set, comp1_var_max, threshold = threshold)
        loglike_set_genes <- set_mixture$loglike
        set_parameters <- set_mixture$param
        ll_trace <- set_mixture$ll.vector
        
        
        
        # compare the two models
        aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 1)
        BIC_all<- 6*log(nrow(postdata))-2*loglike_all_genes
        BIC_set<- 7*log(nrow(postdata))-2*(loglike_set_genes)
        parameters_h1 <- set_parameters
        nonDE_criterion[j]<- set_parameters$alpha[1] < set_parameters$alpha[4]

        
         
          
          
          if(CI_interval == TRUE){
            # hessian(f, pars, method=“Richardson”, method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE)).
            p <- c(set_parameters$alpha[1], set_parameters$alpha[4], set_parameters$alpha[5]/(1-set_parameters$alpha[4]), set_parameters$mu[2], set_parameters$mu[3], set_parameters$var[2], set_parameters$var[3])
            hess <- hessian(func = .ll_set, x = p, lfc = postdata$effect, var_gene = postdata$variance, threshold = threshold, set = set, method="Richardson", method.args=list(eps=1e-8, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=6, v=2, show.details=FALSE))
            CI <- 1.96*sqrt(diag(inv(hess)))
            Upper_CI_nonDE_inset[j] <-(set_parameters$alpha[4] - CI[2]) - (set_parameters$alpha[1] + CI[1])
            if(is.nan(Upper_CI_nonDE_inset[j]) == TRUE){
              print("BOOTSTRAPPING") 
              data <- data.frame("effect" = effect_set, "variance" = variance_set, "set" = set)
              boot_out <- boot(data = data,threshold = threshold,statistic = .get_se, R = 200)
              Upper_CI_nonDE_inset[j] <- (set_parameters$alpha[4] - 1.96*sd(boot_out$t[,2])) - (set_parameters$alpha[1] + 1.96*sd(boot_out$t[,1]))
            }
          }
        
        v[j]<- aa$p.value[[1]]
        tol[j]<- length(raw.gs[[j]])
        BIC[j]<- BIC_set < BIC_all
        progress(j, max.value = length(raw.gs))
      } else if(DF == 6) {
        
        ### run 6DF approach
        set_specific_post<- postdata
        set_specific_post$set<- ifelse(pull(postdata,1) %in% raw.gs[[j]], 1, 0)
        set_specific_post_in<- set_specific_post %>% filter(set==1)
        effect_inset <- pull(set_specific_post_in,2)
        variance_inset<- pull(set_specific_post_in,3)
        
        # fit the gmm to the genes in the gene set
        inset_mixture <- .EM_7FP_fixed(effect_inset, variance_inset, comp1_var_max =  comp1_var_max, threshold = threshold)
        loglike_Inset_genes <- inset_mixture$loglike
        inset_parameters <- inset_mixture$param
        
        # get all genes in the outset
        set_specific_post_out<- set_specific_post %>% filter(set==0)
        effect_outset <- pull(set_specific_post_out,2)
        variance_outset<- pull(set_specific_post_out,3)
        
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
        progress(j, max.value = length(raw.gs))
        
        
      }
    }, error = function(e){
      cat("Error in ", names(raw.gs[j]),"set\n")
      cat(e)
      nonDE_criterion[j]<- NA
      v[j]<- NA
      tol[j]<- NA
      BIC[j]<- NA  
      progress(j, max.value = length(raw.gs))
      
    })
    
  }
  
  
  p_adj <- p.adjust(v, method = "BH", n = length(v))
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"= nonDE_criterion, "size"=tol, "BIC_value"=BIC, "Adj-Pval" = p_adj)
  if(CI_interval == TRUE & DF == 1){ detected_gene_sets <- cbind(detected_gene_sets,  Upper_CI_nonDE_inset)}
  toc()
  
  if(is.null(params) == TRUE) return(detected_gene_sets) else return(parameters_h1)
  
}

### function for bootstrapping
.get_se <- function(data, indices,threshold) {
  comp1_var_max <- seq(0,1, by = 0.00001)
  comp1_var_max <- comp1_var_max[which(pnorm(log2(threshold), 0, sqrt(comp1_var_max)) < 0.975)[1]]
  d <- data[indices, ]
  p <- .EM_1FP_fixed(d$effect, d$variance, d$set,comp1_var_max,threshold = threshold )
  p <- c(p$param$alpha[1], p$param$alpha[4], p$param$alpha[5]/(1-p$param$alpha[4]), p$param$mu[2], p$param$mu[3], p$param$var[2], p$param$var[3])
  p
  
}

## log likelihood function for hessian matrix
.ll_set <- function(parameters, lfc, var_gene, threshold, set){
  
  comp1_var_max <- seq(0,1, by = 0.00001)
  comp1_var_max <- comp1_var_max[which(pnorm(log2(threshold), 0, sqrt(comp1_var_max)) < 0.975)[1]]
  
  comp1 <- dnorm(lfc, 0, sqrt(comp1_var_max+var_gene)) * parameters[1]*set
  comp2 <- dnorm(lfc, max(parameters[4],qnorm(0.75, log2(threshold), sqrt(parameters[6]))) , sqrt(parameters[6]+var_gene)) * (1-parameters[1])*parameters[3]*set
  comp3 <- dnorm(lfc, min(parameters[5],qnorm(0.25, -log2(threshold), sqrt(parameters[7]))) , sqrt(parameters[7]+var_gene)) * (1-parameters[1])*(1-parameters[3])*set
  
  comp4 <- dnorm(lfc, 0, sqrt(comp1_var_max+var_gene)) * parameters[2]*(1-set)
  comp5 <- dnorm(lfc, max(parameters[4],qnorm(0.75, log2(threshold), sqrt(parameters[6]))) , sqrt(parameters[6]+var_gene)) * (1-parameters[2])*parameters[3]*(1-set)
  comp6 <- dnorm(lfc, min(parameters[5],qnorm(0.25, -log2(threshold), sqrt(parameters[7]))), sqrt(parameters[7]+var_gene)) * (1-parameters[2])*(1-parameters[3])*(1-set)
  
  l <- comp1 + comp2 + comp3 + comp4 + comp5 + comp6
  
  -sum(log(
    l
  ))
  
}


# fitting the ggm with all parameters (apart from nonDE component mean) free 
.EM_7FP_fixed <- function(effect, variance, effect_summary_df, comp1_var_max, threshold){
  n<- 1000
  m<- 1e-6
  iter<- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step <- .e_step_iter(effect, variance, c(0, 1, -1), c(0.005, 0.005, 0.005), c(.80, 0.1, 0.1))
      m.step <- .m_step_iter_fixed(effect,variance, c(0.005, 0.005, 0.005),iter, e.step[["posterior_df"]], comp1_var_max, threshold)
      cur.loglik <- e.step[["loglik"]]
      loglik.vector <- e.step[["loglik"]]
    } else {
      # Repeat E and M steps till convergence
      e.step <- .e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], m.step[["alpha"]])
      m.step <- .m_step_iter_fixed(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]], comp1_var_max, threshold)
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

## # fitting the ggm 1DF approach
.EM_1FP_fixed <- function(effect_set, variance_set, set, comp1_var_max, threshold){
  n<- 1000
  m<- 1e-6
  iter<- 10
  for (i in 1:n) {
    if (i == 1) {
      # Initialization
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, c(0, 1, -1), c(0.005, 0.005, 0.005), c(.80, 0.1, 0.1, 0.8, 0.1, 0.1))
      
      m.step.set <- .m_step_set_iter_fixed(effect_set, variance_set, set, c(0.005, 0.005, 0.005),iter, e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold)
      cur.loglik.set <- e.step.set[["loglik"]]
      loglik.vector.set <- e.step.set[["loglik"]]
    } else {
      
      # Repeat E and M steps till convergence
      e.step.set <- .e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]],m.step.set[["alpha"]])
      m.step.set <- .m_step_set_iter_fixed(effect_set,variance_set,set,m.step.set[["var"]],iter,  e.step.set[["posterior_df"]], comp1_var_max, threshold = threshold)
      
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
  parameters <- list("loglike" = loglike_set_genes, "param" = m.step.set, "ll.vector" = loglik.vector.set)
  return(parameters)
}




# the e_step and m_step functions in the EM algorithm
.e_step_iter <- function(x,lfc_var, mu_vector, component_var, alpha_vector) {
  # both lfc_se and component_var contribute to total variance 
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_var)) * alpha_vector[1]
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_var)) * alpha_vector[2]
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_var)) * alpha_vector[3]
  
  sum_of_comps <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps[which(sum_of_comps == 0)] <- 1e-200
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
.e_step_set_iter <- function(x,lfc_var,set, mu_vector, component_var, alpha_vector) {
  # both sd_vector and sigma are variance 
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_var)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_var)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_var)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[1], sqrt(component_var[1]+lfc_var)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(component_var[2]+lfc_var)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(component_var[3]+lfc_var)) * alpha_vector[6]*(1-set)
  
  sum_of_comps1 <- comp1_prod + comp2_prod + comp3_prod
  sum_of_comps2 <- comp4_prod + comp5_prod + comp6_prod
  sum_of_comps<- sum_of_comps1 + sum_of_comps2
  sum_of_comps[which(sum_of_comps == 0)] <- 1e-200
  
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





.m_step_iter_fixed <- function(x, lfc_var, component_var, t, posterior_df, comp1_var_max, threshold) {
  
  
  comp1_n <- sum(posterior_df[, 1]) 
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  
  
  ###########################
  for (i in 1:t) {     
    if (i == 1) { 
      # Initialization
      comp2_var <- component_var[2]
      ## weights
      w2_i <- 1/(comp2_var+lfc_var)
      
      # gets the mininum mean allowed when variance is comp2_var to keep 75% of the distribution above threshold
      comp2_mean_min <- qnorm(0.75, log2(threshold), sqrt(comp2_var))
      ## use either minimum value or mean estimate if bigger
      comp2_mu <- max(comp2_mean_min, sum(posterior_df[,2] * w2_i * x)/sum(posterior_df[,2] * w2_i))
      
      comp3_var <- component_var[3]
      ## weights
      w3_i <- 1/(comp3_var+lfc_var)
      # gets the maximum mean allowed when variance is comp3_var to keep 75% of the distribution below threshold
      comp3_mean_max <- qnorm(0.25, -log2(threshold), sqrt(comp3_var))
      ## use either maximum value or mean estimate if smaller
      comp3_mu <- min(comp3_mean_max,sum(posterior_df[,3] * w3_i * x)/sum(posterior_df[,3] * w3_i))
      comp1_var <- comp1_var_max
      
      comp1_mu<- 0
      
      
      
    } else {
      
      comp2_var <- max(0.005,sum(posterior_df[, 2] * ((w2_i)^2) * (((x - comp2_mu)^2)-lfc_var))/(sum(posterior_df[, 2]*w2_i^2)))
      w2_i <- 1/(comp2_var+lfc_var)
      comp2_mean_min <- qnorm(0.75, log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min,sum(posterior_df[, 2] * w2_i * x)/sum(posterior_df[, 2] * w2_i))
      
      comp3_var <- max(0.005,sum(posterior_df[, 3] * ((w3_i)^2) * (((x - comp3_mu)^2)-lfc_var))/(sum(posterior_df[, 3]*w3_i^2)))
      w3_i <- 1/(comp3_var+lfc_var)
      comp3_mean_max <- qnorm(0.25, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(posterior_df[, 3] * w3_i * x)/sum(posterior_df[, 3] * w3_i))
      
      comp1_var <-  comp1_var_max
      
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


.m_step_set_iter_fixed <- function(x, lfc_var, set, component_var, t, posterior_df, comp1_var_max, threshold) {
  
  
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
  
  
  
  # the following for loop iterates between mean and variance, the number of iteration is t, 
  comp1_pd <- posterior_df[,1] + posterior_df[,4]
  comp2_pd <- posterior_df[,2] + posterior_df[,5]
  comp3_pd <- posterior_df[,3] + posterior_df[,6]
  ###########################
  for (i in 1:t) {
    
    if (i == 1) {       
      comp2_var <- component_var[2]
      w2_i <- 1/(comp2_var+lfc_var)
      comp2_mean_min <- qnorm(0.75, log2(threshold), sqrt(comp2_var))
      comp2_mu <- max(comp2_mean_min, sum(comp2_pd * w2_i * x)/sum(comp2_pd * w2_i))
      
      comp3_var <- component_var[3]
      w3_i <- 1/(comp3_var+lfc_var)
      comp3_mean_max <- qnorm(0.25, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x)/sum(comp3_pd * w3_i))
      
      comp1_var <- comp1_var_max
      
      comp1_mu<- 0
      
      
      
    } else {
      comp2_var <- max(0.005,sum(comp2_pd * ((w2_i)^2) * (((x - comp2_mu)^2)-lfc_var))/(sum(comp2_pd*w2_i^2)))
      w2_i <- 1/(comp2_var+lfc_var)
      comp2_mean_min <- qnorm(0.75, log2(threshold), sqrt(comp2_var))
      comp2_mu <-max(comp2_mean_min, sum(comp2_pd * w2_i * x)/sum(comp2_pd * w2_i))
      
      comp3_var <- max(0.005,sum(comp3_pd * ((w3_i)^2) * (((x - comp3_mu)^2)-lfc_var))/(sum(comp3_pd*w3_i^2)))
      w3_i <- 1/(comp3_var+lfc_var)
      comp3_mean_max <- qnorm(0.25, -log2(threshold), sqrt(comp3_var))
      comp3_mu <- min(comp3_mean_max, sum(comp3_pd * w3_i * x)/sum(comp3_pd * w3_i))
      
      comp1_var <-  comp1_var_max
      comp1_mu<- 0
      
    }
  }  
  
  
  comp1_alpha <- max(comp1_alpha,lower_bound)
  comp2_alpha <- max(comp2_alpha,lower_bound)
  comp3_alpha <- max(comp3_alpha,lower_bound)
  comp4_alpha <- max(comp4_alpha,lower_bound)
  comp5_alpha <- max(comp5_alpha,lower_bound)
  comp6_alpha <- max(comp6_alpha,lower_bound)
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha ))
}


lower_bound<- 0.000000005

## plot the inset and the outset for a given set
## parameters returned by mrema() function
## x_lim: c(-x,x) specify your own x-axis
## y_lim: c(0,y) specify your own y-axis

plot_distributions <- function(parameters, x_lim = NULL, y_lim = NULL){
  inset <- parameters[c(1,2,3)]
  outset <- parameters[c(4,5,6)]
  ## get min and max values for the x-axis
  min <- min((inset$mu[3] - 2*inset$var[3]),(outset$mu[3] - 2*outset$var[3]))
  max <- max((inset$mu[2] + 2*inset$var[2]),(outset$mu[2] + 2*outset$var[2]))
  limit <- max(max, abs(min))
  if(is.null(x_lim) == TRUE) g <- seq(-limit,limit, by= .01) else g <- seq(x_lim[1], x_lim[2], by = 0.01)
  inset_value <- inset$alpha[1]*dnorm(g, inset$mu[1], sqrt(inset$var[1])) + inset$alpha[2]*dnorm(g, inset$mu[2], sqrt(inset$var[2])) + inset$alpha[3]*dnorm(g, inset$mu[3], sqrt(inset$var[3]))
  outset_value <- outset$alpha[1]*dnorm(g, outset$mu[1], sqrt(outset$var[1])) + outset$alpha[2]*dnorm(g, outset$mu[2], sqrt(outset$var[2])) + outset$alpha[3]*dnorm(g, outset$mu[3], sqrt(outset$var[3]))
  value <- c(inset_value, outset_value)
  gg <- c(g,g)
  group <- c(rep("Inset", length(g)), rep("Outset", length(g)))
  data <- dplyr::data_frame(Set = group, Density = value, LFC = gg)
  if(is.null(y_lim) == TRUE) y_limits = c(0,max(value)) else y_limits = y_lim
  plot <- ggplot2::ggplot(data = data, ggplot2::aes(x = LFC, y = Density, group = Set, color = Set))+
    ggplot2::geom_line() + ggplot2::ylim(y_limits) + ggplot2::xlab("Log Fold Change") +
    ggplot2::geom_area(ggplot2::aes( group = Set, fill = Set), alpha = 0.5, position = "identity") +
    #scale_fill_manual(values = c("#FC8D62", "#66C2A5")) +
    #scale_color_manual(values = c("#FC8D62", "#66C2A5")) +
    ggplot2::theme_minimal()
  
  return(plot)
}
