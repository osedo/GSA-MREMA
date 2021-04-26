# we define three different models coded as three functions to perform gsa (and cell type specefic gsa).
# each of these functions takes a list of gene sets and a SummarizedExperiment object 
# (the code could be easily modified to take any expression matrix) and
# it returns the lists of gene sets along with the p-value corresponding to the LRT test and BIC value in each model 


#source("packages.R")

source("functions.R")



mrema<- function(sum_exp_raw_count, raw.gs){
 # the followings are the EM criterion
   n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  
  # the first step, which is fitting regression
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[4,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  # we use kmeans for initialized values, any other methods or values could be used  
  
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  
  
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
 
  
  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
    
  
  # we are done with fitting the mixture model to the whole data, in the following step we
  # fit the mixture model considering the set membership   
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
   
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    # I want to have the second cluster + and the third one negative -
    
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    
   
    
      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.inset <- e_step_iter(effect_inset, variance_inset, effect_summary_df_inset[["mu"]], effect_summary_df_inset[["variance"]],
                                      effect_summary_df_inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,effect_summary_df_inset[["variance"]],iter, e.step.inset[["posterior_df"]])
          cur.loglik.inset <- e.step.inset[["loglik"]]
          loglik.vector.inset <- e.step.inset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.inset <- e_step_iter(effect_inset, variance_inset, m.step.inset[["mu"]], m.step.inset[["var"]], 
                                      m.step.inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,m.step.inset[["var"]],iter, e.step.inset[["posterior_df"]])
          loglik.vector.inset <- c(loglik.vector.inset, e.step.inset[["loglik"]])
          
          loglik.diff.inset <- abs((cur.loglik.inset - e.step.inset[["loglik"]]))
          
          if(loglik.diff.inset < m) {
            break
          } else {
            cur.loglik.inset <- e.step.inset[["loglik"]]
          }
        }
      }  
      
      loglike_Inset_genes<-tail(loglik.vector.inset, n=1)
      
      
    
    
    
   
    #################
    # do the same proccess for genes outside of the gene set    
    
    set_specific_post_out<- set_specific_post %>% filter(set==0)
    
    effect_outset <- set_specific_post_out$GROUP_effect
    variance_outset<- set_specific_post_out$GROUP_variance 
    
    
      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.outset <- e_step_iter(effect_outset, variance_outset, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                       effect_summary_df[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,effect_summary_df[["variance"]],iter, e.step.outset[["posterior_df"]])
          cur.loglik.outset <- e.step.outset[["loglik"]]
          loglik.vector.outset <- e.step.outset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.outset <- e_step_iter(effect_outset, variance_outset, m.step.outset[["mu"]], m.step.outset[["var"]], 
                                       m.step.outset[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,m.step.outset[["var"]],iter, e.step.outset[["posterior_df"]])
          loglik.vector.outset <- c(loglik.vector.outset, e.step.outset[["loglik"]])
          
          loglik.diff.outset <- abs((cur.loglik.outset - e.step.outset[["loglik"]]))
          if(loglik.diff.outset < m) {
            break
          } else {
            cur.loglik.outset <- e.step.outset[["loglik"]]
          }
        }
      }  
      
      loglike_Outset_genes<-tail(loglik.vector.outset, n=1)
      
    
    
    
   
    
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 7)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 14*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
    
    
    
    
    v[j]<- aa$p.value[[1]]
    v_mstep[j]<- m.step.inset$alpha[1]<m.step.outset$alpha[1]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
    
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC )
  
  return(detected_gene_sets) 
  
  
  
  
  
  
  
  
  
}


#######################
#########################
w_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[4,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  
  # generate the initial mixing weights
  
  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
 
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
    
  
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  initial_alpha<- c( effect_summary_df[["alpha"]][1],effect_summary_df[["alpha"]][2],
                     effect_summary_df[["alpha"]][3],effect_summary_df[["alpha"]][1],
                     effect_summary_df[["alpha"]][2],effect_summary_df[["alpha"]][3])
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                        initial_alpha)
          m.step.set <- m_step_set_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                        m.step.set[["alpha"]])
          m.step.set <- m_step_set_iter(effect_set,variance_set,set,m.step.set[["var"]],iter,  e.step.set[["posterior_df"]])
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
      
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
   
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 2)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 9*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}




#######################
#########################
wm_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[4,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  
  # generate the initial mixing weights
  
  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    
    
    # I want to have the second cluster + and the third one negative -
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    
    initial_alpha<- c(effect_summary_df_inset[["alpha"]], effect_summary_df[["alpha"]])
    initial_mu<- c(effect_summary_df_inset[["mu"]], effect_summary_df[["mu"]])
    initial_var<- effect_summary_df[["variance"]]
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, initial_mu, initial_var,
                                          initial_alpha)
          m.step.set <- m_step_set_m_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                          m.step.set[["alpha"]])
          m.step.set <- m_step_set_m_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
      
    
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
   
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 4)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 11*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}






#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
#####################################################################################
# correct for purity: we add purity as a covariate to regress out purity 
# 1. estimate:
estimate_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr +cr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[5,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    # I want to have the second cluster + and the third one negative -
    
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.inset <- e_step_iter(effect_inset, variance_inset, effect_summary_df_inset[["mu"]], effect_summary_df_inset[["variance"]],
                                      effect_summary_df_inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset, effect_summary_df_inset[["variance"]],iter, e.step.inset[["posterior_df"]])
          cur.loglik.inset <- e.step.inset[["loglik"]]
          loglik.vector.inset <- e.step.inset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.inset <- e_step_iter(effect_inset, variance_inset, m.step.inset[["mu"]], m.step.inset[["var"]], 
                                      m.step.inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,m.step.inset[["var"]],iter, e.step.inset[["posterior_df"]])
          loglik.vector.inset <- c(loglik.vector.inset, e.step.inset[["loglik"]])
          
          loglik.diff.inset <- abs((cur.loglik.inset - e.step.inset[["loglik"]]))
          
          if(loglik.diff.inset < m) {
            break
          } else {
            cur.loglik.inset <- e.step.inset[["loglik"]]
          }
        }
      }  
      
      loglike_Inset_genes<-tail(loglik.vector.inset, n=1)
      
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    # do the same proccess for genes outside of the gene set    
    
    set_specific_post_out<- set_specific_post %>% filter(set==0)
    
    effect_outset <- set_specific_post_out$GROUP_effect
    variance_outset<- set_specific_post_out$GROUP_variance 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.outset <- e_step_iter(effect_outset, variance_outset, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                       effect_summary_df[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,effect_summary_df[["variance"]],iter, e.step.outset[["posterior_df"]])
          cur.loglik.outset <- e.step.outset[["loglik"]]
          loglik.vector.outset <- e.step.outset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.outset <- e_step_iter(effect_outset, variance_outset, m.step.outset[["mu"]], m.step.outset[["var"]], 
                                       m.step.outset[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,m.step.outset[["var"]],iter, e.step.outset[["posterior_df"]])
          loglik.vector.outset <- c(loglik.vector.outset, e.step.outset[["loglik"]])
          
          loglik.diff.outset <- abs((cur.loglik.outset - e.step.outset[["loglik"]]))
          if(loglik.diff.outset < m) {
            break
          } else {
            cur.loglik.outset <- e.step.outset[["loglik"]]
          }
        }
      }  
      
      loglike_Outset_genes<-tail(loglik.vector.outset, n=1)
      
    
    
    
    
    #outside_genes_plot<- data.frame(x = effect_outset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                            lam = m.step.inset$alpha[2]),
    #                colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                           lam = m.step.inset$alpha[3]),
    #              colour = "green", lwd = 1.5) +
    #  ylab("Density") +
    #  xlab("Values") +
    #  ggtitle("Outset genes GMM Fit")
    
    
    
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 7)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 14*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
    
    
    
    v[j]<- aa$p.value[[1]]
    v_mstep[j]<- m.step.inset$alpha[1]<m.step.outset$alpha[1]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
    
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC )
  
  return(detected_gene_sets) 
  
  
  
  
  
  
  
  
  
}
#######################
#########################
estimate_w_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr +cr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[5,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  initial_alpha<- c( effect_summary_df[["alpha"]][1],effect_summary_df[["alpha"]][2],
                     effect_summary_df[["alpha"]][3],effect_summary_df[["alpha"]][1],
                     effect_summary_df[["alpha"]][2],effect_summary_df[["alpha"]][3])
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                        initial_alpha)
          m.step.set <- m_step_set_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                        m.step.set[["alpha"]])
          m.step.set <- m_step_set_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 2)
    
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 9*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}




#######################
#########################
estimate_wm_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr +cr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[5,]^2)
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    
    
    # I want to have the second cluster + and the third one negative -
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    
    initial_alpha<- c(effect_summary_df_inset[["alpha"]], effect_summary_df[["alpha"]])
    initial_mu<- c(effect_summary_df_inset[["mu"]], effect_summary_df[["mu"]])
    initial_var<- effect_summary_df[["variance"]]
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, initial_mu, initial_var,
                                          initial_alpha)
          m.step.set <- m_step_set_m_iter(effect_set, variance_set, set,initial_var,iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                          m.step.set[["alpha"]])
          m.step.set <- m_step_set_m_iter(effect_set,variance_set,set, m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
   
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 4)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 11*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}




###############################
###############################
###############################
###############################
###############################
# the following functions are used for cell type specific gene set analysis
# the interaction term between purity and group membership is added

cancer_w_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10

  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  initial_alpha<- c( effect_summary_df[["alpha"]][1],effect_summary_df[["alpha"]][2],
                     effect_summary_df[["alpha"]][3],effect_summary_df[["alpha"]][1],
                     effect_summary_df[["alpha"]][2],effect_summary_df[["alpha"]][3])
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                        initial_alpha)
          m.step.set <- m_step_set_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                        m.step.set[["alpha"]])
          m.step.set <- m_step_set_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
   
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 2)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 9*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}


###############################
###############################
###############################
###############################
###############################

cancer_wm_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    
    
    # I want to have the second cluster + and the third one negative -
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    
    initial_alpha<- c(effect_summary_df_inset[["alpha"]], effect_summary_df[["alpha"]])
    initial_mu<- c(effect_summary_df_inset[["mu"]], effect_summary_df[["mu"]])
    initial_var<- effect_summary_df[["variance"]]
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, initial_mu, initial_var,
                                          initial_alpha)
          m.step.set <- m_step_set_m_iter(effect_set, variance_set, set,initial_var,iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                          m.step.set[["alpha"]])
          m.step.set <- m_step_set_m_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
      
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 4)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 11*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  #detected_gene_sets<- length(which(v<0.05))
  return(detected_gene_sets) 
  
  
  
  
  
}




###############################
###############################
###############################
###############################
###############################
cancer_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  
  cr<-colData(sum_exp_raw_count)$purity
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    # I want to have the second cluster + and the third one negative -
    
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    
    # 980
    
      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.inset <- e_step_iter(effect_inset, variance_inset, effect_summary_df_inset[["mu"]], effect_summary_df_inset[["variance"]],
                                      effect_summary_df_inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,effect_summary_df_inset[["variance"]],iter, e.step.inset[["posterior_df"]])
          cur.loglik.inset <- e.step.inset[["loglik"]]
          loglik.vector.inset <- e.step.inset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.inset <- e_step_iter(effect_inset, variance_inset, m.step.inset[["mu"]], m.step.inset[["var"]], 
                                      m.step.inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset, m.step.inset[["var"]],iter, e.step.inset[["posterior_df"]])
          loglik.vector.inset <- c(loglik.vector.inset, e.step.inset[["loglik"]])
          
          loglik.diff.inset <- abs((cur.loglik.inset - e.step.inset[["loglik"]]))
          
          if(loglik.diff.inset < m) {
            break
          } else {
            cur.loglik.inset <- e.step.inset[["loglik"]]
          }
        }
      }  
      
      loglike_Inset_genes<-tail(loglik.vector.inset, n=1)
      
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    # do the same proccess for genes outside of the gene set    
    
    set_specific_post_out<- set_specific_post %>% filter(set==0)
    
    effect_outset <- set_specific_post_out$GROUP_effect
    variance_outset<- set_specific_post_out$GROUP_variance 
    
      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.outset <- e_step_iter(effect_outset, variance_outset, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                       effect_summary_df[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,effect_summary_df[["variance"]],iter, e.step.outset[["posterior_df"]])
          cur.loglik.outset <- e.step.outset[["loglik"]]
          loglik.vector.outset <- e.step.outset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.outset <- e_step_iter(effect_outset, variance_outset, m.step.outset[["mu"]], m.step.outset[["var"]], 
                                       m.step.outset[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,m.step.outset[["var"]],iter, e.step.outset[["posterior_df"]])
          loglik.vector.outset <- c(loglik.vector.outset, e.step.outset[["loglik"]])
          
          loglik.diff.outset <- abs((cur.loglik.outset - e.step.outset[["loglik"]]))
          if(loglik.diff.outset < m) {
            break
          } else {
            cur.loglik.outset <- e.step.outset[["loglik"]]
          }
        }
      }  
      
      loglike_Outset_genes<-tail(loglik.vector.outset, n=1)
      
    
    
    
    #outside_genes_plot<- data.frame(x = effect_outset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                            lam = m.step.inset$alpha[2]),
    #                colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                           lam = m.step.inset$alpha[3]),
    #              colour = "green", lwd = 1.5) +
    #  ylab("Density") +
    #  xlab("Values") +
    #  ggtitle("Outset genes GMM Fit")
    
    
    
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 7)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 14*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
    
    
    v[j]<- aa$p.value[[1]]
    v_mstep[j]<- m.step.inset$alpha[1]<m.step.outset$alpha[1]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
    
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
  
  
  
  
}



###############################
###############################
###############################
###############################
###############################
# the only difference with the last three functions is that here purity is the proportion of cancer cells
normal_w_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  cr<- 1-cr
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
 
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance,m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  initial_alpha<- c( effect_summary_df[["alpha"]][1],effect_summary_df[["alpha"]][2],
                     effect_summary_df[["alpha"]][3],effect_summary_df[["alpha"]][1],
                     effect_summary_df[["alpha"]][2],effect_summary_df[["alpha"]][3])
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                        initial_alpha)
          m.step.set <- m_step_set_iter(effect_set, variance_set, set,effect_summary_df[["variance"]],iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                        m.step.set[["alpha"]])
          m.step.set <- m_step_set_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 2)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 9*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}



###############################
###############################
###############################
###############################
###############################

normal_wm_o_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  cr<-colData(sum_exp_raw_count)$purity
  cr<- 1-cr
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  ####################3
  #to check the distribution:
  #sum_exp_scaled_counts<- t(sum_exp_normalized_counts)
  #mat0<- sum_exp_scaled_counts[,c(1:4)]
  #mat1<- sum_exp_scaled_counts[,c(5:8)]
  #ex1<- apply(mat1,1,mean)
  #beta0<- apply(mat0,1,mean)
  #beta1<- ex1-beta0 
  #hist(beta1,breaks=10, freq = FALSE)
  #normal_dist <- fitdist(bbeta1, "norm")
  # dist of beta1, for genes inside the gene set:
  #p<- which(rownames(sum_exp_scaled_counts) %in% raw.gs[[1]] )
  #out<- sum_exp_scaled_counts[p,]
  #out0<- out[,c(1:4)]
  #out1<- out[,c(5:8)]
  #eex1<- apply(out1,1,mean)
  #bbeta0<- apply(out0,1,mean)
  #bbeta1<- eex1-beta0[p] 
  #hist(bbeta1,breaks=10, freq = FALSE)
  
  ###################3
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
  
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    
    effect_set <- set_specific_post$GROUP_effect
    variance_set<- set_specific_post$GROUP_variance 
    set<- set_specific_post$set
    
    
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    
    
    # I want to have the second cluster + and the third one negative -
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    
    initial_alpha<- c(effect_summary_df_inset[["alpha"]], effect_summary_df[["alpha"]])
    initial_mu<- c(effect_summary_df_inset[["mu"]], effect_summary_df[["mu"]])
    initial_var<- effect_summary_df[["variance"]]
    
    
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, initial_mu, initial_var,
                                          initial_alpha)
          m.step.set <- m_step_set_m_iter(effect_set, variance_set, set,initial_var,iter, e.step.set[["posterior_df"]])
          cur.loglik.set <- e.step.set[["loglik"]]
          loglik.vector.set <- e.step.set[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.set <- e_step_set_m_iter(effect_set, variance_set, set, m.step.set[["mu"]], m.step.set[["var"]], 
                                          m.step.set[["alpha"]])
          m.step.set <- m_step_set_m_iter(effect_set,variance_set,set,m.step.set[["var"]],iter, e.step.set[["posterior_df"]])
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
      
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    
    
    
    aa<-lr.test(-loglike_all_genes, -loglike_set_genes,alpha = 0.05, 4)
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 11*log(nrow(postdata))-2*(loglike_set_genes)
    
    
    v_mstep[j]<- m.step.set$alpha[1]<m.step.set$alpha[4]
    v[j]<- aa$p.value[[1]]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
}


###############################
###############################
###############################
###############################
###############################

normal_mrema<- function(sum_exp_raw_count, raw.gs){
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  
  cr<-colData(sum_exp_raw_count)$purity
  cr<- 1-cr
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr+cr+gr:cr))$coefficients)
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[6,]^2)
  
  
  
  
  effect <- postdata$GROUP_effect
  variance<- postdata$GROUP_variance 
  
  effect_kmeans <- kmeans(effect, 3)
  effect_kmeans_cluster <- effect_kmeans$cluster
  
  effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
  
  
  #cluster_plot<- effect_df %>%
  #  mutate(num = row_number()) %>%
  #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
  #  geom_point() +
  #  ylab("Values") +
  #  ylab("Data Point Number") +
  #  scale_color_discrete(name = "Cluster") +
  #  ggtitle("K-means Clustering")
  
  # generate the initial mixing means and sd:
  
  effect_summary_df <- effect_df %>%
    group_by(cluster) %>%
    summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
  # i want the first cluster be the one with mean 0
  effect_summary_df$abs<- abs(effect_summary_df$mu)
  effect_summary_df<- effect_summary_df %>% arrange(abs)
  effect_summary_df$mu[1]<- 0
  
  # I want to have the second cluster + and the third one negative -
  
  if (effect_summary_df$mu[2] < effect_summary_df$mu[3]) {
    m2<- effect_summary_df[2,]
    m3<- effect_summary_df[3,]
    effect_summary_df[2,]<- m3
    effect_summary_df[3,]<- m2
  }
  # quick look
  ### effect_summary_df %>%
  ###  select(cluster, mu, variance, std)
  

  effect_summary_df <- effect_summary_df %>%
    mutate(alpha = size / sum(size))
  
  effect_summary_df<- replace(effect_summary_df,is.na(effect_summary_df),0.005)
  
  # quick look:
  ### effect_summary_df %>%
  ###  select(cluster, size, alpha) 
 
    for (i in 1:n) {
      if (i == 1) {
        # Initialization
        e.step <- e_step_iter(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                              effect_summary_df[["alpha"]])
        m.step <- m_step_iter(effect,variance,effect_summary_df[["variance"]],iter, e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                              m.step[["alpha"]])
        m.step <- m_step_iter(effect,variance, m.step[["var"]],iter, e.step[["posterior_df"]])
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
    
  
  
  
  
  #all_genes_plot<- data.frame(x = effect) %>%
  #  ggplot() +
  #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
  #                 fill = "white") +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[1], sqrt(m.step$var[1]), 
  #                            lam = m.step$alpha[1]),
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[2], sqrt(m.step$var[2]), 
  #                            lam = m.step$alpha[2]),
  #               colour = "blue", lwd = 1.5) +
  #  stat_function(geom = "line", fun = plot_mix_comps,
  #                args = list(m.step$mu[3], sqrt(m.step$var[3]), 
  #                            lam = m.step$alpha[3]),
  #                colour = "green", lwd = 1.5) +
  #  ylab("Density") +
  #  xlab("Values") +
  # ggtitle("all genes GMM Fit")
  
  v<- c()
  v_mstep<- c()
  tol<- c()
  BIC<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    
    effect_inset <- set_specific_post_in$GROUP_effect
    variance_inset<- set_specific_post_in$GROUP_variance 
    
    effect_kmeans_inset <- kmeans(effect_inset, 3)
    effect_kmeans_cluster_inset <- effect_kmeans_inset$cluster
    
    effect_df_inset <- data_frame(x = effect_inset, cluster = effect_kmeans_cluster_inset)
    
    
    # cluster_plot_inset<- effect_df_inset %>%
    #    mutate(num = row_number()) %>%
    #    ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #    geom_point() +
    #    ylab("Values") +
    #    ylab("Data Point Number") +
    #    scale_color_discrete(name = "Cluster") +
    #    ggtitle("K-means Clustering Genes Inside Geneset")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_inset <- effect_df_inset %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    effect_summary_df_inset$abs<- abs(effect_summary_df_inset$mu)
    effect_summary_df_inset<- effect_summary_df_inset %>% arrange(abs)
    effect_summary_df_inset$mu[1]<- 0
    
    # I want to have the second cluster + and the third one negative -
    
    if (effect_summary_df_inset$mu[2] < effect_summary_df_inset$mu[3]) {
      m2<- effect_summary_df_inset[2,]
      m3<- effect_summary_df_inset[3,]
      effect_summary_df_inset[2,]<- m3
      effect_summary_df_inset[3,]<- m2
    }
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    

    effect_summary_df_inset <- effect_summary_df_inset %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_inset<- replace(effect_summary_df_inset,is.na(effect_summary_df_inset),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    

      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.inset <- e_step_iter(effect_inset, variance_inset, effect_summary_df_inset[["mu"]], effect_summary_df_inset[["variance"]],
                                      effect_summary_df_inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,effect_summary_df_inset[["variance"]],iter, e.step.inset[["posterior_df"]])
          cur.loglik.inset <- e.step.inset[["loglik"]]
          loglik.vector.inset <- e.step.inset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.inset <- e_step_iter(effect_inset, variance_inset, m.step.inset[["mu"]], m.step.inset[["var"]], 
                                      m.step.inset[["alpha"]])
          m.step.inset <- m_step_iter(effect_inset,variance_inset,m.step.inset[["var"]],iter, e.step.inset[["posterior_df"]])
          loglik.vector.inset <- c(loglik.vector.inset, e.step.inset[["loglik"]])
          
          loglik.diff.inset <- abs((cur.loglik.inset - e.step.inset[["loglik"]]))
          
          if(loglik.diff.inset < m) {
            break
          } else {
            cur.loglik.inset <- e.step.inset[["loglik"]]
          }
        }
      }  
      
      loglike_Inset_genes<-tail(loglik.vector.inset, n=1)
      
      
    
    
    
    
    # inside_genes_plot<- data.frame(x = effect_inset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                          lam = m.step.inset$alpha[2]),
    #              colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #               args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                          lam = m.step.inset$alpha[3]),
    #             colour = "green", lwd = 1.5) +
    #ylab("Density") +
    #xlab("Values") +
    #ggtitle("Inset genes GMM Fit")
    
    #################
    # do the same proccess for genes outside of the gene set    
    
    set_specific_post_out<- set_specific_post %>% filter(set==0)
    
    effect_outset <- set_specific_post_out$GROUP_effect
    variance_outset<- set_specific_post_out$GROUP_variance 
    
      for (i in 1:n) {
        if (i == 1) {
          # Initialization
          e.step.outset <- e_step_iter(effect_outset, variance_outset, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                                       effect_summary_df[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,effect_summary_df[["variance"]],iter, e.step.outset[["posterior_df"]])
          cur.loglik.outset <- e.step.outset[["loglik"]]
          loglik.vector.outset <- e.step.outset[["loglik"]]
        } else {
          # Repeat E and M steps till convergence
          e.step.outset <- e_step_iter(effect_outset, variance_outset, m.step.outset[["mu"]], m.step.outset[["var"]], 
                                       m.step.outset[["alpha"]])
          m.step.outset <- m_step_iter(effect_outset,variance_outset,m.step.outset[["var"]],iter, e.step.outset[["posterior_df"]])
          loglik.vector.outset <- c(loglik.vector.outset, e.step.outset[["loglik"]])
          
          loglik.diff.outset <- abs((cur.loglik.outset - e.step.outset[["loglik"]]))
          if(loglik.diff.outset < m) {
            break
          } else {
            cur.loglik.outset <- e.step.outset[["loglik"]]
          }
        }
      }  
      
      loglike_Outset_genes<-tail(loglik.vector.outset, n=1)
      
      
    
    
    
    
    #outside_genes_plot<- data.frame(x = effect_outset) %>%
    #  ggplot() +
    #  geom_histogram(aes(x, ..density..), binwidth = 1, colour = "black", 
    #                 fill = "white") +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[1], sqrt(m.step.inset$var[1]), 
    #                            lam = m.step.inset$alpha[1]),
    #                colour = "red", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[2], sqrt(m.step.inset$var[2]), 
    #                            lam = m.step.inset$alpha[2]),
    #                colour = "blue", lwd = 1.5) +
    #  stat_function(geom = "line", fun = plot_mix_comps,
    #                args = list(m.step.inset$mu[3], sqrt(m.step.inset$var[3]), 
    #                           lam = m.step.inset$alpha[3]),
    #              colour = "green", lwd = 1.5) +
    #  ylab("Density") +
    #  xlab("Values") +
    #  ggtitle("Outset genes GMM Fit")
    
    
    
    
    
    
    
    aa<-lr.test(-loglike_all_genes, -(loglike_Inset_genes+loglike_Outset_genes),alpha = 0.05, 7)
    
    
    BIC_all<- 7*log(nrow(postdata))-2*loglike_all_genes
    BIC_set<- 14*log(nrow(postdata))-2*(loglike_Inset_genes+loglike_Outset_genes)
    
    
    v[j]<- aa$p.value[[1]]
    v_mstep[j]<- m.step.inset$alpha[1]<m.step.outset$alpha[1]
    tol[j]<- length(raw.gs[[j]])
    BIC[j]<- BIC_set < BIC_all
  }
  
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v, "mid_dist"=v_mstep, "size"=tol, "BIC_value"=BIC)
  
  return(detected_gene_sets) 
  
  
  
  
  
  
  
  
  
}



###############################
###############################
###############################
###############################
###############################
# the self-contained test, just coded for standard gsa but could be easily modified to have the 
# self-contained cell type specific test

Self_mrema<- function(sum_exp_raw_count, raw.gs){
  # first of all data needs to be normalized using deseq2
  n<- 1000
  m<- 1e-6
  iter<- 10
  
  gr<- colData(sum_exp_raw_count)$GROUP
  # the first step, which is fitting regression
  reg1<- apply(sum_exp_normalized_counts,1, function(x) summary(lm(x~ gr))$coefficients)
  
  
  
  
  
  
  
  postdata<- tibble("genesEnsembl"=rownames(sum_exp_normalized_counts),"GROUP_effect"= reg1[2,], "GROUP_variance"=reg1[4,]^2)
  
  
  v<- c()
  for(j in 1:length(raw.gs)){
    
    set_specific_post<- postdata
    set_specific_post$set<- ifelse(postdata$genesEnsembl %in% raw.gs[[j]], 1, 0)
    set_specific_post_in<- set_specific_post %>% filter(set==1)
    postdataa<- set_specific_post_in
    effect <- postdataa$GROUP_effect
    variance<- postdataa$GROUP_variance 
    
    ################# first self contained function
    for (i in 1:t) {     
      if (i == 1) {       
        # Initialization
        
        comp1_var <- 0.05
        comp1_mu<- 0
        
        
        
      } else {
        
        
        comp1_var <-  min(max(0,sum(((1/(comp1_var+sigma))^2)* (((x - comp1_mu)^2)-sigma))  *  (1/(sum(1/(comp1_var+sigma)^2)  ))       ),0.05)
        comp1_mu<- 0
        
        
        
        
      }
    }  
    
    
    
    mu<- comp1_mu
    var<- comp1_var
    loglike_first<- sum(log(dnorm(effect,mu,sqrt(var+variance))))
    
    ################# second self contained function
    effect_kmeans_second <- kmeans(effect, 2)
    effect_kmeans_cluster_second  <- effect_kmeans_second$cluster
    
    effect_df_second <- data_frame(x = effect, cluster = effect_kmeans_cluster_second)
    
    
    #cluster_plot<- effect_df %>%
    #  mutate(num = row_number()) %>%
    #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #  geom_point() +
    #  ylab("Values") +
    #  ylab("Data Point Number") +
    #  scale_color_discrete(name = "Cluster") +
    #  ggtitle("K-means Clustering")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df_second  <- effect_df_second  %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    # i want the first cluster be the one with mean 0
    effect_summary_df_second$abs<- abs(effect_summary_df_second$mu)
    effect_summary_df_second <- effect_summary_df_second %>% arrange(abs)
    effect_summary_df_second$mu[1]<- 0
    
    
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    
    # generate the initial mixing weights
    
    effect_summary_df_second  <- effect_summary_df_second %>%
      mutate(alpha = size / sum(size))
    effect_summary_df_second <- replace(effect_summary_df_second,is.na(effect_summary_df_second),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    for (p in 1:n) {
      if (p == 1) {
        # Initialization
        e.step <- e_step_self_second_iter(effect, variance, effect_summary_df_second[["mu"]], effect_summary_df_second[["variance"]],
                                          effect_summary_df_second[["alpha"]])
        m.step <- m_step_self_second_iter(effect,variance, effect_summary_df_second[["variance"]],iter,e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step_self_second_iter(effect, variance, m.step[["mu"]], m.step[["var"]], 
                                          m.step[["alpha"]])
        m.step <- m_step_self_second_iter(effect,variance, m.step[["var"]],iter,e.step[["posterior_df"]])
        loglik.vector <- c(loglik.vector, e.step[["loglik"]])
        
        loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
        if(loglik.diff < m) {
          break
        } else {
          cur.loglik <- e.step[["loglik"]]
        }
      }
    }  
    
    loglike_second<-tail(loglik.vector, n=1)
    
    ################# third self contained function
    effect_kmeans <- kmeans(effect, 3)
    effect_kmeans_cluster <- effect_kmeans$cluster
    
    effect_df <- data_frame(x = effect, cluster = effect_kmeans_cluster)
    
    
    #cluster_plot<- effect_df %>%
    #  mutate(num = row_number()) %>%
    #  ggplot(aes(y = num, x = x, color = factor(cluster))) +
    #  geom_point() +
    #  ylab("Values") +
    #  ylab("Data Point Number") +
    #  scale_color_discrete(name = "Cluster") +
    #  ggtitle("K-means Clustering")
    
    # generate the initial mixing means and sd:
    
    effect_summary_df <- effect_df %>%
      group_by(cluster) %>%
      summarize(mu = mean(x), variance = var(x), std = sd(x), size = n())
    # i want the first cluster be the one with mean 0
    effect_summary_df$abs<- abs(effect_summary_df$mu)
    effect_summary_df<- effect_summary_df %>% arrange(abs)
    effect_summary_df$mu[1]<- 0
    
    
    # quick look
    ### effect_summary_df %>%
    ###  select(cluster, mu, variance, std)
    
    # generate the initial mixing weights
    
    effect_summary_df <- effect_summary_df %>%
      mutate(alpha = size / sum(size))
    effect_summary_df <- replace(effect_summary_df,is.na(effect_summary_df),0.005)
    
    # quick look:
    ### effect_summary_df %>%
    ###  select(cluster, size, alpha) 
    for (t in 1:n) {
      if (t == 1) {
        # Initialization
        e.step <- e_step(effect, variance, effect_summary_df[["mu"]], effect_summary_df[["variance"]],
                         effect_summary_df[["alpha"]])
        m.step <- m_step(effect,variance, effect_summary_df[["variance"]],iter,e.step[["posterior_df"]])
        cur.loglik <- e.step[["loglik"]]
        loglik.vector <- e.step[["loglik"]]
      } else {
        # Repeat E and M steps till convergence
        e.step <- e_step(effect, variance, m.step[["mu"]], m.step[["var"]], 
                         m.step[["alpha"]])
        m.step <- m_step(effect,variance, m.step[["var"]],iter,e.step[["posterior_df"]])
        loglik.vector <- c(loglik.vector, e.step[["loglik"]])
        
        loglik.diff <- abs((cur.loglik - e.step[["loglik"]]))
        if(loglik.diff < m) {
          break
        } else {
          cur.loglik <- e.step[["loglik"]]
        }
      }
    }  
    
    loglike_third<-tail(loglik.vector, n=1)
    
    
    BIC_first<- 1*log(nrow(postdataa))-2*loglike_first
    BIC_second<- 4*log(nrow(postdataa))-2*loglike_second
    BIC_third<- 7*log(nrow(postdataa))-2*loglike_third
    
    
    
    
    
    
    
    bic<- min(c(BIC_second , BIC_third))
    if(BIC_first < bic){
      aa<- "not enriched"
    }else{
      aa<- "enriched"
    }
    
    
    
    
    
    
    v[j]<- aa
  }
  detected_gene_sets<- tibble("GeneSets"=names(raw.gs),"p-values"=v)
  
  return(detected_gene_sets)
}

















