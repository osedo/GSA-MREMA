source("packages.R")
lower_bound<- 0.000000005
# similar to sim.counts function from ssizeRNA package, not used in the current package,
# https://rdrr.io/cran/ssizeRNA/man/sim.counts.html
# but could be useful when we extend the current project, so you can ignore it now

sim.count.cell<- function (nGenes , pi0 , pi0_prime, pi0_dprime, m, mu, mu_prime, disp, fc, fc_prime,
                           up = 0.5, up_prime, cell_type, replace) 
{
  
  
  
  arg <- list(nGenes = nGenes, pi0 = pi0, pi0_prime = pi0_prime, pi0_dprime = pi0_dprime,
              group = rep(c(1, 2), each = m))
  
  FP <- round(nGenes * pi0)
  TP <- nGenes - FP
  DeS<- round(nGenes *  pi0_prime)
  DeSIn<- round(nGenes *  pi0_dprime)
  
  TP_up <- round(TP * up)
  TP_down <- TP - TP_up
  DeSIn_up<- round(DeSIn * up_prime)
  DeSIn_down <- DeSIn - DeSIn_up
  
  
  
  de <- c(rep(0, FP), rep(1, TP_up), rep(-1, TP_down))
  de <- de[sample.int(length(de))]
  
  
  # creating sl vector  
  t1<- which(de==1)
  tt1<- sample(t1, DeS*up)
  sl<- rep(0,nGenes)
  sl<- replace(sl, tt1, 1)
  t2<- which(de==-1)
  tt2<- sample(t2, DeS*(1-up))
  sl<- replace(sl, tt2, 1)
  #
  
  s1<- which(sl!=0 & de==1)
  ss1<- sample(s1, DeSIn*up)
  sss1<- sample(ss1, DeSIn_up*up)
  inter<- rep(0,nGenes)
  inter<- replace(inter, sss1, 1)
  
  
  `%notin%` <- Negate(`%in%`)
  
  r1<- ss1 %notin% sss1
  p1<- ss1[r1]
  inter<- replace(inter, p1, -1)
  
  
  s2<- which(sl!=0 & de==-1)
  ss2<- sample(s2, DeSIn*(1-up))
  sss2<- sample(ss2, DeSIn_up*(1-up))
  #inter<- rep(0,nGenes)
  inter<- replace(inter, sss2, 1)
  
  
  `%notin%` <- Negate(`%in%`)
  
  r2<- ss2 %notin% sss2
  p2<- ss2[r2]
  inter<- replace(inter, p2, -1)
  
  
  h <- rep(TRUE, nGenes)
  counts <- matrix(0, nrow = nGenes, ncol = 2 * m)
  delta <- rep(0, nGenes)
  slope <- rep(0, nGenes)
  interact <- rep(0, nGenes)
  
  if (is.function(fc)) {
    lfc <- log(fc(TP))
  }
  else {
    lfc <- log(fc)
    lfc_prime<-log(fc_prime) 
  }
  delta[de != 0] <- lfc * de[de != 0]
  slope[sl != 0] <- log(mu_prime) * sl[sl != 0]
  interact[inter != 0] <- lfc_prime * inter[inter != 0]
  
  selected_genes <- true_means <- true_disps <- rep(0, nGenes)
  left_genes <- 1:length(mu)
  lambda <- phi <- matrix(0, nrow = nGenes, ncol = 2 * m)
  while (any(h)) {
    temp <- sample.int(length(left_genes), sum(h), replace)
    temp <- temp[order(temp)]
    selected_genes[h] <- left_genes[temp]
    if (replace == FALSE) {
      left_genes <- left_genes[-temp]
    }
    true_means[h] <- log(mu[selected_genes[h]])
    true_disps[h] <- disp[selected_genes[h]]
    # lambda[h, ] <- matrix(true_means[h], ncol = 1) %*% matrix(rep(1, 
    #     2 * m), nrow = 1) + cbind(matrix(rep(0, sum(h) * 
    #                    m), ncol = m), matrix(rep(delta[h], m), ncol = m))
    
    lambda[h, ] <- matrix(true_means[h], ncol = 1) %*% matrix(rep(1, 
                                                                  2 * m), nrow = 1) +cbind(matrix(slope[h], ncol=1) %*% matrix(cell_type[1:m], nrow=1), 
                                                                                           (matrix(slope[h], ncol=1) %*% matrix(cell_type[(m+1):(2*m)], nrow=1)+
                                                                                              matrix(rep(delta[h], m), ncol = m)+matrix(interact[h], ncol=1) %*% matrix(cell_type[(m+1):(2*m)], nrow=1)) )                             
    
    
    
    
    
    phi[h, ] <- matrix(rep(true_disps[h], 2 * m), ncol = 2 * 
                         m)
    counts[h, ] <- rnegbin(sum(h) * 2 * m, exp(lambda[h, ]), 1/phi[h, 
                                                                   ])
    h <- (rowSums(cpm(counts) > 2) < 3)
  }
  if (any(rowSums(cpm(counts) > 2) < 3)) 
    print("Error: Failed to simulate data: some genes are not expressed.")
  delta <- delta/log(2)
  list(counts = counts, group = arg$group, lambda0 = lambda[, 
                                                            1], phi0 = phi[, 1], de = de, delta = delta, 
       sl=sl, inter = inter, slope=slope, interact=interact, cell_type=cell_type)
}




################3
################
################
# ignore this one too
# I am following the link :
# https://stat.ethz.ch/pipermail/r-help/2013-July/356936.html

scale1<- function (x, center = TRUE, scale = TRUE) 
{
  x <- as.matrix(x)
  nc <- ncol(x)
  if (is.logical(center)) {
    if (center) {
      center <- colMeans(x, na.rm = TRUE)
      x <- sweep(x, 2L, center, check.margin = FALSE)
    }
  }
  else if (is.numeric(center) && (length(center) == nc)) 
    x <- sweep(x, 2L, center, check.margin = FALSE)
  else stop("length of 'center' must equal the number of columns of 'x'")
  if (is.logical(scale)) {
    if (scale) {
      f <- function(v) {
        v <- v[!is.na(v)]
        sqrt(sum(v^2)/max(1, length(v) - 1L))
      }
      scale <- apply(x, 2L, f)
      x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
    }
  }
  else if (is.numeric(scale) && length(scale) == nc) 
    x <- sweep(x, 2L, scale, "/", check.margin = FALSE)
  else stop("length of 'scale' must equal the number of columns of 'x'")
  #if (is.numeric(center)) 
  #    attr(x, "scaled:center") <- center
  #if (is.numeric(scale)) 
  #    attr(x, "scaled:scale") <- scale
  x
}
#####################
###################
#################


# We want to generate n gene sets of equal size, the defined ?makeExampleData, generates 
# gene lists, but it cant generate lists of equal size, so i defined my own function
# see EnrichmentBrowser:::makeExampleData then EnrichmentBrowser:::.makeExmplGS
# hence i am writing my own function. Note the sets here might overlap
makeExampleGeneSets<- function(gnames, n, size){
  if (is.null(gnames)) 
    gnames <- paste0("g", 1:100)
  pre.gs<- replicate(n, sample(gnames, size))
  gs<- lapply(seq_len(ncol(pre.gs)), function(x) pre.gs[,x])
  names(gs) <- paste0("gs", seq_len(n))
  return(gs)
  
}
################
######################
#######################
# the above function generates gene sets but they might overlap, here we define a function 
# to generate non overlaping gene sets: here vector is a vector of gene names 
# like "    gnames <- paste0("g", 1:50)", x= the number of gene sets and y= the number of genes 
# in each gene set.
DisSample<- function(vector, x,y){
  sett<- list()
  for(i in 1:x){
    sett[[i]]<- sample(vector, y)
    vector<- vector[!vector %in% sett[[i]]]
    
  }
  return(sett)
}

# for example
# gnames <- paste0("g", 1:50)
# DisSample(gnames, 10,5)


#############
###########
# the e_step and m_step functions in the EM algorithm for w_o_mrema (weights only)
e_step_set_iter <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  # both sd_vector and sigma are variance 
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[6]*(1-set)
  
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

#### m-step

# m_step functions in the EM algorithm for w_o_mrema (weights only)
m_step_set_iter <- function(x, sigma,set,sd_vector, t, posterior_df) {
  

  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4]) 
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])
  
  
  # the following for loop iterates between mean and variance, the number of iteration is t, 
  
  ###########################
  for (i in 1:t) {     
    if (i == 1) {       
      # Initialization
      comp2_var <- sd_vector[2]
      comp2_mu <- (1/sum((posterior_df[, 2]+posterior_df[, 5]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 2]+posterior_df[, 5]) * x * 1/(comp2_var+sigma))
      
      comp3_var <- sd_vector[3]
      comp3_mu <- (1/sum((posterior_df[, 3]+posterior_df[, 6]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 3]+posterior_df[, 6]) * x * 1/(comp3_var+sigma))
      
      comp1_var <- min(sd_vector[1],0.05)
      comp1_mu<- 0
      
      
      
    } else {
      
      comp2_var <- max(0,sum((posterior_df[, 2]+posterior_df[, 5]) *((1/(comp2_var+sigma))^2)* (((x - comp2_mu)^2)-sigma))  * (1/sum((posterior_df[, 2]+posterior_df[, 5]) * 1/(comp2_var+sigma)^2)))
      
      comp2_mu <- (1/sum((posterior_df[, 2]+posterior_df[, 5]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 2]+posterior_df[, 5]) * x * 1/(comp2_var+sigma))
      
      comp3_var <- max(0,sum((posterior_df[, 3]+posterior_df[, 6]) *((1/(comp3_var+sigma))^2)* (((x - comp3_mu)^2)-sigma))  * (1/sum((posterior_df[, 3]+posterior_df[, 6]) * 1/(comp3_var+sigma)^2)))
      comp3_mu <- (1/sum((posterior_df[, 3]+posterior_df[, 6]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 3]+posterior_df[, 6]) * x * 1/(comp3_var+sigma))
      
      comp1_var <-  min(max(0,sum((posterior_df[, 1]+posterior_df[, 4]) *((1/(comp1_var+sigma))^2)* (((x - comp1_mu)^2)-sigma))  * (1/sum((posterior_df[, 1]+posterior_df[, 4]) * 1/(comp1_var+sigma)^2))),0.05)
      comp1_mu<- 0
      
      
      
      
    }
  }  
  
  
  
  
  
  
  comp1_alpha <- max(comp1_n / sum(set),lower_bound)
  comp2_alpha <- max(comp2_n / sum(set), lower_bound)
  comp3_alpha <- max(comp3_n / sum(set),lower_bound)
  comp4_alpha <- max(comp4_n / (length(set)-sum(set)),lower_bound)
  comp5_alpha <- max(comp5_n / (length(set)-sum(set)),lower_bound)
  comp6_alpha <- max(comp6_n / (length(set)-sum(set)),lower_bound)
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha ))
}


##########
#########
#############
###########
# the e_step and m_step functions in the EM algorithm for wm_o_mrema (weight and means only)
e_step_set_m_iter <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[4], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[5], sqrt(sd_vector[2]+sigma)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[6], sqrt(sd_vector[3]+sigma)) * alpha_vector[6]*(1-set)
  
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

#### m-step

# m_step functions in the EM algorithm for wm_o_mrema (weight and means only)
m_step_set_m_iter <- function(x, sigma,set,sd_vector, t, posterior_df) {
  
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  comp3_n <- sum(posterior_df[, 3])
  comp4_n <- sum(posterior_df[, 4])
  comp5_n <- sum(posterior_df[, 5])
  comp6_n <- sum(posterior_df[, 6])
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  ###########################
  for (i in 1:t) {     
    if (i == 1) {       
      # Initialization
      
      
      comp2_var <- sd_vector[2]
      comp2_mu <- (1/sum((posterior_df[, 2]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 2]) * x * 1/(comp2_var+sigma))
      comp5_mu <- (1/sum((posterior_df[, 5]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 5]) * x * 1/(comp2_var+sigma))

      
      comp3_var <- sd_vector[3]
      comp3_mu <-(1/sum((posterior_df[, 3]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 3]) * x * 1/(comp3_var+sigma))
      comp6_mu <-(1/sum((posterior_df[, 6]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 6]) * x * 1/(comp3_var+sigma))
      
      
      comp1_var <- min(sd_vector[1],0.05)
      comp1_mu<- 0
      comp4_mu <- 0
      
      
      
    } else {
      comp2_var <- max(0,(sum((posterior_df[, 2])* 1/(comp2_var+sigma)^2 * (((x - comp2_mu)^2)-sigma))+sum((posterior_df[, 5]) * 1/(comp2_var+sigma)^2* (((x - comp5_mu)^2)-sigma)))  *  (1/sum((posterior_df[, 2]+posterior_df[, 5]) * 1/(comp2_var+sigma)^2)))
      comp2_mu <- (1/sum((posterior_df[, 2]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 2]) * x * 1/(comp2_var+sigma))
      comp5_mu <- (1/sum((posterior_df[, 5]) * 1/(comp2_var+sigma))) * sum((posterior_df[, 5]) * x * 1/(comp2_var+sigma))
      
      
      comp3_var <- max(0,(sum((posterior_df[, 3])* 1/(comp3_var+sigma)^2 * (((x - comp3_mu)^2)-sigma))+sum((posterior_df[, 6]) * 1/(comp3_var+sigma)^2* (((x - comp6_mu)^2)-sigma)))  *  (1/sum((posterior_df[, 3]+posterior_df[, 6]) * 1/(comp3_var+sigma)^2)))
      comp3_mu <-(1/sum((posterior_df[, 3]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 3]) * x * 1/(comp3_var+sigma))
      comp6_mu <-(1/sum((posterior_df[, 6]) * 1/(comp3_var+sigma))) * sum((posterior_df[, 6]) * x * 1/(comp3_var+sigma))
      
      comp1_var <-  min(max(0,(sum((posterior_df[, 1])* 1/(comp1_var+sigma)^2 * (((x - comp1_mu)^2)-sigma))+sum((posterior_df[, 4]) * 1/(comp1_var+sigma)^2* (((x - comp4_mu)^2)-sigma)))  *  (1/sum((posterior_df[, 1]+posterior_df[, 4]) * 1/(comp1_var+sigma)^2))),0.05)
      comp1_mu<- 0
      comp4_mu <- 0
      
      
      
      
    }
  }  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  comp1_alpha <- max(comp1_n / sum(set),lower_bound)
  comp2_alpha <- max(comp2_n / sum(set),lower_bound)
  comp3_alpha <- max(comp3_n / sum(set),lower_bound)
  comp4_alpha <- max(comp4_n / (length(set)-sum(set)),lower_bound)
  comp5_alpha <- max(comp5_n / (length(set)-sum(set)),lower_bound)
  comp6_alpha <- max(comp6_n / (length(set)-sum(set)),lower_bound)
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu, comp6_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha, comp6_alpha ))
}


##########
#########
# plot
plot_mix_comps <- function(x, mu, sigma, lam) {
  lam * dnorm(x, mu, sigma)
}



##########
#########
# the e_step and m_step functions in the EM algorithm, mrema
e_step_iter <- function(x,sigma, mu_vector, sd_vector, alpha_vector) {
  # both sd_vector and sigma are variance 
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]
  
  sum_of_comps <- comp1_prod + comp2_prod + comp3_prod
  comp1_post <- comp1_prod / sum_of_comps
  comp2_post <- comp2_prod / sum_of_comps
  comp3_post <- comp3_prod / sum_of_comps
  
  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)
  
  list("loglik" = sum_of_comps_ln_sum,
       "posterior_df" = cbind(comp1_post, comp2_post, comp3_post))
}

#### m-step
m_step_iter <- function(x, sigma,sd_vector, t, posterior_df) {
  
  
  comp1_n <- sum(posterior_df[, 1]) 
  comp2_n <- sum(posterior_df[, 2])
    comp3_n <- sum(posterior_df[, 3])
  
  
  ###########################
  for (i in 1:t) {     
    if (i == 1) {       
      # Initialization
      comp2_var <- sd_vector[2]
      #comp2_mu <- (1/(max(sum(posterior_df[, 2] *   1/(comp2_var+sigma))     ,1e-308))) * sum(posterior_df[, 2] *(1/(comp2_var+sigma))* x)
      comp2_mu <- (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma))  )) * sum(posterior_df[, 2] *(1/(comp2_var+sigma))* x)
      
      comp3_var <- sd_vector[3]
      comp3_mu <- (1/(sum(posterior_df[, 3] *   1/(comp3_var+sigma))  )) * sum(posterior_df[, 3] *(1/(comp3_var+sigma))* x)
      
      comp1_var <- min(sd_vector[1],0.05)
      comp1_mu<- 0
      
      
      
    } else {
      comp2_var <-  max(0,sum(posterior_df[, 2] *((1/(comp2_var+sigma))^2)* (((x - comp2_mu)^2)-sigma))  *  (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma)^2)  ))       )
      
      comp2_mu <- (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma))  )) * sum(posterior_df[, 2] *(1/(comp2_var+sigma))* x)
      
      comp3_var <-  max(0,sum(posterior_df[, 3] *((1/(comp3_var+sigma))^2)* (((x - comp3_mu)^2)-sigma))  *  (1/(sum(posterior_df[, 3] *   1/(comp3_var+sigma)^2)  ))       )
      comp3_mu <- (1/(sum(posterior_df[, 3] *   1/(comp3_var+sigma))  )) * sum(posterior_df[, 3] *(1/(comp3_var+sigma))* x)
      
      comp1_var <-  min(max(0,sum(posterior_df[, 1] *((1/(comp1_var+sigma))^2)* (((x - comp1_mu)^2)-sigma))  *  (1/(sum(posterior_df[, 1] *   1/(comp1_var+sigma)^2)  ))       ),0.05)
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


##########
#########
# estep for the second self contained version
e_step_self_second_iter <- function(x,sigma, mu_vector, sd_vector, alpha_vector) {
  # both sd_vector and sigma are variance 
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]

  sum_of_comps <- comp1_prod + comp2_prod 
  comp1_post <- comp1_prod / sum_of_comps
  comp2_post <- comp2_prod / sum_of_comps

  sum_of_comps_ln <- log(sum_of_comps, base = exp(1))
  sum_of_comps_ln_sum <- sum(sum_of_comps_ln)
  
  list("loglik" = sum_of_comps_ln_sum,
       "posterior_df" = cbind(comp1_post, comp2_post))
}

#### m-step
m_step_self_second_iter <- function(x, sigma,sd_vector, t, posterior_df) {
  comp1_n <- sum(posterior_df[, 1])
  comp2_n <- sum(posterior_df[, 2])
  ###########################
  for (i in 1:t) {     
    if (i == 1) {       
      # Initialization
      comp2_var <- sd_vector[2]
      comp2_mu <- (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma))  )) * sum(posterior_df[, 2] *(1/(comp2_var+sigma))* x)
      
      
      comp1_var <- min(sd_vector[1],0.05)
      comp1_mu<- 0
      
      
      
    } else {
      comp2_var <-  max(0,sum(posterior_df[, 2] *((1/(comp2_var+sigma))^2)* (((x - comp2_mu)^2)-sigma))  *  (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma)^2)  ))       )
      
      comp2_mu <- (1/(sum(posterior_df[, 2] *   1/(comp2_var+sigma))  )) * sum(posterior_df[, 2] *(1/(comp2_var+sigma))* x)
      
      
      comp1_var <-  min(max(0,sum(posterior_df[, 1] *((1/(comp1_var+sigma))^2)* (((x - comp1_mu)^2)-sigma))  *  (1/(sum(posterior_df[, 1] *   1/(comp1_var+sigma)^2)  ))       ),0.05)
      comp1_mu<- 0
      
      
      
      
    }
  }  
  
  
  
  
  
  comp1_alpha <- max(comp1_n / length(x),lower_bound)
  comp2_alpha <- max(comp2_n / length(x),lower_bound)

  list("mu" = c(comp1_mu, comp2_mu),
       "var" = c(comp1_var, comp2_var),
       "alpha" = c(comp1_alpha, comp2_alpha))
}


##########
#########
##############################
################### Directionals ################################################


# the e_step and m_step functions in the EM algorithm for right directional (the negative part is fixed, genes are mostly up-regulated or not) w_o_mrema (weights only)
e_step_set_dir_right <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for right directional w_o_mrema (weights only)
m_step_set_dir_right <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <- 1/comp2_5_n * sum((posterior_df[, 2]+posterior_df[, 5]) * x)
  comp3_mu <- 1/comp3_6_n * sum((posterior_df[, 3]+posterior_df[, 6]) * x)
  
  comp1_var <- min(max(0,sum((posterior_df[, 1]+posterior_df[, 4]) * (((x - comp1_mu)^2)-sigma)) * 1/comp1_4_n),0.05)
  comp2_var <- max(0,sum((posterior_df[, 2]+posterior_df[, 5]) * (((x - comp2_mu)^2)-sigma))  * 1/comp2_5_n)
  comp3_var <- max(0,sum((posterior_df[, 3]+posterior_df[, 6]) * (((x - comp3_mu)^2)-sigma))  * 1/comp3_6_n)
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp3_n + n_i*comp6_n^{2} - n_o * n_i * comp6_n + n_o * comp3_n * comp6_n - n_o * n_i * 
    comp3_n - n_o^{2} * comp6_n + n_i * comp3_n * comp6_n + n_o * comp6_n^{2}  
  
  b2<- n_i * comp6_n + 3* n_o * comp3_n + n_o * comp6_n - comp6_n^{2} + n_i * comp3_n - comp3_n^{2}
  
  b3<- -2*comp3_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp3_n)*lambda1) / (n_o - comp6_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp3_n)*lambda2) / (n_o - comp6_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_n / mu)*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp5_n / lambda)*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_n / mu
  comp3_alpha <- comp3_6_n / (lambda+mu)
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp5_n / lambda
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}



########################################################



# the e_step and m_step functions in the EM algorithm for left directional (the positive part is fixed, genes are mostly down-regulated or not) w_o_mrema (weights only)
e_step_set_dir_left <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[5]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for left directional w_o_mrema (weights only)
m_step_set_dir_left <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <- 1/comp2_5_n * sum((posterior_df[, 2]+posterior_df[, 5]) * x)
  comp3_mu <- 1/comp3_6_n * sum((posterior_df[, 3]+posterior_df[, 6]) * x)
  
  comp1_var <- min(max(0,sum((posterior_df[, 1]+posterior_df[, 4]) * (((x - comp1_mu)^2)-sigma)) * 1/comp1_4_n),0.05)
  comp2_var <- max(0,sum((posterior_df[, 2]+posterior_df[, 5]) * (((x - comp2_mu)^2)-sigma))  * 1/comp2_5_n)
  comp3_var <- max(0,sum((posterior_df[, 3]+posterior_df[, 6]) * (((x - comp3_mu)^2)-sigma))  * 1/comp3_6_n)
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp2_n + n_i*comp5_n^{2} - n_o * n_i * comp5_n + n_o * comp2_n * comp5_n - n_o * n_i * 
    comp2_n - n_o^{2} * comp5_n + n_i * comp2_n * comp5_n + n_o * comp5_n^{2}  
  
  b2<- n_i * comp5_n + 3* n_o * comp2_n + n_o * comp5_n - comp5_n^{2} + n_i * comp5_n - comp5_n^{2}
  
  b3<- -2*comp2_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp2_n)*lambda1) / (n_o - comp5_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp2_n)*lambda2) / (n_o - comp5_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_n / mu)*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp6_n / lambda)*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_5_n / (lambda+mu)
  comp3_alpha <- comp3_n / mu
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp6_n / lambda
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}





########################################################



# the e_step and m_step functions in the EM algorithm for right directional (the negative part is fixed, genes are mostly up-regulated or not) for wm_o_mrema (weight and means only)
e_step_set_m_dir_right <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[4], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[5], sqrt(sd_vector[2]+sigma)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for right directional wm_o_mrema (weight and means only)
m_step_set_m_dir_right <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <- 1/comp2_n * sum((posterior_df[, 2]) * x)
  comp3_mu <- 1/comp3_6_n * sum((posterior_df[, 3]+posterior_df[, 6]) * x)
  
  comp4_mu <- 0
  comp5_mu <- 1/comp5_n * sum((posterior_df[, 5]) * x)

  
  
  
  comp1_var <- min(max(0,(sum((posterior_df[, 1]) * (((x - comp1_mu)^2)-sigma))+sum((posterior_df[, 4]) * (((x - comp4_mu)^2)-sigma))) * 1/comp1_4_n),0.05)
  comp2_var <- max(0,(sum((posterior_df[, 2]) * (((x - comp2_mu)^2)-sigma))+sum((posterior_df[, 5]) * (((x - comp5_mu)^2)-sigma)))  * 1/comp2_5_n)
  comp3_var <- max(0,(sum((posterior_df[, 3]) * (((x - comp3_mu)^2)-sigma))+sum((posterior_df[, 6]) * (((x - comp3_mu)^2)-sigma)))  * 1/comp3_6_n)
  
  
  
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp3_n + n_i*comp6_n^{2} - n_o * n_i * comp6_n + n_o * comp3_n * comp6_n - n_o * n_i * 
    comp3_n - n_o^{2} * comp6_n + n_i * comp3_n * comp6_n + n_o * comp6_n^{2}  
  
  b2<- n_i * comp6_n + 3* n_o * comp3_n + n_o * comp6_n - comp6_n^{2} + n_i * comp3_n - comp3_n^{2}
  
  b3<- -2*comp3_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp3_n)*lambda1) / (n_o - comp6_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp3_n)*lambda2) / (n_o - comp6_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu , sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_n / mu)*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp4_mu, sqrt(comp1_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp5_mu, sqrt(comp2_var+sigma)) * (comp5_n / lambda)*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_n / mu
  comp3_alpha <- comp3_6_n / (lambda+mu)
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp5_n / lambda
  
  
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}




########################################################



# the e_step and m_step functions in the EM algorithm for left directional (the positive part is fixed, genes are mostly down-regulated or not) for wm_o_mrema (weight and means only)
e_step_set_m_dir_left <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
 
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[4], sqrt(sd_vector[1]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[5], sqrt(sd_vector[3]+sigma)) * alpha_vector[5]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for left directional wm_o_mrema (weight and means only)
m_step_set_m_dir_left <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <-1/comp2_5_n * sum((posterior_df[, 2]+posterior_df[, 5]) * x)

  comp3_mu <- 1/comp3_n * sum((posterior_df[, 3]) * x)
  comp4_mu <- 0
  comp5_mu <- 1/comp6_n * sum((posterior_df[, 6]) * x)

  
  
  
  comp1_var <- min(max(0,(sum((posterior_df[, 1]) * (((x - comp1_mu)^2)-sigma))+sum((posterior_df[, 4]) * (((x - comp4_mu)^2)-sigma))) * 1/comp1_4_n),0.05)
  comp2_var <- max(0,(sum((posterior_df[, 2]) * (((x - comp2_mu)^2)-sigma))+sum((posterior_df[, 5]) * (((x - comp2_mu)^2)-sigma)))  * 1/comp2_5_n)
  comp3_var <- max(0,(sum((posterior_df[, 3]) * (((x - comp3_mu)^2)-sigma))+sum((posterior_df[, 6]) * (((x - comp5_mu)^2)-sigma)))  * 1/comp3_6_n)
  
  
  
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp2_n + n_i*comp5_n^{2} - n_o * n_i * comp5_n + n_o * comp2_n * comp5_n - n_o * n_i * 
    comp2_n - n_o^{2} * comp5_n + n_i * comp2_n * comp5_n + n_o * comp5_n^{2}  
  
  b2<- n_i * comp5_n + 3* n_o * comp2_n + n_o * comp5_n - comp5_n^{2} + n_i * comp2_n - comp2_n^{2}
  
  b3<- -2*comp2_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp2_n)*lambda1) / (n_o - comp5_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp2_n)*lambda2) / (n_o - comp5_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_n / mu)*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp4_mu, sqrt(comp1_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp5_mu, sqrt(comp3_var+sigma)) * (comp6_n / lambda)*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_5_n / (lambda+mu)
  comp3_alpha <- comp3_n / mu
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp6_n / lambda
  
  
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu),
       "var" = c(comp1_var, comp2_var, comp3_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}





########################################################



# the e_step and m_step functions in the EM algorithm for right directional (the negative part is fixed, genes are mostly up-regulated or not) for mrema 
e_step_set_all_dir_right <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[4], sqrt(sd_vector[4]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[5], sqrt(sd_vector[5]+sigma)) * alpha_vector[5]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for right directional mrema 
m_step_set_all_dir_right <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <- 1/comp2_n * sum((posterior_df[, 2]) * x)
  comp3_mu <- 1/comp3_6_n * sum((posterior_df[, 3]+posterior_df[, 6]) * x)
  
  comp4_mu <- 0
  comp5_mu <- 1/comp5_n * sum((posterior_df[, 5]) * x)

  
  
  
  comp1_var <- min(max(0,(sum((posterior_df[, 1]) * (((x - comp1_mu)^2)-sigma))) * 1/comp1_n),0.05)
  comp2_var <- max(0,(sum((posterior_df[, 2]) * (((x - comp2_mu)^2)-sigma)))  * 1/comp2_n)
  comp3_var <- max(0,(sum((posterior_df[, 3]) * (((x - comp3_mu)^2)-sigma))+sum((posterior_df[, 6]) * (((x - comp3_mu)^2)-sigma)))  * 1/comp3_6_n)
  
  
  
  comp4_var <- min(max(0,(sum((posterior_df[, 4]) * (((x - comp4_mu)^2)-sigma))) * 1/comp4_n),0.05)
  comp5_var <- max(0,(sum((posterior_df[, 5]) * (((x - comp5_mu)^2)-sigma)))  * 1/comp5_n)
  
  
  
  
  
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp3_n + n_i*comp6_n^{2} - n_o * n_i * comp6_n + n_o * comp3_n * comp6_n - n_o * n_i * 
    comp3_n - n_o^{2} * comp6_n + n_i * comp3_n * comp6_n + n_o * comp6_n^{2}  
  
  b2<- n_i * comp6_n + 3* n_o * comp3_n + n_o * comp6_n - comp6_n^{2} + n_i * comp3_n - comp3_n^{2}
  
  b3<- -2*comp3_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp3_n)*lambda1) / (n_o - comp6_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp3_n)*lambda2) / (n_o - comp6_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu , sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_n / mu)*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp4_mu, sqrt(comp4_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp5_mu, sqrt(comp5_var+sigma)) * (comp5_n / lambda)*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_6_n / (lambda+mu))*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_n / mu
  comp3_alpha <- comp3_6_n / (lambda+mu)
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp5_n / lambda
  
  
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu),
       "var" = c(comp1_var, comp2_var, comp3_var, comp4_var, comp5_var),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}




########################################################


# the e_step and m_step functions in the EM algorithm for left directional (the positive part is fixed, genes are mostly down-regulated or not) for mrema 
e_step_set_all_dir_left <- function(x,sigma,set, mu_vector, sd_vector, alpha_vector) {
  
  
  comp1_prod <- dnorm(x, mu_vector[1], sqrt(sd_vector[1]+sigma)) * alpha_vector[1]*set
  comp2_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*set
  comp3_prod <- dnorm(x, mu_vector[3], sqrt(sd_vector[3]+sigma)) * alpha_vector[3]*set
  # posteiror dist for genes outside of the set
  
  comp4_prod <- dnorm(x, mu_vector[4], sqrt(sd_vector[4]+sigma)) * alpha_vector[4]*(1-set)
  comp5_prod <- dnorm(x, mu_vector[2], sqrt(sd_vector[2]+sigma)) * alpha_vector[2]*(1-set)
  comp6_prod <- dnorm(x, mu_vector[5], sqrt(sd_vector[5]+sigma)) * alpha_vector[5]*(1-set)
  
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

#### m-step
# m_step functions in the EM algorithm for left directional mrema 
m_step_set_all_dir_left <- function(x, sigma,set, posterior_df) {
  

  
  
  comp1_n <- max(sum(posterior_df[, 1]),1e-308) 
  comp2_n <- max(sum(posterior_df[, 2]), 1e-308)
  comp3_n <- max(sum(posterior_df[, 3]), 1e-308)
  comp4_n <- max(sum(posterior_df[, 4]),1e-308) 
  comp5_n <- max(sum(posterior_df[, 5]), 1e-308)
  comp6_n <- max(sum(posterior_df[, 6]), 1e-308)
  
  comp1_4_n<- comp1_n + comp4_n
  comp2_5_n<- comp2_n + comp5_n
  comp3_6_n<- comp3_n + comp6_n
  
  comp1_mu <- 0
  comp2_mu <-1/comp2_5_n * sum((posterior_df[, 2]+posterior_df[, 5]) * x)

  comp3_mu <- 1/comp3_n * sum((posterior_df[, 3]) * x)
  comp4_mu <- 0
  comp5_mu <- 1/comp6_n * sum((posterior_df[, 6]) * x)

  
  
  
  
  
  
  comp1_var <- min(max(0,(sum((posterior_df[, 1]) * (((x - comp1_mu)^2)-sigma))) * 1/comp1_n),0.05)
  comp2_var <- max(0,(sum((posterior_df[, 2]) * (((x - comp2_mu)^2)-sigma))+sum((posterior_df[, 5]) * (((x - comp2_mu)^2)-sigma)))  * 1/comp2_5_n)
  comp3_var <- max(0,(sum((posterior_df[, 3]) * (((x - comp3_mu)^2)-sigma)))  * 1/comp3_n)
  
  
  
  
  
  comp4_var <- min(max(0,(sum((posterior_df[, 4]) * (((x - comp4_mu)^2)-sigma))) * 1/comp4_n),0.05)
  comp5_var <- max(0,(sum((posterior_df[, 6]) * (((x - comp5_mu)^2)-sigma)))  * 1/comp6_n)
  
  
  
  
  
  
  
  
  
  n_i<- sum(set)
  n_o<- (length(set)-sum(set))
  
  b1<- -n_o^{2}*comp2_n + n_i*comp5_n^{2} - n_o * n_i * comp5_n + n_o * comp2_n * comp5_n - n_o * n_i * 
    comp2_n - n_o^{2} * comp5_n + n_i * comp2_n * comp5_n + n_o * comp5_n^{2}  
  
  b2<- n_i * comp5_n + 3* n_o * comp2_n + n_o * comp5_n - comp5_n^{2} + n_i * comp2_n - comp2_n^{2}
  
  b3<- -2*comp2_n
  
  lambda1<- (-b2 + sqrt(b2^{2} - 4*b1*b3))/2*b3
  lambda2<- (-b2 - sqrt(b2^{2} - 4*b1*b3))/2*b3
  
  mu1<- (lambda1^{2}-(n_o + comp2_n)*lambda1) / (n_o - comp5_n - lambda1)
  mu2<- (lambda2^{2}-(n_o + comp2_n)*lambda2) / (n_o - comp5_n - lambda2)
  
  
  
  
  ex_likelihood<- function(lambda , mu){
    
    part1 <- (posterior_df[, 1])*log(dnorm(x, comp1_mu, sqrt(comp1_var+sigma)) * (comp1_n / mu)*set)
part2 <- (posterior_df[, 2])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*set)
part3 <- (posterior_df[, 3])*log(dnorm(x, comp3_mu, sqrt(comp3_var+sigma)) * (comp3_n / mu)*set)
# posteiror dist for genes outside of the set

part4 <- (posterior_df[, 4])*log(dnorm(x, comp4_mu, sqrt(comp4_var+sigma)) * (comp4_n / lambda)*(1-set))
part5 <- (posterior_df[, 5])*log(dnorm(x, comp2_mu, sqrt(comp2_var+sigma)) * (comp2_5_n / (lambda+mu))*(1-set))
part6 <- (posterior_df[, 6])*log(dnorm(x, comp5_mu, sqrt(comp5_var+sigma)) * (comp6_n / lambda)*(1-set))

pre_res<- part1+part2+part3+part4+part5+part6 
res<- sum(pre_res)  

  }
  
  like1<- ex_likelihood(lambda1 , mu1)
  like2<-  ex_likelihood(lambda2 , mu2) 
  
  if(like1 < like2){
    lambda<- lambda2
    mu<- mu2
  }else{
    lambda<- lambda1
    mu<- mu1 
  }
  
  
  
  comp1_alpha <- comp1_n / mu
  comp2_alpha <- comp2_5_n / (lambda+mu)
  comp3_alpha <- comp3_n / mu
  comp4_alpha <- comp4_n / lambda
  comp5_alpha <- comp6_n / lambda
  
  
  
  
  list("mu" = c(comp1_mu, comp2_mu, comp3_mu, comp4_mu, comp5_mu),
       "var" = c(comp1_var, comp2_var, comp3_var, comp4_var, comp5_var ),
       "alpha" = c(comp1_alpha, comp2_alpha, comp3_alpha, comp4_alpha, comp5_alpha ))
}






































