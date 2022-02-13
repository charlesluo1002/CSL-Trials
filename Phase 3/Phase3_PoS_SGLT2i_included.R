# phase 2b/3 probability of success table, based on outcomes of log ACR in phase 2a study for both the overall group and the SGLT2i subgroup.

library(mvtnorm)
posterior = read.csv('Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Data Background and Supplementary files\\PosteriorSample.csv')

########## functions
# PoS simulation function for one pair of ACR values
# Not to be used standalone, but as part of the PoS_table function.
pos_simulation <- function(gamma_hatall=log(0.75),gamma_hat1=log(0.85),t1=0.5,t2=0.7,posterior=posterior,SD1=0.75,event_total = 1004,k1=0.85,k2=0.85,alpha1=0.005, alpha2=0.025, epsilon = 0.1, N = 75, mu_prior = NULL, sigma2gamma = 100){
  pos_list = rep(0,nrow(posterior))
  gamma_hat2 = (gamma_hatall-t1*gamma_hat1) / (1-t1)
  delta_hat = gamma_hat1 - gamma_hat2
  T_matrix = matrix(c(1,1,(1-t1),-t1),2,byrow = F)
  S_matrix = matrix(c(t2,1,1-t2,0),2,byrow = F)
  epsilon2 = epsilon^2
  sigma_samp = matrix(c(SD1^2*2/N, 0, 0, SD1^2*2/N*(1/(t1*(1-t1)))), 2)
  sigma2_thetaallhat = 4/event_total
  sigma2_theta1hat = 4/(event_total*t2)
  sigma2_theta2hat = 4/(event_total*(1-t2))
  for (i in 1:nrow(posterior)){
    # select from posterior sample and define all the prior and posterior mean, covariance.
    if (is.null(mu_prior)){
      mu_prior = c(posterior$MuGamma[i],0)
    }
    alpha = posterior$intercept[i]
    beta = posterior$beta[i]
    lambda_2 = posterior$Lambda2[i]
    sigma_prior = matrix(c(sigma2gamma, (2*t1-1)*epsilon2/2, (2*t1-1)*epsilon2/2, epsilon2), 2)
    sigma_post = solve(solve(sigma_prior) + solve(sigma_samp))
    # simulate posterior gammas
    mu_post = sigma_post%*%solve(sigma_samp)%*%c(gamma_hatall, delta_hat) + sigma_post%*%solve(sigma_prior)%*%mu_prior
    mu_post_gamma = T_matrix%*%mu_post
    sigma_post_gamma = T_matrix%*%sigma_post%*%t(T_matrix)
    gammas = rmvnorm(1, mu_post_gamma, sigma_post_gamma)
    # define PoS thresholds for overall and SGLT2i subgroup
    u1 = min(sqrt(sigma2_thetaallhat)*qnorm(alpha1), log(k1))
    u2 = min(sqrt(sigma2_theta1hat)*qnorm(alpha2), log(k2))
    # define parameters for joint distribution of theta_hat_all and theta_hat_1
    Mu = c(alpha+beta*S_matrix%*%t(gammas))
    Sigma = S_matrix%*%matrix(c(lambda_2+sigma2_theta1hat,0,0,lambda_2+sigma2_theta2hat), 2)%*%t(S_matrix)
    # calculate and record PoS for this particular parameter sample
    pos_list[i] = pmvnorm(upper = c(u1,u2), mean=Mu, sigma=Sigma)
  }
  # return average of PoS for all 10000 samples
  return(mean(pos_list))
}

# pos_simulation(gamma_hatall=log(1),gamma_hat1=log(0.6),t1=0.5,t2=0.7,posterior=posterior,SD1=0.75,event_total = 1004,k1=0.85,k2=0.85,alpha1=0.005, alpha2=0.025, epsilon = 0.05, N = 75, mu_prior = NULL)


# function to produce 2 pos tables for a grid of ACR values. The first table has rows as non-SGLT2i ACR, columns as SGLT2i ACR. The second table has rows as overall ACR, columns as SGLT2i ACR.
PoS_table <- function(epsilon = 0.2, t1=0.5,t2=0.7,posterior=posterior,SD1=0.75,event_total = 844, N = 75, k1=0.8,k2=0.85,alpha1=0.005, alpha2=1,mu_prior = NULL, acrall_list = seq(1,0.5,-0.05), acr1_list = seq(1,0.5,-0.05), acr2_list =seq(1,0.5,-0.05)){
  # turn into model inputs gamma1_hat and gamma2_hat (sglt2i and non-sglt2i subgroup outcomes in phase 2a)
  gamma_hatall_list = log(acrall_list)
  gamma_hat1_list = log(acr1_list)
  gamma_hat2_list = log(acr2_list)
  
  # pre allocate pos table and names
  pos_table1 = data.frame(matrix(0,length(acrall_list),length(acr1_list)))
  pos_table2 = data.frame(matrix(0,length(acr2_list),length(acr1_list)))
  dimnames(pos_table1) = list(acrall_list,acr1_list)
  dimnames(pos_table2) = list(acr2_list,acr1_list)
  pos_table1 = data.frame(matrix(0,length(acrall_list),length(acr1_list)))
  pos_table2 = data.frame(matrix(0,length(acr2_list),length(acr1_list)))
  dimnames(pos_table1) = list(acrall_list,acr1_list)
  dimnames(pos_table2) = list(acr2_list,acr1_list)
  for (i in 1:length(gamma_hatall_list)){
    for (j in 1:length(gamma_hat1_list)){
      pos_table1[i,j]= pos_simulation(gamma_hatall=gamma_hatall_list[i],gamma_hat1=gamma_hat1_list[j],t1=t1,t2=t2,posterior=posterior,SD1=SD1,event_total = event_total,k1=k1,k2=k2,alpha1=alpha1, alpha2=alpha2, epsilon = epsilon, N = N, mu_prior = mu_prior)
    }
    cat(paste0('acrall_list = ',acrall_list[i],' done.\n'))
  }
  # run simulation of PoS for each value on grid with rows as non-SGLT2i acr, columns as SGLT2i acr
  for (i in 1:length(gamma_hat2_list)){
    for (j in 1:length(gamma_hat1_list)){
      pos_table2[i,j]= pos_simulation(gamma_hatall=(gamma_hat2_list[i] + gamma_hat1_list[j])/2,gamma_hat1=gamma_hat1_list[j],t1=t1,t2=t2,posterior=posterior,SD1=SD1,event_total = event_total,k1=k1,k2=k2,alpha1=alpha1, alpha2=alpha2, epsilon = epsilon, N = N, mu_prior = mu_prior)
    }
    cat(paste0('acr2_list = ',acr2_list[i],' done.\n'))
  }
  return(list(pos_table1,pos_table2))
}









########## produce tables for one setting

#### Produce tables (overall x SGLT2i   and    non-SGLT2i x SGLT2i) for epsilon = 0.2 and save to file
# epsilon = variability of delta
# t1 = % of SGLT2i events in phase 2a
# t2 = % of SGLT2i events in phase 3
# posterior = posterior sample from Tom
# SD1 = SD of log ACR in phase 2a 
# event_total = number of events in phase 3 
# N = number of patients in one arm in phase 2a
# k1 = hazard ratio cutoff in phase 3 for overall group
# k2 = hazard ratio cutoff in phase 3 for SGLT2i subgroup
# alpha1 = alpha level in phase 3 for overall group
# alpha2 = alpha level in phase 3 for SGLT2i subgrou
# mu_prior = prior mean log ACR for 2 subgroups (not needed)
# acrall_list = list of acrs to consider for overall group
# acr1_list = list of acrs to consider for SGLT2i subgroup
# acr2_list = list of acrs to consider for non-SGLT2 subgroup
epsilon = 1
t1 = 0.5
t2 = 0.7
# returns a list with 2 tables, first one non-SGLT2i x SGLT2i, second one overall x SGLT2i.
pos_tables = PoS_table(epsilon = epsilon, t1=t1,t2=t2,posterior=posterior,SD1=0.75,event_total = 844, N = 75, k1=0.8,k2=0.85,alpha1=0.005, alpha2=1,mu_prior = NULL, acrall_list = seq(1,0.5,-0.05), acr1_list = seq(1,0.5,-0.05), acr2_list =seq(1,0.5,-0.05))

# folder path to save the overall x SGLT2i PoS table
folder_path1 = 'Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\pos_table_2subgroup'
# folder path to save the non-SGLT2i x SGLT2i PoS table 
folder_path2 = 'Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\pos_table_overall_subgroup'


# save overall x SGLT2i table
write.csv(pos_tables[[1]], paste0(folder_path1, '\\pos_table_overall_subgroup_t1=',t1,'_t2=',t2,'_epsilon=',epsilon,'.csv'))
# save non-SGLT2i x SGLT2i table
write.csv(pos_tables[[2]], paste0(folder_path2, '\\pos_table_2subgroup_t1=',t1,'_t2=',t2,'_epsilon=',epsilon,'.csv'))






########## Produce a set of tables 
# produce PoS table for a list of epsilons and save them to a specified folder.
# inputs are the same as above
# epsilon_list = list of epsilon values to consider
epsilon_list = c(0.05,0.2,1,100)
t1 = 0.5
t2 = 0.7
folder_path1 = 'Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\pos_table_2subgroup'
folder_path2 = 'Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\pos_table_overall_subgroup'


for (i in 1:length(epsilon_list)){
  pos_tables = PoS_table(epsilon = epsilon_list[i], t1=0.5,t2=0.7,posterior=posterior,SD1=0.75,event_total = 844, N = 75, k1=0.8,k2=0.85,alpha1=0.005, alpha2=1,mu_prior = NULL, acrall_list = seq(1,0.5,-0.05), acr1_list = seq(1,0.5,-0.05), acr2_list =seq(1,0.5,-0.05))
  write.csv(pos_tables[[1]], paste0(folder_path1, '\\pos_table_overall_subgroup_t1=',t1,'_t2=',t2,'_epsilon=',epsilon,'.csv'))
  # non-SGLT2i x SGLT2i table
  write.csv(pos_tables[[2]], paste0(folder_path2, '\\pos_table_2subgroup_t1=',t1,'_t2=',t2,'_epsilon=',epsilon,'.csv'))
}







