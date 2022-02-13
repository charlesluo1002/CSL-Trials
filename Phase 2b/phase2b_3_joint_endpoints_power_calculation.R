###### Phase2b/3 boundary and power calculation using eGFR slope and log ACR

############# Read in tables ##########

# load library
library(mvtnorm)

# read in 10000 posterior samples with Omega and lambda2_gamma1 included.
posterior <-  read.csv('H:\\My Documents\\CSL work file\\further work\\posterior_samples_3endpoints.csv')
# read in SEs of the mean table for different study designs
SEs <- read.csv('H:\\My Documents\\CSL work file\\further work\\completeSlopeSEs.csv')

# change the file location to where the file is located in your computer

#######################################





############# Functions ###############
# function to extract the SEs for a specific study design
# inputs: 
# bgfr = base GFR
# YrChronic = year chronic
# Tint = time interval
# ngroup = number of patients per group
# gamma2_SD = assumed SD of log ACR.
SE_extract <- function(bgfr = 40, YrChronic = 1.5, Tint = 0.0833, ngroup = 200, gamma2_SD = 0.75){
  SE1 = SEs[SEs$bgfr==bgfr & SEs$YrChronic==YrChronic & SEs$Tint==Tint & SEs$ngroup==ngroup, ]$StdErr[1]
  SE2 = gamma2_SD*sqrt(2/ngroup)
  return(c(SE1,SE2))
}


# function to simulate a list of posterior thetas for one value pair of surrogate endpoints, used as part of make_table function, can also be used standalone.
# inputs:
# runs = number of thetas to simulate for each posterior sample
# posterior = posterior sample table
# gamma1_hat = estimated eGFR slope trt eff
# gamma2_hat = estimated log ACR
# SE = standard errors of the mean for the 2 surrogates.
theta_simulation <- function(runs = 100, posterior = posterior, gamma1_hat, gamma2_hat, SE = SE){
  SE1 = SE[1]
  SE2 = SE[2]
  thetas = rep(0, runs*nrow(posterior))
  for (j in 1:nrow(posterior)){
    ps = posterior[j,]
    sigma_prior = matrix(c(100*ps$omega^2 + ps$Lambda2_gamma1,100*ps$omega,100*ps$omega,100), 2)
    sigma_samp = matrix(c(SE1^2,-0.2*SE1*SE2,-0.2*SE1*SE2,SE2^2), 2)
    sigma_post = solve(solve(sigma_prior) + solve(sigma_samp))
    mu_post = sigma_post %*% solve(sigma_samp) %*% c(gamma1_hat, gamma2_hat) + sigma_post %*% solve(sigma_prior) %*% c(ps$mu2, ps$mu3)
    gammas0 = rmvnorm(1, mu_post, sigma_post)
    thetas[((j - 1)*runs + 1):(j*runs)] = rnorm(runs, ps$alphaClinonSur1Sur2 + ps$beta1ClinonSur1Sur2*gammas0[1] + ps$beta2ClinonSur1Sur2 * gammas0[2], sqrt(ps$sigSqClinonSur1Sur2))
  }
  return(thetas)
}


# function to produce a conf_lvl upper bound for HR on clinical outcomes.
# inputs:
# slope_trt_eff_values = vectors of estimated eGFR slope values
# acr_gmr_values = vectors of estimated ACR GMR values.
# conf_lvl = confidence level for HR. 0.9 refers to 90% upper bound HR. Can be a range of values, e.g. 0.9 (90%), 0.8, 0.7, 0.6, 0.5(median)
# runs, posterior, SE = parameters passed on to theta_simulation function
# Note: can take a few hours to run depending on how fine the grids one wants.
make_table <- function(slope_trt_eff_values=seq(0,1.2,0.02), acr_gmr_values=seq(1,0.6,-0.01), conf_lvl = c(0.9, 0.8, 0.7, 0.6, 0.5), runs = 100, posterior = posterior, SE=SE){
  gamma1_hat_values = sort(slope_trt_eff_values)
  acr_gmr_values = sort(acr_gmr_values,decreasing = TRUE)
  gamma2_hat_values = log(acr_gmr_values)
  table_list = replicate(length(conf_lvl),matrix(rep(0,length(gamma1_hat_values)*length(gamma2_hat_values)), length(gamma1_hat_values)),simplify = FALSE)
  table_list <- lapply(table_list, function(x) {dimnames(x) <- list(gamma1_hat_values, acr_gmr_values);return(x)})
  names(table_list) <- paste0('conf_lvl=',conf_lvl)

  for (i in 1:length(gamma1_hat_values)){
    for (j in 1:length(gamma2_hat_values)){
      simulated_theta = theta_simulation(runs = runs, posterior = posterior, gamma1_hat = gamma1_hat_values[i], gamma2_hat = gamma2_hat_values[j], SE=SE)
      for (k in 1:length(conf_lvl)){
        table_list[[k]][i,j] = exp(quantile(simulated_theta, conf_lvl[k]))
      }
    }
    cat(paste0('gamma1_hat = ', gamma1_hat_values[i],' done.\n'))
  }
  print('One setting done.\n')
  return(table_list)
}


# function to calculate the power for a list of diffeent assumed means.
# inputs:
# table_extract = one of the outputed table from make_table function
# thresh = threshold for defining the boundary, e.g. 1 for clinical benifit
# mu_gamma1_hat = vectors of assumed eGFR slope mean values
# mu_gamma2_hat = vectors of assumed log ACR values.
# approximation = methods for approximating the decision boundary, linear or quadratic
# runs = number of pairs of surrogate endpoints to simulate for calculating power
# SE = standard errors of the means, used for simulating gammas.
# slope_trt_eff_values, acr_gmr_values = vectors of estimated endpoint values as in make_table function. NEED TO MATCH WITH what's used in make_table function.

power_calc <- function(table_extract = tables[[1]], thresh = 1, mu_gamma1_hat = c(0,0.5,0.9), mu_gamma2_hat = log(c(1,0.85,0.73)),  approximation = 'quadratic', runs = 1000000, SE = SE, slope_trt_eff_values=seq(0,1.2,0.02), acr_gmr_values=seq(1,0.6,-0.01)){
  gamma1_vec = rep(0,ncol(table_extract))
  gamma1_seq = slope_trt_eff_values
  gamma2_vec = acr_gmr_values
  SE1 = SE[1]
  SE2 = SE[2]
  power_table <- matrix(rep(0,length(mu_gamma1_hat)*length(mu_gamma2_hat)),length(mu_gamma1_hat))
  dimnames(power_table) <- list(paste('slope = ', mu_gamma1_hat), paste0('acr_gmr = ',exp(mu_gamma2_hat)))
  # find boundary
  for (i in 1:length(gamma2_vec)){
    if (length(which(as.numeric(table_extract[,i])<thresh)) == 0 | length(which(as.numeric(table_extract[,i])<thresh)) == length(gamma1_seq)){
      gamma1_vec[i] = -999
      gamma2_vec[i] = -999
    }else{
      gamma1_vec[i] = gamma1_seq[min(which(as.numeric(table_extract[,i])<thresh))]
    }
  }
  gamma1_vec = gamma1_vec[gamma1_vec!=(-999)]
  gamma2_vec = gamma2_vec[gamma2_vec!=(-999)]
  for (j in 1:length(mu_gamma1_hat)){
    for (k in 1:length(mu_gamma2_hat)){
      if (approximation == 'linear'){
        #linear regression with gamma2
        lm1 <- lm(gamma1_vec~gamma2_vec)
        a = lm1$coefficients[1]
        b = lm1$coefficients[2]
        Mu = c(mu_gamma1_hat[j], mu_gamma2_hat[k])
        Sigma = matrix(c(SE1^2,-0.2*SE1*SE2,-0.2*SE1*SE2,SE2^2), 2)
        
        # simulate gammas
        gammas_hat <- rmvnorm(runs,Mu, Sigma)
        power_table[j,k] <- mean(gammas_hat[,1]>a+b*exp(gammas_hat[,2]))
      }
      if (approximation == 'quadratic'){
        #linear regression with second order term in gamma2
        gamma2_vec_sq <- gamma2_vec^2
        lm1 <- lm(gamma1_vec~gamma2_vec+gamma2_vec_sq)
        a = lm1$coefficients[1]
        b = lm1$coefficients[2]
        c = lm1$coefficients[3]
        Mu = c(mu_gamma1_hat[j], mu_gamma2_hat[k])
        Sigma = matrix(c(SE1^2,-0.2*SE1*SE2,-0.2*SE1*SE2,SE2^2), 2)
        
        #simulate gammas
        gammas_hat <- rmvnorm(runs,Mu, Sigma)
        power_table[j,k] <- mean(gammas_hat[,1]>a+b*exp(gammas_hat[,2])+c*exp(gammas_hat[,2])^2)
      }
    }
  }
  return(power_table)
}

# an all in one function that outputs a list of power tables. Each element of the list refers to one confidence level, and the table of powers is produced for different pre-defined assumed means. It is essentially a wrapper of all the previous functions.
# inputs:
# all inputs are passed on to other functions, detailed explanations can be found above.
all_in_one_power_calc <- function(conf_lvl = c(0.9, 0.8, 0.7,0.6,0.5), thresh = 1, mu_gamma1_hat = c(0, 0.5, 0.9), mu_gamma2_hat = log(c(1, 0.85, 0.73)), approximation = 'quadratic', slope_trt_eff_values=seq(0,1.2,0.02), acr_gmr_values=seq(1,0.6,-0.01), table_runs = 100, power_runs = 1000000, bgfr = 40, Yrchronic = 2.5, Tint = 0.0833, ngroup = 200, gamma2_SD = 0.75){
  SE <- SE_extract(bgfr, Yrchronic, Tint, ngroup, gamma2_SD)
  
  tables = make_table(slope_trt_eff_values=slope_trt_eff_values, acr_gmr_values=acr_gmr_values, conf_lvl = conf_lvl,runs = table_runs, posterior = posterior, SE = SE)
  
  power_table_list <- replicate(length(conf_lvl), matrix(rep(0,length(mu_gamma1_hat)*length(mu_gamma2_hat)),length(mu_gamma1_hat)),simplify = FALSE)
  power_table_list <- lapply(power_table_list, function(x) {dimnames(x) <- list(mu_gamma1_hat, exp(mu_gamma2_hat));return(x)})
  names(power_table_list) <- paste0('conf_lvl=',conf_lvl)
  
  for (i in 1:length(conf_lvl)){
    power_table_list[[i]] <- power_calc(table_extract = tables[[i]], thresh = thresh, mu_gamma1_hat = mu_gamma1_hat, mu_gamma2_hat = mu_gamma2_hat,  SE = SE, approximation = approximation, slope_trt_eff_values=slope_trt_eff_values, acr_gmr_values=acr_gmr_values, runs = power_runs)
  }
  return(power_table_list)
}

#######################################






##### Step by Step calculations #######

# 1. run all the previous codes above.

# 2. input study design settings, get SEs of mean.
SE = SE_extract(bgfr = 40, YrChronic = 1.5, Tint = 0.0833, ngroup = 200, gamma2_SD = 0.75)

# 3. produce a table of conf_lvl% upper bound table for HR of the clinical outcomes between patient groups.
set.seed(1)
tables = make_table( slope_trt_eff_values=seq(-1,2,0.05), acr_gmr_values=seq(1.25,0.05,-0.025), conf_lvl = c(0.9), runs = 100, posterior = posterior, SE=SE)

# the output contains a list of tables, each one is an upper bound table for one specific confidence level. In this case, tables[[1]] would refer to the 90% upper bound HR table.

# write them to file for easier access later.
# write.csv(tables[[1]], 'H:\\My Documents\\CSL work file\\further work\\Simulated tables Charles\\reproduced_90UppPredTable__fup1.5_N200.csv')

# 4. calculate powers for a list of assumed means for one specific confidence level upper bound table . The slope_trt_eff_values and acr_gmr_values need to match the ones used for producing the table.
# approximation = quadratic works better for larger scales.
write.csv(tables[[1]], paste0('H:\\My Documents\\CSL work file\\further work\\Simulated tables Charles\\final_90UpperPredTable_fup', YrTotal, '_N', ngroup, '.csv'))

power <- power_calc(table_extract = tables[[1]], thresh = 1, mu_gamma1_hat = seq(0,1.2,0.05), mu_gamma2_hat = log(seq(1,0.6,-0.01)),  SE = SE, approximation = 'quadratic', slope_trt_eff_values=seq(-1,2,0.05), acr_gmr_values=seq(1.25,0.05,-0.025) , runs = 1000000)

########################################






##### All-in-one power calculation #####

# all in one power calculation
# Would suggest running the 3 functions separately above since it gives intermediary outputs and allows for adjustments.
# gives a list of powers for different confidence levels. For example, powers[[1]] would give power table for design with 90% upper bound HR (alpha = 0.1) and various assumed means.
set.seed(1)
powers <- all_in_one_power_calc(conf_lvl = c(0.9,0.5), thresh = 1,mu_gamma1_hat = c(0,0.5,0.9), mu_gamma2_hat = log(c(1,0.85,0.73)), approximation = 'quadratic', slope_trt_eff_values=seq(-1,2,0.05), acr_gmr_values=seq(1.25,0.05,-0.025), table_runs = 100, power_runs = 1000000, bgfr = 40, Yrchronic = 1.5, Tint = 0.0833, ngroup = 200, gamma2_SD = 0.75)


########################################




