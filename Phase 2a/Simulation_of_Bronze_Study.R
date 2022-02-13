#### Simulation of Bronze study.
#### Goal: 4 conditional probabilities and assurance for different cutoff values on hazard ratio.

# Variables:
# 1. 10000 samples of joint posterior distribution of mu_theta, mu_gamma, sigma2gamma,beta,lambda2 and other derived parameters in excel.
# 2. Assumed mu (trt effect of ACR on log scale), sd, alpha,beta and dropout rate: mu_as, sd_as, alpha_as,power_as. Sample size S can be derived from these.
# 3. For simulation iterations, we sample l rows of the posterior samples, for each row, we simulate n set of true parameters gamma0, for each gamma0 simulated, we further simulate m instances of gamma0hat and theta0.
# 4. For the final estimands calculation, we would need the cutoff of hazard ratio for true clinical treatment effect: k.





## read in the posterior sample table
df_post = read.csv('H:/Desktop/R/data/PosteriorSample.csv')


set.seed(1)
# inputs:
# mu_as = assumed mean log ACR
# sd_as = assumed SD of log ACR for each patient
# power_as = power of study
# alpha_as = alpha level of study
# S = one arm sample size
# m = simulation runs
# n = another simulation runs
# l = number of samples from Posterior sample
# k_sequence = sequence of the cutoff values to be considered
quick_simulation <- function(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=100, n=1, l=1000, k_sequence = (500:1000)/1000){
  if (is.null(S)) S = ceiling(2*(qnorm(alpha_as)+qnorm(1-power_as))^2*sd_as^2/(mu_as)^2)
  sigma0 <- sd_as*sqrt(2/S)
  regionA_cutoff = qt(alpha_as, 2*S-2)*sigma0
  estimation_table = data.frame(matrix(0, length(k_sequence), 11), stringsAsFactors=F)
  colnames(estimation_table) <-  colnames(estimation_table) <- c('k','P(HR<=k | reject H0)', 'P(HR>k | not reject H0)', 'P(reject H0 | HR<=k)', 'P(not reject H0 | HR>k)', 't1o1', 't0o0','t1','t0','o1','o0')
  estimation_table$k = k_sequence
  # simulate data
  for (i in 1:l){
    sample_row <- as.numeric(df_post[i,])
    for (j in 1:n){
      gamma0 <- rnorm(1, sample_row[5], sqrt(sample_row[4]))
      gamma0hat <- rnorm(m,gamma0, sigma0)
      theta0 <- sample_row[11] + sample_row[6] * gamma0 + rnorm(m, 0, sqrt(sample_row[3]))
      for (kk in 1:length(k_sequence)){
        k = k_sequence[kk]
        estimation_table$t1o1[kk] = estimation_table$t1o1[kk] + sum(gamma0hat <= regionA_cutoff & theta0 <= log(k))
        estimation_table$t0o0[[kk]] = estimation_table$t0o0[kk] + sum(theta0 > log(k) & gamma0hat > regionA_cutoff)
        o1 = sum(gamma0hat <= regionA_cutoff)
        t1 =sum(theta0 <= log(k))
        estimation_table$o1[kk] = estimation_table$o1[kk] + o1
        estimation_table$o0[kk] = estimation_table$o0[kk] + m - o1
        estimation_table$t1[kk] = estimation_table$t1[kk] + t1
        estimation_table$t0[kk] = estimation_table$t0[kk] + m - t1
      }
    }
    if (i%%50 == 0) cat(i,'rounds done.\n')
  }
  
  estimation_table$`P(HR<=k | reject H0)` = with(estimation_table, t1o1/o1)
  estimation_table$`P(HR>k | not reject H0)` = with(estimation_table,t0o0/o0)
  estimation_table$`P(reject H0 | HR<=k)` = with(estimation_table,t1o1/t1)
  estimation_table$`P(not reject H0 | HR>k)` = with(estimation_table,t0o0/t0)
  return(estimation_table)
}





# estimate operating characteristics for different k.
quick_estimation_table <- quick_simulation(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=1000, n=1, l=10000, k_sequence = (500:1000)/1000)

# calculate assurance
assurance = quick_estimation_table$o1[1]/10000000 # Assurance = 0.6843737

# Assurance = 0.6843737


# quick_estimation_table[quick_estimation_table$k %in% c(0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1),]

# Plot the probabilities
library(ggplot2)
library(reshape2)
table_long <- melt(quick_estimation_table[,1:5], id='k')  # convert to long format

ggplot(data=table_long,
       aes(x=k, y=value, colour = variable)) +
  geom_line(size=1.5) +
  geom_vline(xintercept=c(0.7,0.75,0.8, 0.85), linetype="dotted") + 
  #geom_line(size=2, linetype = rep(1:4,each=length(k_sequence)))
  labs(title = 'probability estimates', x = 'Hazard ratio cutoff k', y = 'Probability', color = 'Legend Title\n') +
  scale_color_manual(values = c("blue", "green", 'red', 'yellow')) +
  theme(
    plot.title = element_text(size=12, face="bold.italic", hjust = 0.5),
    axis.title.x = element_text(size=10, face="bold"),
    axis.title.y = element_text(size=10, face="bold")
  )







## simulation function
# Inputs: i) posterior sample table, ii) mu_as, sd_as, alpha_as, power_as, or simply sample size of each group S, iii) sample and simulation iterations l, n, m iv) hazard ratio cutoff for true clinical trt effect: k
# Outputs: i) the 4 estimands required ii) the simulated table of gamma0hat, theta0 and gamma0 with l*m*n rows

#
# set.seed(1)
# simulation <- function(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=100, n=1, l=1000){
#   # create empty df with columns as gamma0hat, theta0, gamma0,  and the orignal columns of df_post
#   out_df = data.frame(matrix(0,m*n*l, ncol(df_post)+3, dimnames=list(c(), c('gamma0hat','theta0','gamma0',colnames(df_post)))), stringsAsFactors=F)
#   # projected sigma0 based on our sample
#   if (is.null(S)) S = ceiling(2*(qnorm(alpha_as)+qnorm(1-power_as))^2*sd_as^2/(mu_as)^2)
#   sigma0 <- sd_as*sqrt(2/S)
#   # simulate data
#   for (i in 1:l){
#     sample_row <- as.numeric(df_post[i,])
#     for (j in 1:n){
#       gamma0 <- rnorm(1, sample_row[5], sqrt(sample_row[4]))
#       for (k in 1:m){
#         gamma0hat <- rnorm(1,gamma0, sigma0)
#         theta0 <- sample_row[11] + sample_row[6] * gamma0 + rnorm(1, 0, sqrt(sample_row[3]))
#         out_df[(i-1)*n*m+(j-1)*m+k,] <- c(gamma0hat, theta0, gamma0, sample_row)
#       }
#     }
#     if (i%%50 == 0) cat(i,'rounds done.\n')
#   }
#   return(p_table = out_df)
# }
#
#
# simulation_of_1k <-  simulation(df_post = df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=10, n=1, l=100)
#
# write.csv(simulation_of_1k,'H:\\Desktop\\R\\data\\simulation_of_1k.csv')
#
#
# # given simualtion table, estimate all the estimands given one value of k
# estimation <- function(table, S=75, k = 0.8, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9){
#   if (is.null(S)) S = ceiling(2*(qnorm(alpha_as)+qnorm(1-power_as))^2*sd_as^2/(mu_as)^2)
#   sigma0 <- sd_as*sqrt(2/S)
#   regionA_cutoff = qt(alpha_as, 2*S-2)*sigma0
#   p_t1o1 = sum(table$theta0 <= log(k) & table$gamma0hat <= regionA_cutoff)/sum(table$gamma0hat <= regionA_cutoff)
#   p_t0o0 = sum(table$theta0 > log(k) & table$gamma0hat > regionA_cutoff)/sum(table$gamma0hat > regionA_cutoff)
#   p_o1t1 = sum(table$theta0 <= log(k) & table$gamma0hat <= regionA_cutoff)/sum(table$theta0 <= log(k))
#   p_o0t0 = sum(table$theta0 > log(k) & table$gamma0hat > regionA_cutoff)/sum(table$theta0 > log(k))
#   return(c(p_t1o1 = p_t1o1, p_t0o0 = p_t0o0, p_o1t1 = p_o1t1, p_o0t0 = p_o0t0))
# }
#
#
# # estimation given a sequence of cutoffs k
# k_sequence = (650:850)/1000
# estimation_table = data.frame(matrix(vector(), 0, 5), stringsAsFactors=F)
#
# for (k in k_sequence){
#   estimation_table[nrow(estimation_table)+1,] <- c(k,estimation(simulation_of_100k, S=75, k = k, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9))
# }
# colnames(estimation_table) <- c('k','P(HR<=k | reject H0)', 'P(HR>k | not reject H0)', 'P(reject H0 | HR<=k)', 'P(not reject H0 | HR>k)')
# print(estimation_table)
#
# # Key k values of interest: 0.8, 0.75
# prob_estimates_k0.8 <- format(estimation_table[estimation_table$k==0.8,],digits=3)
# prob_estimates_k0.75 <- format(estimation_table[estimation_table$k==0.75,],digits=3)
# prob_estimates_k0.7 <- format(estimation_table[estimation_table$k==0.7,],digits=3)
#




# Assurance2 <- function(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=100, n=1, l=10000, k_sequence = (500:1000)/1000, prior_mu_values = log(1+c(-0.3,-0.27,-0.24,-0.21)), prior_mu_distribution = rep(1/4,4)){
#   res = 0
#   l = length(prior_mu_values)
#   for (i in 1:l){
#     res = res + quick_simulation(df_post, mu_as=-0.315, sd_as=0.75, alpha_as=0.1, power_as=0.9, S=75, m=1000, n=1, l=10000, k_sequence = (500:1000)/1000)$o1[1]/10000000
#   }
#   return(res)
# }

