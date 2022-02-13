# Overall power calculation:
library(mvtnorm)
overallPower <- function(sem_sglt2 = 0.173, min_treat_eff_sglt2=-0.136, proportion_SGLT2_subject=0.5, min_treat_eff_prim = -0.158, mu=-0.315){
  t=proportion_SGLT2_subject
  Sigma = matrix(c(1,t,t,t),2) * (sem_sglt2^2)
  Mu = c(mu, mu)
  res = pmvnorm(upper = c(min_treat_eff_sglt2,min_treat_eff_prim), mean =Mu, sigma=Sigma)
  return(res)
}


## proportion = 0.5
# 1
overallPower(0.173, -0.136, 0.5)

# 2
overallPower(0.173, -0.093, 0.5)

# 3
overallPower(0.173, -0.030, 0.5)

## proportion = 0.6
# 4
overallPower(0.158, -0.152, 0.6)

# 5
overallPower(0.158, -0.113, 0.6)

# 6
overallPower(0.158, -0.055, 0.6)

## proportion = 0.7
# 7
overallPower(0.141, -0.170, 0.75)

# 8
overallPower(0.141, -0.135, 0.75)

# 9
overallPower(0.141, -0.083, 0.75)






