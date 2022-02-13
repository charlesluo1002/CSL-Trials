# extract data points from the simulation that satisfy 90% upper bound of HR < 1, in followup= 1.5 years and patients = 200/group setting.



# table_extract = table_90
# This table comes from the output in the main program.
table_extract = table_extended_90
thresh = 1
gamma1_vec = rep(0,ncol(table_extract))
gamma1_seq = seq(-10,30,1)/10
gamma2_vec = seq(2,0.05,-0.05)
#gamma1_seq = seq(0,1.2,0.02)
#gamma2_vec = seq(1,0.6,-0.01)

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



# ACR boundary

# quadratic decision boundary
plot(gamma1_vec~gamma2_vec,type = 'l',lwd =2 , xlab='ACR effect size', ylab = 'slope effect size', main = 'Quadratic decision boundary (R-squared = 0.999)')
#linear regression
x1 <- gamma2_vec
x2 <- gamma2_vec^2
lm1 <- lm(gamma1_vec~x1+x2)
summary(lm1)
z1 <- seq(min(x1),max(x1),(max(x1)-min(x1))/100)
y1 <- lm1$coefficients[1]+lm1$coefficients[2]*z1+lm1$coefficients[3]*(z1^2)
lines(z1,y1, col = 'blue')

legend("topleft", legend=c("Decision boundary", "Quadratic"),
       col=c("black", "blue"), lty = 1,lwd = c(2,1),cex=0.8)

# linear decision boundary
plot(gamma1_vec~gamma2_vec,type='l', lwd = 2, xlab='ACR effect size', ylab = 'slope effect size', main = 'Linear decision boundary (R-squared = 0.991)')
#linear regression
x1 <- gamma2_vec
lm1 <- lm(gamma1_vec~x1)
summary(lm1)
z1 <- seq(min(x1),max(x1),(max(x1)-min(x1))/100)
y1 <- lm1$coefficients[1]+lm1$coefficients[2]*z1
lines(z1,y1, col = 'blue')

legend("topleft", legend=c("Decision boundary", "Linear"),
       col=c("black", "blue"), lty = 1,lwd = c(2,1),cex=0.8)


# Theoretical line
# take mean of each sampled parameter, plot theoretical line using log(1-x) = -x
pm <- apply(posterior,2,mean)
sigma_prior = matrix(c(100*pm[18]^2 + pm[19],100*pm[18],100*pm[18],100), 2)
sigma_samp = matrix(c(SE1^2,-0.2*SE1*SE2,-0.2*SE1*SE2,SE2^2), 2)
sigma_post = solve(solve(sigma_prior) + solve(sigma_samp))
A1 = pm[c(15,16)] %*% sigma_post %*% solve(sigma_prior) %*% pm[c(3,4)]
A2 = sqrt(pm[17]+pm[c(15,16)]%*%sigma_post%*%pm[c(15,16)])*qnorm(0.9)
B1 = (pm[c(15,16)]%*%sigma_post%*%solve(sigma_samp))[1]
B2 = (pm[c(15,16)]%*%sigma_post%*%solve(sigma_samp))[2]
A = (-pm[14]-A1-A2)/B1
B = -B2/B1
abline(A,-B, col = 'red')




# log ACR boundary
gamma2_vec_log <- log(gamma2_vec)
plot(gamma1_vec~gamma2_vec_log, xlab = 'log ACR', ylab = 'slope effect size', main = 'simulated data points vs theoretical line')
#linear regression
lm2 <- lm(gamma1_vec~gamma2_vec_log)
summary(lm2)
abline(lm2$coefficients, col = 'blue')

# take mean of each sampled parameter, plot theoretical line
pm <- apply(posterior,2,mean)
sigma_prior = matrix(c(100*pm[18]^2 + pm[19],100*pm[18],100*pm[18],100), 2)
sigma_samp = matrix(c(SE1^2,-0.2*SE1*SE2,-0.2*SE1*SE2,SE2^2), 2)
sigma_post = solve(solve(sigma_prior) + solve(sigma_samp))
A1 = pm[c(15,16)] %*% sigma_post %*% solve(sigma_prior) %*% pm[c(3,4)]
A2 = sqrt(pm[17]+pm[c(15,16)]%*%sigma_post%*%pm[c(15,16)])*qnorm(0.9)
B1 = (pm[c(15,16)]%*%sigma_post%*%solve(sigma_samp))[1]
B2 = (pm[c(15,16)]%*%sigma_post%*%solve(sigma_samp))[2]
A = (-pm[14]-A1-A2)/B1
B = -B2/B1
abline(A,B, col='red')

abline(v=0)
abline(v=log(0.6))