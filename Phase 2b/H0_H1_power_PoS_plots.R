# Define Functions
plot_H0_H1 <- function(mean_H0 = log(1),mean_H1 = log(0.73), N = 75, SD = 0.75, alpha = 0.05, filename = 'Y:\\CSLB\\phase 2b_3 power and decision boundary\\Outputs Plots and Reports\\H0_H1_plot.tiff'){
  s = SD*sqrt(2/N) # standard error of the mean
  thresh = s*qnorm(alpha) # threshold for statistical significance, one-sided alpha = 0.1
  tiff(filename, width = 8, height = 8, units = 'in', res = 100)
  x <- seq(min(mean_H0,mean_H1)-3*s,max(mean_H0,mean_H1)+3*s,(abs(mean_H0- mean_H1) + 6*s)/500)
  
  # H0
  dy1 <- dnorm(x,mean_H0,s)
  plot(x,dy1,type = 'l', xaxs="i", yaxs = 'i',ylim = range(0,1.2*max(dy1)), col = 'red',xlab = 'treatment effects', ylab = 'density', main = 'H0 vs H1 distributions, alpha and beta')
  segments(mean_H0,min(dy1),mean_H0,max(dy1),col='red',lty=3)
  
  # H1
  dy2 <- dnorm(x,mean_H1,s)
  lines(x,dy2,col='blue', type = 'l')
  segments(mean_H1,min(dy2),mean_H1,max(dy2),col='blue',lty=3)
  
  # threshold
  abline(v = thresh,lty=5)
  
  # alpha and beta
  polygon(c(x[x<=thresh], thresh),  c(dy1[x<=thresh],min(dy1) ), col="purple")
  polygon(c(x[x>=thresh], thresh),  c(dy2[x>=thresh],min(dy2) ), col="grey")
  alpha = round(pnorm(thresh, mean_H0,s),3)
  beta = round(pnorm(thresh, mean_H1, s, lower.tail = F),3)
  
  # legend
  legend('topleft', inset = 0.02, c('distribution assuming H0','distribution assuming H1'),col = c('red','blue'),lty = 1, cex=1.2)
  legend('topright', inset = 0.02, c(paste0('alpha = ',alpha),paste0('beta = ',beta)), col = c('purple','grey'), pch  = 15 ,cex = 1.2)
  dev.off()
}

plot_power_pos <- function(mean_H0 = log(1),mean_H1 = log(0.73), mu = -0.254, sigma = 0.152, N = 75, SD = 0.75, alpha = 0.1, filename = 'Y:\\CSLB\\clinical\\phase 2b_3 power and decision boundary\\Outputs Plots and Reports\\H0_H1_plot.tiff'){
  s = SD*sqrt(2/N) # standard error of the mean
  s2 = sqrt(s^2 + sigma^2) # standard error for posterior 
  thresh = s*qnorm(alpha) # threshold for statistical significance, one-sided alpha = 0.1
  ######################### power vs PoS plot
  tiff(filename, width = 8, height = 8, units = 'in', res = 100)
  x <- seq(min(mean_H1-3*s, mu - 3*s2),max(mean_H1+3*s, mu + 3*s2),(max(mean_H1+3*s, mu + 3*s2) - min(mean_H1-3*s, mu - 3*s2))/500)
  
  # around H1
  dy1 <- dnorm(x,mean_H1,s)
  plot(x,dy1,type = 'l', xaxs="i", yaxs = 'i',ylim = range(0,1.2*max(dy1)), col = 'blue',xlab = 'treatment effects', ylab = 'density', main = 'power vs probability of success')
  segments(mean_H1,min(dy1),mean_H1,max(dy1),col='blue',lty=3)
  
  
  # posterior
  dy2 <- dnorm(x,mu,s2)
  lines(x,dy2,type = 'l', col = 'green')
  segments(mu,min(dy2),mu,max(dy2),col='green',lty=3)
  
  # threshold
  abline(v = thresh,lty=5)
  
  # legend
  legend('topleft',inset = 0.02, c('distribution in power calculation','distribution in PoS calculation'),col = c('blue','green'),lty = 1,cex=1.2)
  dev.off()
  
}





############# H0 vs H1 plot
# mean_H0 = mean for Null hypothesis
# mean_H1 = mean for alternative hypothesis
# N = one-arm sample size
# SD = individual patient standard deviation
# alpha = alpha level
# beta = 1 - power
# filename = place and name for the picture to be saved
plot_H0_H1(mean_H0 = log(1),mean_H1 = log(0.73), N = 75, SD = 0.75, alpha = 0.05, filename = 'Y:\\CSLB\\clinical\\phase 2b_3 power and decision boundary\\Outputs Plots and Reports\\H0_H1_plot.tiff')





############## Power vs PoS plot
# mu = posterior mean of the mean for PoS calculation
# sigma = posterior standard deviation of the mean for PoS calculation
# N = one arm sample size
# SD = individual patient standard deviation
# alpha = alpha level
# filename = place and name for the picture to be saved
plot_power_pos(mean_H0 = log(1),mean_H1 = log(0.73), mu = -0.254, sigma = 0.152, N = 75, SD = 0.75, alpha = 0.1, filename = 'Y:\\CSLB\\clinical\\phase 2b_3 power and decision boundary\\Outputs Plots and Reports\\power_vs_pos_plot.tiff')











