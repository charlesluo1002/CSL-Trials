
########### log-normal plots for epsilon elicitation

# generate pdf of log normal for different epsilons
# create tiff file
tiff("Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\R_plot_for_different_epsilons.tiff", width = 8, height = 8, units = 'in', res = 600)

x = seq(0,3,0.001)
y1 = dlnorm(x, meanlog = 0, sdlog = 0.1)
y2 = dlnorm(x, meanlog = 0, sdlog = 0.2)
y3 = dlnorm(x, meanlog = 0, sdlog = 0.5)
y4 = dlnorm(x, meanlog = 0, sdlog = 2)

# plot distributions
plot(x,y1,type=  'l', main = 'Distribution of the ratio R for different epsilons', xlab = 'R', ylab = 'density', col = 'red', lwd = 3)
lines(x,y2, type = 'l', col = 'blue', lwd = 3)
lines(x,y3, type = 'l', col ='green', lwd = 3)
lines(x,y4, type = 'l', col ='orange', lwd = 3)
segments(1,0,1,100, lty = 3, col = 'grey',lwd=2)
legend('topright', c('epsilon = 0.1','epsilon = 0.2','epsilon = 0.5','epsilon = 2'),lty = 1, col = c('red','blue','green', 'orange'), lwd=3)
dev.off()




# generate pdf for one epsilon
# create tiff file
tiff("Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\R_plot_for_epsilon=0.2.tiff", width = 8, height = 8, units = 'in', res = 600)

x = seq(0,3,0.001)
y1 = dlnorm(x, meanlog = 0, sdlog = 0.2)
plot(x, y1, type=  'l', main = 'Distribution of the ratio R when epsilon = 0.2', xlab = 'R', ylab = 'density', col = 'blue', lwd = 3)
axis(1, at = c(0.73, 1/0.73), labels = c('',''))
segments(1,0,1,100, lty = 3, col = 'grey',lwd=2)

# fill tail areas
polygon(c(x[x>=1/0.73], 1/0.73),  c(y1[x>=1/0.73],min(y1) ), col="green")
polygon(c(x[x<=0.73], 0.73),  c(y1[x<=0.73],min(y1) ), col="green")

# add comments
text(1.7,0.25,'~5%')
text(0.4,0.25,'~5%')
text(1,0.25, '~90%')
# text(1.5,0.7, 'R=1/0.73')
# text(0.55,0.7, 'R=0.73')
dev.off()





########### normal plots for epsilon elicitation

# generate pdf of log normal for different epsilons
# create tiff file
tiff("Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\R_plot_for_different_epsilons_normal.tiff", width = 8, height = 8, units = 'in', res = 600)

x = seq(-2,2,0.001)
y1 = dnorm(x, mean = 0, sd = 0.1)
y2 = dnorm(x, mean = 0, sd = 0.2)
y3 = dnorm(x, mean = 0, sd = 0.5)
y4 = dnorm(x, mean = 0, sd = 1)

# plot distributions
plot(x,y1,type=  'l', xaxt='n', ylim = range(0.15,4.1),main = 'Distribution of the log(R) for different epsilons', xlab = 'R', ylab = 'density for log(R)', col = 'red', lwd = 3)
axis(1, at = c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2), labels = round(exp(c(-2,-1.5,-1,-0.5,0,0.5,1,1.5,2)),2))
lines(x,y2, type = 'l', col = 'blue', lwd = 3)
lines(x,y3, type = 'l', col ='green', lwd = 3)
lines(x,y4, type = 'l', col ='orange', lwd = 3)
segments(0,0,0,100, lty = 3, col = 'grey',lwd=2)
legend('topright', c('epsilon = 0.1','epsilon = 0.2','epsilon = 0.5','epsilon = 1'),lty = 1, col = c('red','blue','green', 'orange'), lwd=3)
dev.off()




# generate pdf for one epsilon
# create tiff file
tiff("Y:\\CSLB\\clinical\\phase 3 SGLT2i PoS\\Outputs Plots and Reports\\R_plot_for_epsilon=0.2_normal.tiff", width = 8, height = 8, units = 'in', res = 600)

x = seq(-0.6,0.6,0.001)
y1 = dnorm(x, mean= 0, sd = 0.2)
plot(x, y1, type=  'l', xaxt='n',ylim = range(0.1,2), main = 'Distribution of the log(R) when epsilon = 0.2', xlab = 'R', ylab = 'density of log(R)', col = 'blue', lwd = 3)
axis(1, at = log(c(0.73,1,1/0.73)), labels = c('0.73', '1','1/0.73'))
segments(0,0,0,100, lty = 3, col = 'grey',lwd=2)
legend('topright','region of extreme differences',col='green',pch = 15,cex=  1.2)
# fill tail areas
polygon(c(x[x>=log(1/0.73)], log(1/0.73)),  c(y1[x>=log(1/0.73)],min(y1) ), col="green")
polygon(c(x[x<=log(0.73)], log(0.73)),  c(y1[x<=log(0.73)],min(y1) ), col="green")

# add comments
text(log(1.7),0.25,'~5%',cex = 1.5)
text(log(0.6),0.25,'~5%',cex = 1.5)
text(0,0.25, '~90%',cex = 1.5)
# text(1.5,0.7, 'R=1/0.73')
# text(0.55,0.7, 'R=0.73')
dev.off()





##### Shrinkage of 2 subgroup treatment effects demonstration

x <- c(0.7, 0.9)
plot(0.7, xlim = c(0.7, 0.9),ylim = range(0.01,0.3), axes=FALSE, type = "n", xlab = "non-SGLT2i ACR GMR                                                SGLT2i ACR GMR", main = "Observed phase 2a outcome")
axis(1, at = x, labels = x)
segments(x[1],0,x[1],0.3, col = 'red',lwd = 3)
segments(x[2],0,x[2],0.3, col = 'blue',lwd = 3)


x <- c(0.72, 0.88)
plot(0, xlim = c(0.7, 0.9),ylim = range(0.01,0.3), axes=FALSE, type = "n", xlab = "non-SGLT2i ACR GMR                                                SGLT2i ACR GMR", main = "Posterior ACR means with mild similarity assumption")
axis(1, at = c(0.7,x, 0.9), labels = c(0.7,x,0.9))
segments(x[1],0,x[1],0.3, col = 'red',lwd = 3)
segments(x[2],0,x[2],0.3, col = 'blue',lwd = 3)


x <- c(0.79, 0.81)
plot(0, xlim = c(0.7, 0.9),ylim = range(0.01,0.3), axes=FALSE, type = "n", xlab = "non-SGLT2i ACR GMR                                                SGLT2i ACR GMR", main = "Posterior ACR means with high similarity assumption")
axis(1, at = c(0.7,x, 0.9), labels = c(0.7,x,0.9))
segments(x[1],0,x[1],0.3, col = 'red',lwd = 3)
segments(x[2],0,x[2],0.3, col = 'blue',lwd = 3)