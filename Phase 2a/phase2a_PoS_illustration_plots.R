# plots that demonstrate how the PoS is calculated
m = -0.249
s = 0.151
t = 0.75*sqrt(2/75)

x <- seq(m-3*s,m+3*s,6*s/500)


######### 1
# piror distribution
dx <- dnorm(x,m,s)
plot(x,dx,type = 'l', ylim = range(0,1.2*max(dx)), xlab = 'mean log ACR', ylab = 'density', main = 'prior distribution of the means')
segments(m,0,m,max(dx),lty=3)

######### 2
# distributions around points on the prior distribution
dx <- dnorm(x,m,s)
plot(x,dx,type = 'l', ylim = range(0,1.2*max(dx)), xlab = 'mean log ACR', ylab = 'density', main = 'Simulation of estimated log ACR')
segments(m,0,m,max(dx),lty=3)


mult = -1
my <- m+mult*t
index = x>my - 2*t & x < my + 2*t
dy <- dnorm(x,mean = my,sd=t)/3
lines(x[index],dy[index], col = 'blue')
segments(my,min(dy[index]),my,max(dy), lty=3, col = 'blue')


mult = -0.1
my <- m+mult*t
index = x>my - 2*t & x < my + 2*t
dy <- dnorm(x,mean = my,sd=t)/3
lines(x[index],dy[index], col = 'blue')
segments(my,min(dy[index]),my,max(dy), lty=3, col = 'blue')


mult = 0.6
my <- m+mult*t
index = x>my - 2*t & x < my + 2*t
dy <- dnorm(x,mean = my,sd=t)/3
lines(x[index],dy[index], col = 'blue')
segments(my,min(dy[index]),my,max(dy), lty=3, col = 'blue')


mult = 0.8
my <- m+mult*t
index = x>my - 2*t & x < my + 2*t
dy <- dnorm(x,mean = my,sd=t)/3
lines(x[index],dy[index], col = 'blue')
segments(my,min(dy[index]),my,max(dy), lty=3, col = 'blue')


###### 3
# posterior distribution
dfinal <- dnorm(x, mean = m, sd = sqrt(s^2+t^2))
plot(x, dfinal, col='blue', type = 'l',ylim = range(0,1.2*max(dx)), xlab = 'estimated log ACR', ylab = 'density', main = 'posterior distribution of the estimated log ACR')
segments(m,min(dfinal),m,max(dfinal), lty=3,col='blue')



######### 4
# PoS calculated based on thresholds
dfinal <- dnorm(x, mean = m, sd = sqrt(s^2+t^2))
plot(x, dfinal, col='blue', type = 'l',ylim = range(0,1.2*max(dx)), xlab = 'estimated log ACR', ylab = 'density', main = 'Calculate PoS based on threshold')
segments(m,min(dfinal),m,max(dfinal), lty=3,col='blue')

thresh = qnorm(0.1)*t

segments(thresh,min(dfinal),thresh,dfinal[which.min(abs(x-thresh))] ,lwd=3,col='blue')
polygon(c(x[x<=thresh], thresh),  c(dfinal[x<=thresh],min(dfinal) ), col="red")
text(thresh, 2.5,'threshold = -0.157')


######### joint distribution (2d)

ellipse <- read.csv('H:\\Desktop\\R\\data\\simulation_of_1k.csv')
x1 <- ellipse$gamma0hat

y1 <- ellipse$theta0


plot(x1,y1,pch=20, col='blue',xlim = range(-1,1),ylim=range(-1,1),xlab = 'estimated log ACR', ylab = 'log HR', main = 'conditional joint distribution of log HR and log ACR')
abline(v=thresh,h=log(0.8),lty=3)
text(-0.85,c(-0.7,0.45),c('A','B'), cex = 2)
text(0.7,c(-0.7,0.45),c('D','C'), cex = 2)

