## AR likelihood

library(mvtnorm)

H <- matrix(c(0,1,2,1,0,1,2,1,0),3,3)
n = 100
H = as.matrix(dist(1:n,diag=TRUE,upper=TRUE))

g = rep(1,n)
g= 1:n

lnlik <- function(theta){
  beta  <- theta[1]
  sigma <- theta[2]
  rho   <- theta[3]
  SIGMA <- sigma/(1-rho^2)*rho^H
  
  L = dmvnorm(g,rep(beta,n),SIGMA,log=TRUE)
  
  return(L)
  
}

lnlik(c(1,1,0))

rho = seq(0,0.99999,length=100)
L = rep(NA,100)
for(i in 1:100){
  L[i] = lnlik(c(1,1,rho[i])) 
}

plot(rho,-L)

plot(rho,(1-rho)/(1+rho)*-L[1])
lines(rho,-L)
abline(h=-dnorm(0,log=TRUE))

-dnorm(0,log=TRUE)
min(-L)
neff = L/L[1]*n
plot(rho,neff)
plot(rho,neff,log='y')
lines(rho,(1-rho)/(1+rho)*n,col=2)
abline(h=1,col=3)
rho[which.min((neff-1)^2)]
ar(g)

i = matrix(1,1,1)
dmvnorm(1,1,i,log=TRUE)

fit <- optim(c(0,0,1,rho),lnlik)

dat = read.csv("~/Dropbox/GE375_2013/Labs/GE375 Lab 6 - Model Assessment/AMF_USWCr_2002_L2_WG_V004.csv",skip=18,na.strings=c("-9999","-6999"))
dat = read.csv("~/stats/US-WCr-Clean-2012-L0-vMar2013.csv",skip=126,na.strings=c("-999","-9999"))
dat[dat==-999] = NA

hist(dat$par_30)
plot(dat$par_30,dat$par_2)

sel = which(dat$par_30 > 700)
fpar = dat$par_2/dat$par_30
fpar[sel] = NA
plot(fpar,ylim=c(0,1))
plot(exp(-fpar))
plot(dat$par_2)
ar(fpar,FALSE,order.max=1,na.action=na.exclude)
hist(fpar)

LI = -log(fpar)  ## unscaled Beer's Law inversion to LAI ('leaf index')
fit.ar = ar(LI,FALSE,order.max=1,na.action=na.exclude)
(1-fit.ar$ar)/(1+fit.ar$ar)

(1-fit.ar$ar)/(1+fit.ar$ar)*sum(!is.na(fpar))/366  ## effective meas/day w/o nugget correction

## need to see if we can fit the latent AR term out of a noisy TS (nugget correction)


## need to see what happens to the likelihood if we hold the length in time constant
## but change the sampling frequency (which should also increase rho)
xlim = c(0,100)
x = list()
n.seq = c(2,3,4,6,8,12,16,24,32,48,64,128,256,512)
for(i in 1:length(n.seq)){x[[i]] = seq(xlim[1],xlim[2],length=n.seq[i])}
rho = sapply(x,function(x){ar(x,FALSE,order.max=1)$ar})
plot(n.seq,rho,type='l')
plot(n.seq,rho,type='l',log='x')
plot(n.seq,1-rho,type='l',log='xy')   ## rho approaches 1 linearly on a log-log scale!!
lm(log(1-rho)~log(n.seq))             ## intercept is 1.1, slope is -1


lnlik.sub <- function(theta,x,rho){
  #  beta  <- theta[1]
  sigma <- theta[1]
  #  rho   <- theta[3]
  SIGMA <- sigma/(1-rho^2)*rho^H
  
  L = dmvnorm(x,rep(mean(xlim),length(x)),SIGMA,log=TRUE)
  
  return(L)
  
}


Lsub.1 = NA
Lsub.v = NA
Lsub.r = NA
for(i in 1:length(n.seq)){
  H = as.matrix(dist(1:n.seq[i],diag=TRUE,upper=TRUE))
  Lsub.1[i] = lnlik.sub(1,x=x[[i]],rho=rho[i])
  Lsub.v[i] = lnlik.sub(var(x[[i]]),x=x[[i]],rho=rho[i])
  Lsub.r[i] = lnlik.sub(1,x=x[[i]],rho=0)
}
plot(n.seq,Lsub.1) ## linear?
plot(n.seq,-Lsub.r,log='xy') ## linear
lines(n.seq,-Lsub.1)


plot(n.seq,Lsub.r/n.seq,log='x')
lines(n.seq,Lsub.1/n.seq)


plot(n.seq,(1-rho)/(1+rho)*n.seq)
plot(n.seq,(1-rho)/(1+rho)*n.seq,log="xy")  
## this doesn't make sense intuitively the Neff should rise from n to an asymptote, not decline
rho2 = 1-exp(1.1)/n.seq
lines(n.seq,(1-rho2)/(1+rho2)*n.seq)

#what it x is a sine wave?
xlim = c(0,2*pi)
x = y = list()
n.seq = c(6,8,12,16,24,32,48,64,128,256,512)
for(i in 1:length(n.seq)){x[[i]] = seq(xlim[1],xlim[2],length=n.seq[i])}
rho = sapply(x,function(x){ar(sin(x),FALSE,order.max=1)$ar})
plot(n.seq,rho,type='l')
plot(n.seq,rho,type='l',log='x')
plot(n.seq,1-rho,type='l',log='xy')   ## rho approaches 1 linearly on a log-log scale!!
lm(log(1-rho)~log(n.seq))             ## intercept is 3.251, slope is -2.05
plot(n.seq,(1-rho)/(1+rho)*n.seq,log="xy") 