
#### ERRORS IN VARIABLES EXAMPLE: TDR + GROWTH

n <- 100
TDR.sd = 0.05
theta <- runif(n,0.12,0.40)
theta.calib <- seq(0.12,0.40,length=20)
TDR   <- rnorm(n,2*theta - 0.2,TDR.sd)
TDR.calib <- rnorm(20,2*theta.calib - 0.2,TDR.sd)

plot(TDR,theta)
points(TDR.calib,theta.calib,col=2,pch=18)
calib.curve = lm(theta.calib ~ TDR.calib)
abline(coef(calib.curve))
abline(0.1,0.5,col=2)
summary(calib.curve)

tmin=0.12
Egrow = 10*(theta-tmin)/((theta-tmin)+(0.2-tmin))
plot(theta,Egrow)
y <- rnorm(n,Egrow,1)
plot(theta,y)
points(theta,Egrow,col=2,pch=3)

plot(TDR,y)
library(rjags)
library(R2WinBUGS)
### With latent X
model <- function(){
  
  ### calibration curve
  for(i in 1:2) { alpha[i] ~ dnorm(0,0.001)}
  sigma ~ dgamma(2,0.005)
  for(i in 1:20){
    ESMc[i] <- alpha[1] + alpha[2]*TDRc[i]
    SMc[i] ~ dnorm(ESMc[i],sigma)
    PSMc[i] ~ dnorm(ESMc[i],sigma)
  }
  
  ## priors
  beta[1] ~ dlnorm(2.3,0.01)
  beta[2] ~ dlnorm(-2,0.01)
  beta[3] ~ dlnorm(-3.5,0.01)
  tau ~ dgamma(2,2)
  #  r ~ dunif(0,100)    
  for(i in 1:n){
    ESM[i] <-  alpha[1] + alpha[2]*TDR[i]
    SM[i] ~ dnorm(ESM[i],sigma)
    mu[i] <- beta[1]*(SM[i]-beta[2])/(SM[i]+beta[3])
#    log(mu[i]) <- beta[1]+beta[2]*SM[i]
    #    p[i] <- r/(mu[i]+r)
    #    y[i] ~ dnegbin(p[i],r)
    #    py[i] ~ dnegbin(p[i],r)
    y[i] ~ dnorm(mu[i],tau)
    py[i] ~ dnorm(mu[i],tau)
  }
}

write.model(model,"EIV.jags")

mod <- jags.model("EIV.jags",data=list(TDR=TDR,y=y,n=n,TDRc=TDR.calib,
                  SMc=theta.calib),n.adapt=1000,n.chains=3,
                  init=list(beta=c(10,0.1,0.02),tau=1,sigma=1/TDR.sd^2))
jdat <- coda.samples(mod,variable.names=c("beta"),n.iter=3000) ## burnin
plot(jdat)
jdat <- coda.samples(mod,variable.names=c("alpha","beta","tau","sigma","mu","py","SM","ESMc","PSMc"),
                     n.iter=30000) ## samples
dat <- as.matrix(jdat)
pi <- apply(dat,2,quantile,c(0.025,0.5,0.975))
#dic <- dic.samples(mod,5000,type="pD")
plot(as.mcmc(dat[,n+40+1:5]))
mu <- pi[,5+n+1:n+40]
py <- pi[,5+2*n+1:n+40]
SM <- pi[,1:n+40]
ESMc <- pi[,1:20]
PSMc <- pi[,20+1:20]

ciEnvelope <- function(x,ylo,yhi,...){
  polygon(cbind(c(x, rev(x), x[1]), c(ylo, rev(yhi),
                                      ylo[1])), border = NA,...) 
}
par(lwd=2,cex=1.2)

alpha <- pi[2,n+1:2+40]
SMbar <- alpha[1] + alpha[2]*TDR
#plot(SM[2,],y)
plot(SMbar,y,pch=2,type='n')
legend("topleft",legend=c("Obs","Est","True",'Median',"CI","PI"),
       pch=c(2,18,NA,NA,NA,NA),lty=c(NA,NA,1,1,1,1),col=c(1,2,1,2,"lightblue","lightpink"),
       lwd=c(1,1,3,3,16,16))
ord <- order(SM[2,])
ciEnvelope(SM[2,ord],py[1,ord],py[3,ord],col="lightpink")
ciEnvelope(SM[2,ord],mu[1,ord],mu[3,ord],col="lightblue")
points(SMbar,y,pch=2)
lines(theta[order(theta)],Egrow[order(theta)],lwd=3)
points(SM[2,],y,col=2,pch=18)
lines(SM[2,ord],mu[2,ord],col=2,lwd=3)


#par(mfrow=c(2,1),mar=c(3.5,4,2,1))
ord2 <- order(ESMc[2,])
plot(TDR.calib,theta.calib,type='n',ylim=range(c(range(ESMc),theta.calib)),
     xlab="TDR",ylab="Soil moisture",main="Calibration",
     cex.lab=1.5,cex=2,cex.main=2,mgp=c(2.4,0.9,0))
legend("topleft",legend=c("Obs",'Median',"CI","PI"),
       pch=c(18,NA,NA,NA),lty=c(NA,1,1,1),col=c(1,1,"darkgrey","lightgrey"),
       lwd=c(1,3,16,16))
ciEnvelope(TDR.calib[ord2],PSMc[1,ord2],PSMc[3,ord2],col="lightgrey")
ciEnvelope(TDR.calib[ord2],ESMc[1,ord2],ESMc[3,ord2],col="darkgrey")
lines(TDR.calib[ord2],ESMc[2,ord2],lwd=3)
points(TDR.calib,theta.calib,pch=18)
#abline(0.1,0.5,col=2)


### Without latent X
model <- function(){
  
  ## priors
  beta[1] ~ dlnorm(2.3,0.01)
  beta[2] ~ dlnorm(-2,0.01)
  beta[3] ~ dlnorm(-3.5,0.01)
  tau ~ dgamma(2,2)
  #  r ~ dunif(0,100)    
  for(i in 1:n){
    mu[i] <- beta[1]*(SM[i]-beta[2])/(SM[i]+beta[3])
    y[i] ~ dnorm(mu[i],tau)
    py[i] ~ dnorm(mu[i],tau)
  }
}

write.model(model,"Monod.jags")
SMc = predict(lm(theta.calib~TDR.calib),newdata=data.frame(TDR.calib=TDR))
mod <- jags.model("Monod.jags",data=list(y=y,n=n,SM = SMc),n.adapt=1000,n.chains=3,
                  init=list(beta=c(10,0.1,0.02),tau=1))
jdat <- coda.samples(mod,variable.names=c("beta"),n.iter=3000) ## burnin
plot(jdat)
jdat <- coda.samples(mod,variable.names=c("beta","tau","mu","py"),
                     n.iter=30000) ## samples
dat2 <- as.matrix(jdat)
pi2 <- apply(dat2[10000:nrow(dat2),],2,quantile,c(0.025,0.5,0.975))
mu2 <- pi2[,3+1:n]
py2 <- pi2[,3+n+1:n]

plot(SMbar,y,pch=2,type='n',main="Growth Response to Moisture",
     xlab=expression(paste("Soil Moisture ",(m^3/m^3))),ylab="Growth (cm/yr)",
     cex.lab=1.5,cex=2,cex.main=2,mgp=c(2.4,0.9,0))
legend("topleft",legend=c("Obs","True","no EIV","EIV"),
       pch=c(18,NA,NA,NA),lty=c(NA,1,1,1),col=c(1,1,4,2),
       lwd=c(1,4,4,4))
ord <- order(SM[2,])
ciEnvelope(SM[2,ord],mu[1,ord],mu[3,ord],col="lightpink")
ciEnvelope(SMc[order(SMc)],mu2[1,order(SMc)],mu2[3,order(SMc)],col="lightblue")
points(SMbar,y,pch=18)
lines(theta[order(theta)],Egrow[order(theta)],lwd=4)
lines(SM[2,ord],mu[2,ord],col=2,lwd=4)
lines(SMc[order(SMc)],mu2[2,order(SMc)],col=4,lwd=4)

save(theta,theta.calib,TDR,TDR.calib,y,pi,pi2,file="Ch9.Rdata")

#########################################################################

####    MISSING DATA EXAMPLE

## Missing Data
a = 10
b = -0.3
v = 0.5
n = 25

xtrue <- runif(n,0,10)
ey <- a + b *xtrue
ytrue <- rnorm(n,ey,v)

x <- c(xtrue,NA)
y <- c(ytrue,7.5)
mis <- n+1
n <- n+1


par(lwd=2,cex=1.5)
plot(x,y)
abline(a,b)
#for(i in 1:n){lines(rep(x[i],2),c(y[i],ybar[i]),lty=2)}


### JAGS
library(rjags)
library(R2WinBUGS)
model <- function(){
  ## priors
  for(i in 1:2) { beta[i] ~ dnorm(0,0.001)}  
  sigma ~ dgamma(0.1,0.1)
  for(i in mis) { x[i] ~ dunif(0,10)}
  
  for(i in 1:n){
    mu[i] <- beta[1]+beta[2]*x[i]
    mup[i] <- beta[1]+beta[2]*xp[i]
    y[i] ~ dnorm(mu[i],sigma)
    yp[i] ~ dnorm(mup[i],sigma)
  }
}

write.model(model,"mis.jags")
xp = seq(0,10,length=n)
mod <- jags.model("mis.jags",data=list(x=x,y=y,mis=mis,n=n,xp=xp),n.adapt=1000,n.chains=3)
jdat <- coda.samples(mod,variable.names=c("beta","sigma","x[26]"),n.iter=1000) ## burnin
plot(jdat)

jdat <- coda.samples(mod,variable.names=c("beta","sigma","x[26]","mup","yp"),n.iter=30000) ## samples
dat <- as.matrix(jdat)
pi <- apply(dat,2,quantile,c(0.025,0.5,0.975))
#dic <- dic.samples(mod,5000,type="pD")
CI = pi[,3:28]
PI = pi[,31:56]
x26 = pi[,30]
  
par(lwd=2,cex=1.2)
plot(x,y,xlim=c(0,10),ylim=range(CI),type='n',cex.main=2,main="Missing Data Model",cex.lab=1.5,mgp=c(2.4,0.9,0))
ciEnvelope(xp,PI[1,],PI[3,],col="lightgrey")
ciEnvelope(xp,CI[1,],CI[3,],col="darkgrey")
lines(xp,CI[2,],lwd=4)
#abline(pi[2,1],pi[2,2],lwd=4)
#abline(a,b)
abline(h=7.5,col=3,lty=2,lwd=3)
xmis <- density(dat[,30])
xmis$y[xmis$x > 10] = 0
lines(xmis$x,xmis$y*4 + 7.5,col=3,lwd=8)
points(x,y,pch=18)