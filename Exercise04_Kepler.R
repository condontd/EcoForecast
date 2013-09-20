
## pull data
URL = "http://phenocam.sr.unh.edu/data/archive/uiefprairie/ROI/uiefprairie_canopy_0001_gcc.csv"
dat <- read.csv(URL,skip=6)
dat$date <- as.Date(as.character(dat$date))
fname <- strsplit(URL,"/",fixed=TRUE)[[1]]
fname = fname[length(fname)]  ## grab the last part of the file
fname = sub("csv","RData",fname)
save(dat,file=fname)
#n.samp <- dat$nimage
#gcc <- dat$gcc_mean

## extreme value distribution
## order statistics: f_z = n*f_x(z)*F_x^(n-1)(z)
onorm <- function(z,mu,sd,n,log=FALSE){
  if(log){
    x = log(n)+dnorm(z,mu,sd,log=TRUE)+(n-1)*pnorm(z,mu,sd,log.p=TRUE)
  }else{
    x = n*dnorm(z,mu,sd)*pnorm(z,mu,sd)^(n-1) 
  }
  return(x)
}
z = seq(-1,4,length=1000)
plot(z,onorm(z,0,1,1),lwd=3,ylab="density",ylim=c(0,.8),type='l')
n = c(1,3,6,9,18)
for(i in 2:5){
  lines(z,onorm(z,0,1,n[i]),col=i,lwd=3)
}
legend("topright",legend=n,lwd=3,col=1:5)

sd.LL <- function(sd,mean,min,max,n){
  -onorm(max-mean,0,sd,n,log=TRUE)-onorm(mean-min,0,sd,n,log=TRUE)
}

sd.mle <- function(mean,min,max,n){
  k = length(n)
  sd = rep(NA,k)
  for(i in 1:k){
    if(all(!is.na(c(mean[i],min[i],max[i],n[i])))){
      fit = optimize(sd.LL,c(0,0.25),mean=mean[i],min=min[i],max=max[i],n=n[i])   
      sd[i]=  fit$minimum
    }
  }
  return(sd)
}
sd = sd.mle(dat$gcc_mean,dat$gcc_min,dat$gcc_max,dat$nimage)
hist(sd)
plot(dat$nimage,sd)
plot(date,sd)
plot(dat$gcc_mean,sd)
het = lm(sd~dat$gcc_mean)
abline(het,col=2,lwd=4)
sd_bar = mean(sd,na.rm=TRUE)
sd_lm = coef(het)

## split/subset data

## plot data
ciEnvelope <- function(x,ylo,yhi,col="lightgrey",...){
  has.na = apply(is.na(cbind(x,ylo,yhi)),1,sum)
  block = cumsum(has.na);block[has.na>0] = NA
  for(i in unique(block)){
    sel = which(block==i)
  polygon(cbind(c(x[sel], rev(x[sel]), x[sel[1]]), c(ylo[sel], rev(yhi[sel]),
                                      ylo[sel[1]])), col=col,border = NA,...) 
  }
}


plot(dat$date,dat$gcc_mean,type='l')
ciEnvelope(dat$date,(1-sd_lm[2])*dat$gcc_mean-sd_lm[1],(1+sd_lm[2])*dat$gcc_mean+sd_lm[1])
lines(date,dat$gcc_mean,lwd=1.5)
#lines(date,dat$gcc_min,col=3)
#lines(date,dat$gcc_max,col=2)


## fit logit