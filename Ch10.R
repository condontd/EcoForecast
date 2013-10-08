### Chapter 10

### Meta-analysis stats
library(MCMCpack)
mu = 5  ## true mean
sigma = 4 ## true standard deviation
N = 10^(1:4)
n = 10000
par(mfrow=c(4,2))
for(i in 1:4){
  T = S = rep(NA,n)
  for(j in 1:n){
    x = rnorm(N[i],mu,sigma)
    T[j] = mean(x)
    S[j] = var(x)
  }
  hist(T,probability=TRUE)
  xseq = seq(min(T),max(T),length=1000)
  lines(xseq,dnorm(xseq,mu,sigma/sqrt(N[i])))
  hist(S,probability=TRUE)
  sseq = seq(min(S),max(S),length=1000)
  lines(sseq,dinvgamma(sseq,N[i]/2,N[i]/2*sigma^2))  
}