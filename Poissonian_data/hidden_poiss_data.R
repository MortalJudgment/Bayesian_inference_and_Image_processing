# Initialization 
y=rep(0,13)
lambda=1
for (i in 1:13)
{while(y[i]<4) y[i]=rpois(1, lambda)}
niter=2000
lambdas=rep(0,niter) #for the Non-Rao-Blackwell estimates
lambdasest=rep(0,niter) #for the Rao-Blackwell estimates
 for(j in 1:niter)
 {
   
 lambda=rgamma(1,shape = 313+sum(y), scale = 1/360) 
 y=rep(0,13)
 for (i in 1:13)
  {while(y[i]<4) y[i]=rpois(1, lambda)}
  lambdahat=(313+sum(y))/360
  lambdasest[j]=lambdahat
  lambdas[j]=lambda
 }
 
 # Convergence 
 
 par(mfrow=c(1,3))
 moyLam=cumsum(lambdasest)/(1:niter) # Empir avg of Rao Blackwell estim
 moylam=cumsum(lambdas)/(1:niter) #Empirical average of lambdas
 plot(moylam, type='l', ylab='', xlab="iterations")
 lines(moyLam, col='red')
 
 # Post burn-in estimates and histograms
 mean(lambdasest[1000:2000]) # Rao-Blackwell estimate
 mean(lambdas[1000:2000]) # Usual estimate
 
 hist(lambdas[1000:2000],main='Empirical Average', xlab="")
 hist(lambdasest[1000:2000],main='Rao-Blackwell', xlab="")
 par(mfrow=c(1,1))
