
#Metropolis Samples (Symmetric Proposal Distribution)
set.seed(575)
#initialize variables and constants
data = data.frame(enc=0:16,freq=c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1))
N = sum(data$freq)
y = rep(data$enc,data$freq)
reject <- 0
alpha = rep(0,N)
beta = rep(0,N)
mu = rep(0,N)
lambda = rep(0,N)
minn = -1
maxx = 1

#functions
log.likelihood <- function(alpha,beta,mu,lambda,x){
  l = 0
  #reparameterization
  alpha = exp(alpha)/(1+exp(alpha))
  beta = exp(beta)/(1+exp(beta))
  mu = exp(mu)/(1+exp(mu))
  lambda = exp(lambda)/(1+exp(lambda))
  for(i in 1:length(x$enc))
  {
    e = x$enc[i]
    n = x$freq[i]
    if(e==0){
      l = l + n*log(alpha + beta*exp(-mu) + (1-alpha-beta)*exp(-lambda))
      print(l)
    }
    else{
      l =l + n*log(beta*(mu^e)*exp(-mu) + (1-alpha-beta)*exp(-lambda)*lambda^e)-log(factorial(e))
      print(l)
    }
  }
  return(l)
}

#initialize sample values i.e chains
alpha[1] = rnorm(1,0,1)
beta[1] = rnorm(1,0,1)
mu[1] = rnorm(1,0,1)
lambda[1] = rnorm(1,0,1)
sigma = 1
i = 1

for(i in 2:N){
  
  #sample from proposal distribution (symmetrical)
  alpha.star <- rnorm(1,mean = alpha[i-1], sd = sigma)
  beta.star <- rnorm(1,mean = beta[i-1], sd = sigma)
  mu.star <- rnorm(1,mean = mu[i-1], sd = sigma)
  lambda.star <- rnorm(1,mean = lambda[i-1], sd = sigma)
  print(c(alpha.star,beta.star,mu.star,lambda.star))
  
  num <- log.likelihood(alpha.star,beta.star,mu.star,lambda.star,data)
  dem <- log.likelihood(alpha[i-1],beta[i-1],mu[i-1],lambda[i-1],data)
  ratio <- exp(num) - exp(dem)
  accept.prob <- min(1,ratio)
  U = runif(1)
  
  if(U <= accept.prob){
    alpha[i] <- alpha.star
    beta[i] <- beta.star
    mu[i] <- mu.star
    lambda[i] <- lambda.star
  }
  else{
    alpha[i] <- alpha[i-1]
    beta[i] <- beta[i-1]
    mu[i] <- mu[i-1]
    lambda[i] <- lambda[i-1]
    reject = reject + 1
  }
}
# alpha = exp(alpha)/(1+exp(alpha))
# beta = exp(beta)/(1+exp(beta))
# mu = exp(mu)
# lambda = exp(lambda)
c(mean(alpha),mean(beta),mean(mu),mean(lambda))
print("Rejection Rate:")
100*(reject/N)

