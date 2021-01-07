
#Metropolis Samples (Symmetric Proposal Distribution)
set.seed(575)
#initialize variables and constants
df = data.frame(e=0:16,f=c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1))
N = sum(df$f)
#N = 50
reject <- 0
alpha = rep(0,N)
beta = rep(0,N)
mu = rep(0,N)
lambda = rep(0,N)

#functions
pi_i <- function(a, b, m, l,i) {
  if(i==0){
    return(a + b*exp(-m) + (1-a-b)*exp(-l))
    }
  else{
    return(b*(m^i)*exp(-m) + (1-a-b)*(l^i)*exp(-l))
    }
}
loglike <- function(a,b,m,l,df){
  output = 0
  a = exp(a)/(1+exp(a))
  b = exp(b)/(1+exp(b))
  m = exp(m)
  l = exp(l)
  for(i in 1:length(df)){
    #print(c(i,pi_i(a,b,m,l,i-1)))
    output[i] = df[i]*(log(pi_i(a,b,m,l,i-1))-log(factorial(i-1)))
    if(is.nan(output[i])){output[i] = 0}
  }
  return(sum(output))
}


#initialize sample values i.e chains
alpha[1] = rnorm(1,0,1)
beta[1] = rnorm(1,0,1)
mu[1] = rnorm(1,0,1)
lambda[1] = rnorm(1,0,1)
sigma = 2
i = 1
U = rnorm(N,0,sigma)
count = 0
reject = 0

for(i in 2:N){
  
  #sample from proposal distribution (symmetrical)
  alpha.star <- rnorm(1,mean = alpha[i-1], sd = sigma)
  beta.star <- rnorm(1,mean = beta[i-1], sd = sigma)
  mu.star <- rnorm(1,mean = mu[i-1], sd = sigma)
  lambda.star <- rnorm(1,mean = lambda[i-1], sd = sigma)
  
  num <- loglike(alpha.star,beta.star,mu.star,lambda.star,df$f)
  dem <- loglike(alpha[i-1],beta[i-1],mu[i-1],lambda[i-1],df$f)
  ratio <- exp(num-dem)
  accept.prob <- min(1,ratio)
  
  if(U[i] < accept.prob){
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
  count = count + 1
}
# alpha = exp(alpha)/(1+exp(alpha))
# beta = exp(beta)/(1+exp(beta))
# mu = exp(mu)
# lambda = exp(lambda)
c(mean(alpha),mean(beta),mean(mu),mean(lambda))
print("Rejection Rate:")
100*(reject/count)

