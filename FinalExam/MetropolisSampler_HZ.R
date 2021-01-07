
# #Metropolis Samples (Symmetric Proposal Distribution)
# set.seed(575)
# #initialize variables and constants
# df = data.frame(e=0:16,f=c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1))
# N = sum(df$f)
# #N = 50
# reject <- 0
# alpha = rep(0,N)
# beta = rep(0,N)
# mu = rep(0,N)
# lambda = rep(0,N)
# 
# #functions
# pi_i <- function(a, b, m, l,i) {
#   if(i==0){
#     return(a + b*exp(-m) + (1-a-b)*exp(-l))
#     }
#   else{
#     return(b*(m^i)*exp(-m) + (1-a-b)*(l^i)*exp(-l))
#     }
# }
# loglike <- function(a,b,m,l,df){
#   output = 0
#   a = exp(a)/(1+exp(a))
#   b = exp(b)/(1+exp(b))
#   m = exp(m)
#   l = exp(l)
#   for(i in 1:length(df)){
#     #print(c(i,pi_i(a,b,m,l,i-1)))
#     output[i] = df[i]*(log(pi_i(a,b,m,l,i-1))-log(factorial(i-1)))
#     if(is.nan(output[i])){output[i] = 0}
#   }
#   return(sum(output))
# }
# 
# 
# #initialize sample values i.e chains
# alpha[1] = rnorm(1,0,1)
# beta[1] = rnorm(1,0,1)
# mu[1] = rnorm(1,0,1)
# lambda[1] = rnorm(1,0,1)
# sigma = 2
# i = 1
# U = rnorm(N,0,sigma)
# count = 0
# reject = 0
# 
# for(i in 2:N){
#   
#   #sample from proposal distribution (symmetrical)
#   alpha.star <- rnorm(1,mean = alpha[i-1], sd = sigma)
#   beta.star <- rnorm(1,mean = beta[i-1], sd = sigma)
#   mu.star <- rnorm(1,mean = mu[i-1], sd = sigma)
#   lambda.star <- rnorm(1,mean = lambda[i-1], sd = sigma)
#   
#   num <- loglike(alpha.star,beta.star,mu.star,lambda.star,df$f)
#   dem <- loglike(alpha[i-1],beta[i-1],mu[i-1],lambda[i-1],df$f)
#   ratio <- exp(num-dem)
#   accept.prob <- min(1,ratio)
#   
#   if(U[i] < accept.prob){
#     alpha[i] <- alpha.star
#     beta[i] <- beta.star
#     mu[i] <- mu.star
#     lambda[i] <- lambda.star
#   }
#   else{
#     alpha[i] <- alpha[i-1]
#     beta[i] <- beta[i-1]
#     mu[i] <- mu[i-1]
#     lambda[i] <- lambda[i-1]
#     reject = reject + 1
#   }
#   count = count + 1
# }
# # alpha = exp(alpha)/(1+exp(alpha))
# # beta = exp(beta)/(1+exp(beta))
# # mu = exp(mu)
# # lambda = exp(lambda)
# c(mean(alpha),mean(beta),mean(mu),mean(lambda))
# print("Rejection Rate:")
# 100*(reject/count)

set.seed(5)
df <- c(379,299,222,145,109,95,73,59,45,30,24,12,4,2,0,1,1)
alpha <- 0.23
beta <- 0.25
mu <- 3.5
lambda <- 3.5
N <- 1000
chain <- matrix(nrow = N, ncol = 4)
chain[1,] <- c(alpha,beta,mu,lambda)
minn = -0.5
maxx = 0.5
reject = 0
pi_i <- function(alpha,beta,mu,lambda,i){
  if(i==0) return(alpha + beta*exp(-mu) + (1-alpha-beta)*exp(-lambda))
  else return(beta*(mu^i)*exp(-mu) + (1-alpha-beta)*(lambda^i)*exp(-lambda))
}
loglike <- function(alpha,beta,mu,lambda){
  sum.out <- 0
  alpha = exp(alpha)/(1+exp(alpha))
  beta = exp(beta)/(1+exp(beta))
  mu = exp(mu)
  lambda = exp(lambda)
  for(i in 1:length(df)){
    sum.out <- sum.out + df[i]*(log(pi_i(alpha,beta,mu,lambda,i-1)) - log(factorial(i-1)))
  }
  return(sum.out)
}
print(pi_i(alpha,beta,mu,lambda,0))
print(loglike(alpha,beta,mu,lambda))
print(exp(loglike(alpha,beta,mu,lambda)))

for(i in 2:N){
  alpha.star <- chain[i-1,1] + runif(1,minn,maxx)
  while(alpha.star < 0 || alpha.star > 1){
    alpha.star <- chain[i-1,1] + runif(1,minn,maxx)
  }
  beta.star <- chain[i-1,2] + runif(1,minn,maxx)
  while(beta.star < 0 || beta.star > 1){
    beta.star <- chain[i-1,2] + runif(1,minn,maxx)
  }
  mu.star <- chain[i-1,3] + runif(1,-0.1,0.1)
  lambda.star <- chain[i-1,4] + runif(1,-0.1,0.1)
  ratio <- exp(loglike(alpha.star,beta.star,mu.star,lambda.star) - loglike(chain[i-1,1],chain[i-1,2],chain[i-1,3],chain[i-1,4]))
  if(runif(1) < ratio){
    chain[i,] <- c(alpha.star,beta.star,mu.star,lambda.star)
  }
  else{
    chain[i,] <- chain[i-1,]
    reject <- reject + 1
  }
}
  


x = rep(0,10000)
x[1] = 3     #initialize; I've set arbitrarily set this to 3

target = function(x){
  return(ifelse(x<0,0,exp(-x)))
}

for(i in 2:10000){
  current_x = x[i-1]
  proposed_x = current_x + rnorm(1,mean=0,sd=1)
  A = target(proposed_x)/target(current_x) 
  if(runif(1)<A){
    x[i] = proposed_x       # accept move with probabily min(1,A)
  } else {
    x[i] = current_x        # otherwise "reject" move, and stay where we are
  }
}


