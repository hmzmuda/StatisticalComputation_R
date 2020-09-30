x<-1
while(x <= 10){
  U2 <- runif(1,0,1)
  print(U2)
  x = x+1
  print(x)
}


set.seed(475)
n=10; 
U1=runif(n,min=0,max=1); 
U2=runif(n,min=0,max=1)
M=2*exp(-1/2)*sqrt(pi/2)
print(M)
Y=tan(pi*(U1-1/2))
fY=exp(-Y^2/2)/sqrt(2*pi)
print(fY)
gY=1/pi/(1+Y^2)
print(gY)
X=Y[U2<=fY/gY/M]
print(X)

hist(X,freq=FALSE,ylim=c(0,0.5), main=paste("Acceptance Rate:",length(X)/n))
tt=seq(from=-4,to=4,by=0.01); 
lines(tt,exp(-tt^2/2)/sqrt(2*pi),col='red') ## add the true normal pdf
lines(tt,1/pi/(1+tt^2)*M,col='green') 

