## -----------------------------------------------------------------------------
ChangePoint<-function (y,v1,v2,e1,e2,N,burnin){
  mu1<-mu2<-sd1<-sd2<-GAMMA<-R<-prob<-cp<-numeric(N)
  n<-length(y) 
  prob[1]<-1
  R[1]<-sample(1:n,1)
  GAMMA[1]<-rgamma(1,1,1)
  m<-0
  u<-runif(N)
  sd1<-e1^2/(1+GAMMA[1]*e1^2*R[1])
  sd2<-e2^2/(1+GAMMA[1]*e2^2*(n-R[1]))
  mean1<-(v1+GAMMA[1]*e1^2*sum(y[1:R[1]]))/(1+GAMMA[1]*e1^2*R[1])
  mean2<-(v2+GAMMA[1]*e2^2*sum(y[(R[1]+1):n]))/(1+GAMMA[1]*e2^2*(n-R[1]))
  mu1[1]<-rnorm(1,mean=mean1,sd=sd1)
  mu2[1]<-rnorm(1,mean=mean2,sd=sd2)
for (i in 2:N) {
#generate the change point r using the Metropolis-Hastings algorithm for the complex conditional 
#distribution of r given other parameters is hard to sample from.
R[i]<-sample(1:n,1)
kt<-R[i-1]
if (kt<R[i]) {prob[i-1]<-exp(GAMMA[i-1]/2*(sum((y[(kt+1):R[i]]-mu2[i-1])^2)-sum((y[(kt+1):R[i]]-mu1[i-1])^2)))}
else {prob[i-1]<-exp(GAMMA[i-1]/2*(sum((y[min(R[i-1]+1,n):kt]-mu1[i-1])^2)-sum((y[min(R[i-1]+1,n):kt]-mu2[i-1])^2)))}

if (u[i]<=min(prob[i-1],1))  {r<-R[i];  m<-m+1;}
else {r<-kt}
cp[i]<-r;


if (r==n){
mean1<-(v1+GAMMA[i-1]*e1^2*sum(y[1:r]))/(1+GAMMA[i-1]*e1^2*r)
mean2<-v2+GAMMA[i-1]*e2^2*y[n]
sd1<-e1^2/(1+GAMMA[i-1]*e1^2*r)
sd2<-e2^2

mu1[i]<-rnorm(1,mean=mean1,sd=sd1)
mu2[i]<-rnorm(1,mean=mean2,sd=sd2)
z<-1+1/2*sum((y[1:r]-mu1[i])^2)+1/2*(y[n]-mu2[i])^2
GAMMA[i]<-rgamma(1,shape=n/2+1,rate=z)
}

else {
mean1<-(v1+GAMMA[i-1]*e1^2*sum(y[1:r]))/(1+GAMMA[i-1]*e1^2*r)
mean2<-(v2+GAMMA[i-1]*e2^2*sum(y[(r+1):n]))/(1+GAMMA[i-1]*e2^2*(n-r))
sd1<-e1^2/(1+GAMMA[i-1]*e1^2*r)
sd2<-e2^2/(1+GAMMA[i-1]*e2^2*(n-r))

mu1[i]<-rnorm(1,mean=mean1,sd=sd1)
mu2[i]<-rnorm(1,mean=mean2,sd=sd2)
z<-1+1/2*sum((y[1:r]-mu1[i])^2)+1/2*sum((y[(r+1):n]-mu2[i])^2)
GAMMA[i]<-rgamma(1,shape=n/2+1,rate=z)
}
}
  return(cp)[(burnin+1):N]
}

## -----------------------------------------------------------------------------
Checkfunction1<-function (y,k,p){
  index1<-seq(k,1,-1)
  index2<-seq(1,k,1)
  weight1<-1/(2^index1)
  weight2<-1/(2^index2)
  weight<-c(weight1,0,weight2)
  n<-length(y)
  averageslope<-numeric(n)
  for (i in 1:n){
    if (i>k & i<n-k ){averageslope[i]<-weighted.mean(x=abs(y[(i-k):(i+k)]-y[i]),wt=weight/sum(weight))}
    else {averageslope[i]<--1 }
  }
  return(match(sort(averageslope,decreasing=TRUE),averageslope)[1:p])
  }

## -----------------------------------------------------------------------------
Checkfunction2<-function (y,k,p){
  n<-length(y)
  deviation<-slope1<-slope2<-ratio<-numeric(n)
  for (i in 1:n){
    if (i>k & i<n-k ) {slope1[i]<-(y[i]-y[i-k])/k
    slope2[i]<-(y[i+k]-y[i])/k
    ratio[i]<-(slope1[i]/slope2[i])
    }
  }
  for (i in 1:n){
    if (ratio[i]<0&abs(slope1[i])>median(abs(slope1))&abs(slope2[i])>median(abs(slope2)) )
    {deviation[i]<-max(ratio[i],1/ratio[i])+max(abs(slope1[i]/median(slope1[i])),abs(slope2[i]/median(slope2[i])))       }
    else if (ratio[i]>0) {deviation[i]<-max(ratio[i],1/ratio[i]) }
    else {deviation[i]<--0.001 }
  }
  return (match(sort(deviation,decreasing=TRUE),deviation)[1:p])
}

## -----------------------------------------------------------------------------
set.seed(1234)
x1<-rnorm(40,5,1)
x2<-rnorm(60,10,1)
y<-c(x1,x2)
plot(y,main="The overview of simulation data")
abline(v=length(x1),col="red",lwd=2,lty=2)

## -----------------------------------------------------------------------------
s<-ChangePoint(y,v1=5,v2=10,e1=1,e2=1,N=10000,burnin=1000)
hist(s,freq=FALSE,main="histogram of the change point",xlab="")
dens<-density(s,bw = "nrd0",kernel='gaussian',window=kernel)
plot(dens,main="kernel estimation of the density of change point")
abline(v=dens$x[which.max(dens$y)],col="red")
abline(v=length(x1),col="blue")
results<-matrix(0,nrow=1,ncol=2)
rownames(results)<-"Change point"
colnames(results)<-c("Theoretical","Computed")
results[,1]<-40
results[,2]<-dens$x[which.max(dens$y)]

library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
Checkfunction1(y,3,2)

## -----------------------------------------------------------------------------
Japan<-c(36.42584 ,36.59264 ,36.79944 ,37.03648 ,37.32840 ,37.67568 ,38.09784 ,38.60000 ,39.22192 ,39.97360 ,40.89944 ,42.04432 ,42.61632 ,45.65936 ,48.24248 ,30.53680 ,51.69000 ,56.83000 ,57.71000 ,61.15200 ,61.58300 ,62.42700 ,63.23400 ,64.00500 ,64.73900 ,65.43700 ,66.09900 ,66.72700 ,67.32200 ,67.88700 ,68.42400 ,68.93500 ,69.42600 ,69.89800 ,70.35100 ,70.78600 ,71.20000 ,71.59600 ,71.98000 ,72.35800 ,72.74000)
plot(x=1930:1970,y=Japan,main="The overview of life expectancy of Japan from 1930 to 1970",xlab="",ylab="Life expectancy (year)")
abline(v=which.min(Japan)+1930-1,col="red",lwd=2,lty=2)

## -----------------------------------------------------------------------------
s<-ChangePoint(Japan,v1=40,v2=65,e1=2,e2=2,N=10000,burnin=1000)
hist(s,freq=FALSE,main="histogram of the change point",xlab="")
dens<-density(s,bw = "nrd0",kernel='gaussian',window=kernel)
plot(dens,xlab="",main="kernel estimation of the density of change point")
abline(v=dens$x[which.max(dens$y)],col="red")
abline(v=which.min(Japan),col="blue")
results<-matrix(0,nrow=1,ncol=2)
rownames(results)<-"Change point"
colnames(results)<-c("Theoretical","Computed")
results[,1]<-which.min(Japan)-1+1930
results[,2]<-round(dens$x[which.max(dens$y)]-1+1930)

library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
Checkfunction1(Japan,2,1)+1930-1

## -----------------------------------------------------------------------------
China<-c(0.1422 ,0.1792 ,0.2197 ,0.2255 ,0.2668 ,0.3112 ,0.3473 ,0.4056 ,0.8193 ,1.1076 ,1.2096 ,0.8505 ,0.6738 ,0.6579 ,0.6504 ,0.6908 ,0.7401 ,0.6035 ,0.6320 ,0.7540 ,0.9763 ,1.0712 ,1.1140 ,1.1336 ,1.1317 ,1.2774 ,1.3010 ,1.3998 ,1.5355 ,1.5471 ,1.4944 ,1.4561 ,1.5629 ,1.6253 ,1.7434 ,1.8577 ,1.9216 ,2.0229 ,2.1322 ,2.1269 ,2.1114 ,2.1847 ,2.2630 ,2.3971 ,2.5203 ,2.7060 ,2.7989 ,2.7814 ,2.6419 ,2.6144 ,2.6649 ,2.7056 ,2.9622 ,3.4555 ,3.9482 ,4.4159 ,4.8481 ,5.1847 ,5.5386 ,5.7939 ,6.2950 ,6.9223 ,7.0623 ,7.1506 ,7.1353 ,7.0003 ,6.8742 ,6.9812 ,7.2077 ,7.3163 ,7.4117)
plot(x=1950:2020,y=China,main="The overview of Annual CO2 emissions per capita of China from 1950 to 2020",xlab="",ylab="Annual CO2 emissions per capita (ton)")
abline(v=2001,col="red",lwd=2,lty=2)

## -----------------------------------------------------------------------------
#Changepoint method
s<-ChangePoint(China,v1=1,v2=5,e1=1,e2=2,N=10000,burnin=1000)
hist(s,freq=FALSE,main="histogram of the change point",xlab="")
dens<-density(s,bw = "nrd0",kernel='gaussian',window=kernel)
plot(dens,main="kerne estimation of the density of changepoint")
abline(v=dens$x[which.max(dens$y)],col="red")
abline(v=50,col="blue")
results<-matrix(0,nrow=1,ncol=2)
rownames(results)<-"Change point"
colnames(results)<-c("Theoretical","Computed")
results[,1]<-2001
results[,2]<-round(dens$x[which.max(dens$y)]+1950)

library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
#Checkfunction2 method
micro<-Checkfunction2(China,3,6)+1950-1
macro<-Checkfunction2(China,10,6)+1950-1
results<-matrix(0,nrow=2,ncol=6)
rownames(results)<-c("micro perspective","macro perspective")
results[1,]<-micro
results[2,]<-macro

library(knitr)
knitr::kable(results,align=rep('c', 5))

## -----------------------------------------------------------------------------
plot(x=1950:2020,y=China,main="The most possible 6 change points of CO2 emissions of China",xlab="",ylab="Annual CO2 emissions per capita (ton)")
abline(v=micro,col="green",lwd=1.5,lty=2)

## -----------------------------------------------------------------------------
plot(x=1950:2020,y=China,main="The most possible 6 change points of CO2 emissions of China",xlab="",ylab="Annual CO2 emissions per capita (ton)")
abline(v=macro,col="blue",lwd=1.5,lty=2)

