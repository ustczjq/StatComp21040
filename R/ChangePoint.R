#' @title A illustration dataset
#' @name China
#' @description A dataset used for change point detection.
#' @examples
#' \dontrun{
#' data(China)
#' }
NULL

#' @title A illustration dataset
#' @name Japan
#' @description A dataset used for change point detection.
#' @examples
#' \dontrun{
#' data(Japan)
#' }
NULL

#' @title The Detection of Change Point using a MCMC method
#' @description returns the generated vector of possible change point,to see the most likely one, please use a histogram or a kernel density estimation.
#' @param y a vector that contains the interested value in time series data.
#' @param v1 the mean of the prior distribution of a Normal distribution(mu1,sigma^2)'s mu1.i.e. mu1 ~N(v1,e1^2).
#' @param v2 the mean of the prior distribution of a Normal distribution(mu2,sigma^2)'s mu2.i.e. mu2 ~N(v2,e2^2).
#' @param e1 the standard error of the prior distribution of a Normal distribution(mu1,sigma^2)'s mu1.i.e. mu1 ~N(v1,e1^2).
#' @param e2 the standard error of the prior distribution of a Normal distribution(mu2,sigma^2)'s mu2.i.e. mu2 ~N(v2,e2^2).
#' @param N the length of chain.
#' @param burnin the burn-in time of the chain.
#' @return a generated sample of possible change point
#' @importFrom stats rgamma rnorm runif
#' @examples
#' \dontrun{
#' result<-ChangePoint(data,v1=5,v2=10,e1=1,e2=1,N=10000,burnin=1000)
#' hist(result,freq=FALSE,main="histogram of the change point",xlab="")
#' dens<-density(result,bw = "nrd0",kernel='gaussian',window=kernel)
#' plot(dens,main="kernel estimation of the density of changepoint")
#' abline(v=dens$x[which.max(dens$y)],col="red")
#' }
#' @export
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
