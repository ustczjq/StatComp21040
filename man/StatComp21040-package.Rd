\name{StatComp21040-package}
\alias{ChangePointDetection-package}
\docType{package}
\title{
  The Detection of Change Point with a MCMC method and an intuitive method
}
\description{
  In the study involved with time series,the sudden and even anomalous change of data is a common phenomenon, which sometimes implies valuable information and therefore cannot be simply regarded as a "outlier". This  calls for the change point analysis. 
  Three functions are included and three numeric experiments are implemented in this package to demonstrate two methods for the detection of change point.
  One is an intuitive method using slope, designed to deal with different types of change points,such as "piece-wise constant","V" type and "cliff" type. 
  And the other method is implemented using the MCMC method,where the data is assumed to follow a Normal distribution with the prior assumption that the location parameters also follow Normal distributions. Under this framework,the posterior distributions of necessary parameters are determined and the parameters are derived by Gibbs Sampler,except the sample of change point, which is finished according to a certain acceptance-criterion of Metropolis-Hastings Sampler. Specific steps for these will be detailed in the introduction. 
  What's more, previous codes of homework in StatComp course and several Rcpp functions are included in this package.
}
\details{
  In brief, the function ChangePoint returns the generated vector of possible change point after the burn-in time. To see the most likely one, please use a histogram or a kernel density estimation. See the details in introduction. The function Checkfunction1 and Checkfunction2 return the possible change point.
  To be more specific,the function ChangePoint involves with the Metropolis-Hastings Sampler within the Gibbs Sampler.The only parameter generated using the Metropolis-Hastings algorithm is the change point, as its conditional distribution is too complex to sample. The deduction of posterior distribution of parameters and the criterion of acceptance are demonstrated in the introduction, too.
  The function Checkfunction1 can be used as a preliminary detection for possible change points. This function is inspired by the intuitive idea that the slope of change point and its neighbor points should be large. To simulate the different degrees of impact, I adopt the idea of kernel density estimation by computing a weighted mean of slopes. This is an easy and feasible method for many cases, like a prominent "V" type change point, the piece-wise constant case with a sudden change and so on. 
  However,for some certain types,for example: a sudden surge/plunge after a stable plateau (like a cliff), function Checkfunction1 fails to detect the most obvious change point. The reason may be that the actual change point' slope might be smaller than the non-change ones' slopes. Therefore, I further improve the original algorithm using the ratio of two near slopes to magnify the change, which are computed with k nearest points. And this manifests priority in detecting possible change point in both macro and micro perspective (i.e the span of interested period, by adjusting the parameter k, which means the number of k nearest points).
  The numeric study consists of one simulation test and two real-data analysis, which reveal the excellence of all three functions in detection of change point. Notably, the ChangePoint is a universal method, Checkfunction1 is easy to implement and Checkfunction2 can detect possible change points in different perspectives. Last but not least, Checkfunction1 and Checkfunction2 can act as auxiliary tools for improvement of ChangePoint. 
}
\author{
Jiaqi Zhang <ustczjq1998@ustc.edu.cn>
Maintainer: Jiaqi Zhang <ustczjq1998@ustc.edu.cn>
}
\references{
 APPLICATION OF MCMC TO CHANGE POINT DETECTION by Jaromir Antoch, Praha, David Legat, Praha
 The real data (China and Japan) comes from "Our world in data"
}
\keyword{ package }
\examples{
  \dontrun{
  s<-ChangePoint(y,v1=5,v2=10,e1=1,e2=1,N=10000,burnin=1000)
  hist(s,freq=FALSE,main="histogram of the changepoint",xlab="")
  dens<-density(s,bw = "nrd0",kernel='gaussian',window=kernel)
  plot(dens,main="kernel estimation of the density of changepoint")
  abline(v=dens$x[which.max(dens$y)],col="red")
  Checkfunction1(Japan,2,1)
  Checkfunction2(China,10,6)
  }
}
