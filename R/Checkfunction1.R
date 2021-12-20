#' @title The Detection of Change Point using an intuitive method
#' @description returns the the possible change point.
#' @param y a vector that contains the interested value in time series data.
#' @param k the k nearest point involved with the computation of slope.
#' @param p the mean of the prior distribution of a Normal distribution(mu2,sigma^2)'s mu2.i.e. mu2 ~N(v2,e2^2).
#' @return possible change point 
#' @importFrom stats weighted.mean
#' @examples
#' \dontrun{
#' Checkfunction1(data,2,1)
#' }
#' @export
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
