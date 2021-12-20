#' @title The Detection of Change Point using an intuitive method
#' @description returns the the possible change point.
#' @param y a vector that contains the interested value in time series data.
#' @param k the k nearest point involved with the computation of slope.
#' @param p the mean of the prior distribution of a Normal distribution(mu2,sigma^2)'s mu2.i.e. mu2 ~N(v2,e2^2).
#' @return possible change point 
#' @importFrom stats median
#' @examples
#' \dontrun{
#' Checkfunction2(data,7,6)
#' }
#' @export
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