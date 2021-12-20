#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp for homework of Chapter 9 
//' @param a the parameter of the given density
//' @param b the parameter of the given density
//' @param n the parameter of the given density
//' @param N the number of samples
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' GC<- GibbsC(1,1,25,10000)
//' plot(GC[,1],main= paste("a=",a,",b=",b,",n=",n,",Condition function:Binomial"))
//' plot(GC[,2],main= paste("a=",a,",b=",b,",n=",n,",Condition function:Beta"))
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix GibbsC (int a,int b,int n,int N){
  NumericMatrix X(N,2); 
  X(0,0)=0;
  X(0,1)=0;
  double x1=0;
  double x2=0;
  for (int i=1;i<N;i++){ 
    x2=X(i-1,1);
    X(i,0)=rbinom(1,n,x2)[0];
    x1=X(i,0);
    X(i,1)=rbeta(1,x1+a,n-x1+b)[0];
  }
  return X;
}