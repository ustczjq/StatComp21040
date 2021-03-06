% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{GibbsC}
\alias{GibbsC}
\title{A Gibbs sampler using Rcpp}
\usage{
GibbsC(a, b, n, N)
}
\arguments{
\item{a}{the parameter of the given density}

\item{b}{the parameter of the given density}

\item{n}{the parameter of the given density}

\item{N}{the number of samples}
}
\value{
a random sample of size \code{n}
}
\description{
A Gibbs sampler using Rcpp for homework of Chapter 9
}
\examples{
\dontrun{
GC<- GibbsC(1,1,25,10000)
plot(GC[,1],main= paste("a=",a,",b=",b,",n=",n,",Condition function:Binomial"))
plot(GC[,2],main= paste("a=",a,",b=",b,",n=",n,",Condition function:Beta"))
}
}
