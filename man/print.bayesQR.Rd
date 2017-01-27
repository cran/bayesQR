\name{print.bayesQR}
\alias{print.bayesQR}
\title{ Prints the contents of \code{bayesQR} object to the console }
\description{
  \code{print.bayesQR} is an S3 method that prints the content of an S3 object of class \code{bayesQR} to the console. 
}
\usage{
  \method{print}{bayesQR}(x, digits=3, ...)
}
\arguments{
  \item{x}{ an output object of the \code{bayesQR} function, i.e. an S3 object of class \code{bayesQR}. }
  \item{digits}{ the number of printed digits of the estimates (default=3). }
  \item{...}{ additional arguments that are passed to the generic \code{print} function }
}
\value{
  Formatted output of a \code{bayesQR} object.
}
\author{
  Dries F. Benoit
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(666)
n <- 200
X <- runif(n=n,min=0,max=10)
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Analyze 0.5 quantile using default prior
out <- bayesQR(y ~ X, ndraw=5000) 

# Print the bayesQR object
print(out)
}
