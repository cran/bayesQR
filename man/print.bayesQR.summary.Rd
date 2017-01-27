\name{print.bayesQR.summary}
\alias{print.bayesQR.summary}
\title{ Prints the contents of \code{bayesQR.summary} object to the console }
\description{
  \code{print.bayesQR.summary} is an S3 method that prints the content of an S3 object of class \code{bayesQR.summary} to the console. 
}
\usage{
  \method{print}{bayesQR.summary}(x, digits=3, ...)
}
\arguments{
  \item{x}{ an output object of the \code{summary.bayesQR} function, i.e. an S3 object of class \code{bayesQR.summary}. }
	\item{digits}{ the number of printed digits of the estimates (default=3). }
  \item{...}{ additional arguments that are passed to the generic \code{print} function }
}
\value{
  Formatted output of a \code{bayesQR.summary} object.
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

# Analyze 0.5 quantile using default prior and adaptive lasso
out <- bayesQR(y ~ X, alasso=TRUE, ndraw=5000) 

# Return Bayes estimates and credible intervals 
sum <- summary(out, burnin=1000)

# Print the bayesQR.summary object
sum
}
