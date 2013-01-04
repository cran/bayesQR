\name{QRseq}
\alias{QRseq}
\title{ Estimate series of Bayesian quantile regressions with \code{QRc} or \code{QRb} }
\description{
  \code{QRseq} is a wrapper function that automates the estimation of a series of quantile
	regressions. This function will be used most often to produce the input for the \code{QRplot} 
	or \code{QRb.pred} functions. It should be noted that, as always, the convergence of each mcmc chain 
	has to be checked by the user. A visual convergence check can be done, for example, with the 
	\code{QRsummary} function.
}
\usage{
  QRseq(Method, Data, Prior, Mcmc)
}
\arguments{
  \item{Method}{ should be \code{QRc}, \code{QRc.AL}, \code{QRb} or \code{QRb.AL} }
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ See documentation of specific method used (i.e. \code{QRc}, \code{QRc.AL}, \code{QRb} or \code{QRb.AL}) }
  \item{Mcmc}{ list(R, keep) }
}
\value{
  A list of lists, where each lower level list contains:
  \item{method}{ the method that was used, i.e. \code{QRc}, \code{QRc.AL}, \code{QRb} or \code{QRb.AL} }
  \item{p}{ a scalar containing the estimated quantile }
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
	And in case of \code{QRc} and \code{QRc.AL} also:
  \item{sigmadraw}{ R/keep vector of sigma draws }
}
\details{
  For more details about the models, see the specific documentation about \code{QRc}, \code{QRc.AL}, \code{QRb} or \code{QRb.AL}.
}
\author{
  Dries F. Benoit
}
\examples{
# Simulate data from heteroskedastic regression
n <- 200
X <- runif(n=n,min=0,max=10)
X <- cbind(1,X)
y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6*X[,2])

# Initialize the inputs for QRc
Data = list(y=y, X=X, p=c(.25,.5,.75))
Mcmc = list(R=5000)

# Analyze 5 quantiles using default prior
out = QRseq(Method="QRc", Data=Data, Mcmc=Mcmc)

# Look at the structure of the output
str(out)
}
