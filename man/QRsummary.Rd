\name{QRsummary}
\alias{QRsummary}
\title{ Summarize the betadraws of \code{QRc}, \code{QRc.AL}, \code{QRb}, \code{QRb.AL} or \code{QRseq} }
\description{
  \code{QRsummary} is a function that summarizes the betadraws of \code{QRc}, \code{QRc.AL}, \code{QRb}, 
  \code{QRb.AL} or \code{QRseq}. For every estimated beta, the Bayes estimate and the posterior credible interval
  is returned. The function also has the option to produce the traceplots of every mcmc chain of the 
  different betas. The latter makes is possible to do a visual check of the convergence of the mcmc chains.
}
\usage{
  QRsummary(QRobj, burnin, credint, traceplot)
}
\arguments{
  \item{QRobj}{ an output object of the \code{QRc}, \code{QRc.AL}, \code{QRb}, \code{QRb.AL} or \code{QRseq} function }
  \item{burnin}{ the number of burnin draws that should be discared (default=1, meaning all draws are included) }
  \item{credint}{ the width of the credible interval (default=c(0.025, 0.975)) }
  \item{traceplot}{ logical value, if TRUE, traceplots will be returned (default=FALSE) }
}
\value{
  If traceplot=FALSE then a matrix of Bayes estimates and credible intervals will be returned. 
  If traceplot=TRUE, for every beta a traceplot will be plotted.
}
\author{
  Dries F. Benoit
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- cbind(1,X)
y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6*X[,2])

# Initialize the inputs for QRc
Data = list(y=y, X=X, p=.5)
Mcmc = list(R=5000)

# Analyze 0.5 quantile using default prior
out = QRc(Data=Data, Mcmc=Mcmc)

# Return Bayes estimates and credible intervals 
QRsummary(out, burnin=1000)

# Produce traceplots of every beta
# QRsummary(out, traceplot=TRUE)
}
