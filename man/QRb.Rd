\name{QRb}
\alias{QRb}
\title{ Bayesian quantile regression with binary dependent variable }
\description{
  \code{QRb} implements a Bayesian method for estimating binary quantile regression 
  parameters. This method uses the asymmetric Laplace distribution for this goal 
  (see references). To improve the speed of the routine, the Markov Chain Monte Carlo 
  (MCMC) part of the algorithm is programmed in FORTRAN and is called from within the 
  R function \code{QRb}.
}
\usage{
  QRb(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(betabar, A) }
  \item{Mcmc}{ list(R, keep, step) }
}
\value{
  A list containing:
  \item{betadraw}{ R x nvar(X) array of beta draws }
  \item{loglike}{ R vector of evaluations of the loglikelihood }
  \item{rejrate}{ Rejection rate of the MH proposals for beta }
}
\details{
  \tabular{ll}{
  Model: \tab y* = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=1, quantile=p)} \cr
         \tab \emph{y = 1, if (y* > 0); y = 0, otherwise} \cr
  Priors: \tab \emph{beta ~ N(betabar, A^(-1))} \cr
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if intercept is desired)
    \item \bold{y} : n x 1 vector of zero-one observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{betabar} : nvar(X) x 1 vector of prior means (default: 0)
    \item \bold{A} : nvar(X) x nvar(X) prior precision matrix (inverse variancecovariance matrix) (default: .01*diag(ncol(X)))
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
    \item \bold{step} : MH step for beta, tune to get \emph{rejrate} between 0.3 and 0.5
  }
}
\author{
  Dries F. Benoit, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is an implementation of: \cr
  Benoit, D.F and Van den Poel, D. (2011). Binary quantile regression: A Bayesian approach
  based on the asymmetric Laplace distribution, \emph{Journal of Applied Econometrics}, (in press).
}
\examples{
## Simulate data from binary regression model
set.seed(66)
n = 200
X = as.matrix(runif(n=n, min=0, max=10))
ystar = 1.5*X + rnorm(n=n, mean=0, sd=2*X)
y <- as.numeric(ystar>0)

## Set input arguments for QRb, with quantile of interest p=0.75
Data = list(y=y, X=X, p=.75)
Prior = list(betabar=c(rep(0,ncol(X))),A=.01*diag(ncol(X)))
Mcmc = list(R=100000, keep=10, step=.2)

out <- QRb(Data=Data, Prior=Prior, Mcmc=Mcmc)

# Check rejection rate of MH algorithm
out$rejrate

# Traceplot of betadraws
matplot(out$betadraw,typ="l")

# Posterior distribution of betadraws
hist(out$betadraw,breaks=100)

# Posterior Bayes estimates and credible intervals
quantile(out$betadraw[,1],c(.025,.5,.975))
}
