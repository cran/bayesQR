\name{QRb}
\alias{QRb}
\title{ Bayesian quantile regression with binary dependent variable }
\description{
  \code{QRb} implements a Bayesian method for estimating binary quantile regression 
  parameters. This method uses the asymmetric Laplace distribution for this goal 
  (see references). To improve the speed of the routine, the Markov Chain Monte Carlo 
  (MCMC) part of the algorithm is programmed in Fortran and is called from within the 
  R function \code{QRb}.
}
\usage{
  QRb(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(beta0, V0) }
  \item{Mcmc}{ list(R, keep) }
}
\value{
  A list containing:
  \item{method}{ a string containing 'QRb' }
  \item{p}{ the quantile that was estimated }
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
}
\details{
  \tabular{ll}{
  Model: \tab y* = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=1, quantile=p)} \cr
         \tab \emph{y = 1, if (y* > 0); y = 0, otherwise} \cr
  Priors: \tab \emph{beta ~ N(beta0, V0)} \cr
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if intercept is desired)
    \item \bold{y} : n x 1 vector of zero-one observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{beta0} : nvar(X) x 1 vector of prior means (default: 0)
    \item \bold{V0} : nvar(X) x nvar(X) prior covariance matrix (default: 100*diag(ncol(X)))
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is an improved (Gibbs instead of Metropolis-Hastings) implementation of: \cr
  Benoit, D.F. and Van den Poel, D. (2012). Binary quantile regression: A Bayesian approach
  based on the asymmetric Laplace distribution, \emph{Journal of Applied Econometrics}, 27(7), 1174-1188.
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
Prior = list(beta0=rep(0,ncol(X)),V0=100*diag(ncol(X)))
Mcmc = list(R=5000, keep=1)

out <- QRb(Data=Data, Prior=Prior, Mcmc=Mcmc)

# Traceplot of betadraws
matplot(out$betadraw,typ="l")

# Posterior distribution of betadraws
hist(out$betadraw,breaks=100)

# Posterior Bayes estimates and credible intervals
quantile(out$betadraw[,1],c(.025,.5,.975))
}
