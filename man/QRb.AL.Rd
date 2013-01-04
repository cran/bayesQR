\name{QRb.AL}
\alias{QRb.AL}
\title{ Bayesian quantile regression with binary dependent variable using adaptive LASSO }
\description{
  \code{QRb.AL} implements a Bayesian method for estimating binary quantile regression 
  parameters with adaptive LASSO penalty (see references). To improve the speed of the 
	routine, the Markov Chain Monte Carlo (MCMC) part of the algorithm is programmed in 
	Fortran and is called from within the R function \code{QRb.AL}.
}
\usage{
  QRb.AL(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(lambdasq_shape, lambdasq_scale) }
  \item{Mcmc}{ list(R, keep) }
}
\value{
  A list containing:
  \item{method}{ a string containing 'QRb.AL' }
  \item{p}{ the quantile that was estimated }
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
}
\details{
  \tabular{ll}{
  Model: \tab y* = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=1, quantile=p)} \cr
         \tab \emph{y = 1, if (y* > 0); y = 0, otherwise} \cr
  Priors: \tab \emph{beta ~ ALD(location=0, scale=1/lambda, p=0.5)} \cr
          \tab \emph{lambda^2 ~ Gamma(shape=lambdasq_shape, scale=lambdasq_scale)}
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if intercept is desired)
    \item \bold{y} : n x 1 vector of zero-one observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{lambdasq_shape} : shape parameter for the gamma prior on lambdasq (default: 0.01)
    \item \bold{lambdasq_scale} : scale parameter for the gamma prior on lambdasq (default: 0.01)
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is a slightly modified implementation of: \cr
  Al-Hamzawi, R., Benoit, D.F. and Yu, K. (2012). Binary quantile regression with adaptive lasso:
	A Bayesian approach, working paper.
}
\examples{
## set seed value
set.seed(66)

## set dimension of problem
n <- 200
k <- 4

## set parameter values
beta <- c(-1,-.5,0,2.5)

## simulate Laplace distributed data
X <- cbind(1,matrix(runif(n=n*(k-1)),nrow=n))
y <- X\%*\%beta
Vo1 <- rgamma(n=n,shape=1,scale=2)
ystar <- y + sqrt(Vo1)*rnorm(n=n)
y <- ifelse(ystar>0,1,0)
table(y)

## prepare the inputs
Data <- list(y=y,X=X,p=.5)
Prior <- list(lambda_shape=.01,lambda_scale=.01)
Mcmc <- list(R=5000,keep=1)

## estimate the parameters
out <- QRb.AL(Data=Data,Prior=Prior,Mcmc=Mcmc)

## calculate posterior median and credible intervals
apply(out$betadraw,FUN=quantile,MARGIN=2,prob=c(.025,.5,.975))
}
