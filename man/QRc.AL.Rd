\name{QRc.AL}
\alias{QRc.AL}
\title{ Bayesian quantile regression with continuous dependent variable using adaptive LASSO }
\description{
  \code{QRc.AL} implements a Bayesian method for variable selection and estimation in quantile
  regression models. The implementation is based on Alhamzawi, Yu and Benoit (2012). This method extends 
  the Bayesian Lasso quantile regression (Li et al., 2010) by allowing different penalization 
  parameters for different regression coefficients. To improve the speed of the routine, 
  the Markov Chain Monte Carlo (MCMC) part of the algorithm is programmed in Fortran and is called 
  from within the R function \code{QRc.AL}.
}
\usage{
  QRc.AL(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(sigma_shape, sigma_scale, etasq_shape, etasq_scale) }
  \item{Mcmc}{ list(R, keep) }
}
\value{
  A list containing:
  \item{method}{ a string containing 'QRc.AL' }
  \item{p}{ the quantile that was estimated }
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
}
\details{
  \tabular{ll}{
  Model: \tab y = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=sigma, quantile=p)} \cr
  Priors: \tab \emph{beta ~ ALD(location=0, scale=sigma/lambda, p=0.5)} \cr
          \tab \emph{sigma ~ InvGamma(shape=sigma_shape, scale=sigma_scale)} \cr
          \tab \emph{(lambda/sigma)^2 ~ Gamma(shape=etasq_shape, scale=etasq_scale)} 
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if an intercept is desired)
    \item \bold{y} : n x 1 vector of observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{sigma_shape} : shape parameter for the inverse gamma prior on sigma (default: 0.01)
    \item \bold{sigma_scale} : scale parameter for the inverse gamma prior on sigma (default: 0.01)
    \item \bold{etasq_shape} : shape parameter for the gamma prior on etasq (default: 0.01)
    \item \bold{etasq_scale} : scale parameter for the gamma prior on etasq (default: 0.01)
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is an implementation (with minor changes) of: \cr
  Alhamzawi, R., Yu, K. and Benoit, D.F. (2012). Bayesian adaptive LASSO 
  quantile regression, \emph{Statistical Modelling}, 12(3), 279-297. \cr
  \cr
  Also see: \cr
  Li, Q., Xi, R. and Lin, N. (2010). Bayesian Regularized Quantile Regression. 
  \emph{Bayesian Analysis}, 5, 1-24.
}
\examples{
# Load the Prostate cancer dataset
data(Prostate)

# Initialize the inputs for QRc.AL
Data <- list(y=as.vector(Prostate[,9]),X=as.matrix(Prostate[,1:8]),p=.95)
Prior <- list(sigma_shape=.01, sigma_scale=.01, etasq_shape=.01, etasq_scale=.01)
Mcmc <- list(R=5000, keep=1)

# Estimate the model parameters 
out <- QRc.AL(Data=Data, Prior=Prior, Mcmc=Mcmc)

# Plot the trace plots of the regression parameters
for (i in 1:ncol(Data$X)){
  plot(out$betadraw[,i],typ="l")
  if (i == 1) par(ask=TRUE)
  if (i == ncol(Data$X)) par(ask=FALSE)
}

# Bayes estimates and credible intervals for the regression parameters
burnin <- 2000
parest <- matrix(NA,ncol=3,nrow=ncol(Data$X))
for (i in 1:ncol(Data$X)){
  parest[i,] <- round(quantile(out$betadraw[burnin:Mcmc$R,i],c(.025,.5,.975)),5)
}
rownames(parest) <- names(Prostate)[1:8]
colnames(parest) <- c(2.5,50,97.5)
parest
}
