\name{QRc.AL}
\alias{QRc.AL}
\title{ Bayesian quantile regression with continuous dependent variable using adaptive LASSO }
\description{
  \code{QRc.AL} implements a Bayesian method for variable selection and estimation in quantile
  regression models, as described in Alhamzawi, Yu and Benoit (forthcoming). This method extends 
  the Bayesian Lasso quantile regression (Li et al., 2010) by allowing different penalization 
  parameters for different regression coefficients. Also, all hyperparameters are treated as unknows 
  and are estimated from the data along with other parameters. To improve the speed of the routine, 
  the Markov Chain Monte Carlo (MCMC) part of the algorithm is programmed in FORTRAN and is called 
  from within the R function \code{QRc.AL}.
}
\usage{
  QRc.AL(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(sig_shape, sig_rate) }
  \item{Mcmc}{ list(R, keep, step_delta) }
}
\value{
  A list containing:
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
  \item{lambda2draw}{ R/keep x nvar(X) matrix of lambda2 draws }
  \item{deltadraw}{ R/keep vector of delta draws }
  \item{taudraw}{ R/keep vector of tau draws }
  \item{rejrate}{ Rejection rate of the MH proposals for delta }

}
\details{
  \tabular{ll}{
  Model: \tab y = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=sigma, quantile=p)} \cr
  Priors: \tab \emph{beta ~ ALD(location=0, scale=sigma/lambda, p=0.5)} \cr
          \tab \emph{lambda^2 ~ InvGamma(shape=delta, rate=tau)} \cr
	  \tab \emph{delta} \eqn{\propto} \emph{1} \cr
          \tab \emph{tau ~ Uniform(tau^(-1)).}
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if an intercept is desired)
    \item \bold{y} : n x 1 vector of observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{sig_shape} : shape parameter for the gamma prior on sigma (default: 0.01)
    \item \bold{sig_rate} : rate parameter for the gamma prior on sigma (default: 0.01)
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
    \item \bold{step_delta} : MH step for delta, tune to get \emph{rejrate} between 0.3 and 0.5
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is an implementation (with minor changes) of: \cr
  Alhamzawi, R., Yu, K. and Benoit, D.F. (2011). Bayesian adaptive LASSO 
  quantile regression, \emph{Statistical Modelling}, forthcoming. \cr
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
Prior <- list(sig_shape=.01, sig_rate=.01)
Mcmc <- list(R=10000, keep=1, step_delta=1)

# Estimate the model parameters 
out <- QRc.AL(Data=Data, Prior=Prior, Mcmc=Mcmc)

# Check rejection rate (should be between .3 and .5)
out$rejrate

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
