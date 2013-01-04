\name{QRc}
\alias{QRc}
\title{ Bayesian quantile regression with continuous dependent variable }
\description{
  \code{QRc} implements a Bayesian method for estimating quantile regression parameters. This
  method uses the asymmetric Laplace distribution for this goal (see references). To improve 
  the speed of the routine, the Markov Chain Monte Carlo (MCMC) part of the algorithm is 
  programmed in Fortran and is called from within the R function \code{QRc}.
}
\usage{
  QRc(Data, Prior, Mcmc)
}
\arguments{
  \item{Data}{ list(y, X, p) }
  \item{Prior}{ list(beta0, V0, shape0, scale0) }
  \item{Mcmc}{ list(R, keep) }
}
\value{
  A list containing:
  \item{method}{ a string containing 'QRc' }
  \item{p}{ the quantile that was estimated }
  \item{betadraw}{ R/keep x nvar(X) matrix of beta draws }
  \item{sigmadraw}{ R/keep vector of sigma draws }
}
\details{
  \tabular{ll}{
  Model: \tab y = X\emph{beta} + e \cr
         \tab \emph{e ~ ALD(location=0, scale=sigma, quantile=p)} \cr
  Priors: \tab \emph{beta ~ N(beta0, V0)} \cr
          \tab \emph{sigma ~ invGamma(shape0,scale0)} \cr
  }
  List arguments contain: \cr
  \itemize{
    \item \bold{X} : n x nvar(X) matrix of predictors (first column should be a vector of ones if intercept is desired)
    \item \bold{y} : n x 1 vector of observations (dependent variable)
    \item \bold{p} : quantile of interest (should be between 0 and 1)
    \item \bold{beta0} : nvar(X) x 1 vector of prior means (default: 0)
    \item \bold{V0} : nvar(X) x nvar(X) prior covariance matrix (default: 100*diag(ncol(X)))
    \item \bold{shape0} : shape parameter for inverse Gamma prior for sigma (default: 0.01)
    \item \bold{scale0} : scale parameter for inverse Gamma prior for sigma (default: 0.01)
    \item \bold{R} : number of MCMC draws
    \item \bold{keep} : thinning parameter, i.e. keep every keepth draw (default: 1)
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{ The algorithm is an improved (Gibbs instead of Metropolis-Hastings) implementation of: \cr
  Yu, K. and Moyeed, R. (2001). Bayesian quantile regression, \emph{Statistics 
  and Probability Letters}, 54, 437-447. \cr
  \cr
  Also see: \cr
  Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression,
  \emph{Journal of Statistical Computation and Simulation}, 81(11), 1565-1578.
}
\examples{
# Simulate data from heteroskedastic regression
n <- 200
X <- runif(n=n,min=0,max=10)
X <- cbind(1,X)
y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6*X[,2])

# Initiate plot
## Plot datapoints
plot(X[,2], y, main="", cex=.6, xlab="X")

# Initialize the inputs for QRc
Data = list(y=y, X=X, p=.05)
Prior = list(beta0=rep(0,ncol(X)),V0=100*diag(ncol(X)), shape0=0.01, scale0=0.01)
Mcmc = list(R=5000, keep=1)

# Write loop to analyze 5 quantiles
for (i in 1:5) {
    if (i==1) p = .05
    if (i==2) p = .25
    if (i==3) p = .50
    if (i==4) p = .75
    if (i==5) p = .95

    Data$p = p
    out = QRc(Data=Data, Prior=Prior, Mcmc=Mcmc)

    ## Add quantile regression lines to the plot (exclude first 500 burn-in draws)
    abline(a=mean(out$betadraw[500:Mcmc$R,1]),b=mean(out$betadraw[500:Mcmc$R,2]),lty=i,col=i)
}

# Estimate and plot OLS model
outOLS = lm(y~X[,2])
abline(outOLS,lty=1,lwd=2,col=6)

# Add legend to plot
legend(x=0,y=max(y),legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")
}
