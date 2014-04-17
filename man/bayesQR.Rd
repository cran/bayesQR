\name{bayesQR}
\alias{bayesQR}
\title{ Bayesian quantile regression }
\description{
  \code{bayesQR} implements a Bayesian method for estimating quantile regression  models (see references).
	To improve the speed of the routine, the Markov Chain Monte Carlo (MCMC) part of the algorithm is programmed in Fortran and is called from within the R function \code{bayesQR}.
}
\usage{
  bayesQR(formula, data, quantile, alasso, ndraw, keep, prior)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ an optional data frame containing the variables in the model. }
  \item{quantile}{ numerical scalar or vector containing quantile(s) of interest (default=0.5). }
	\item{alasso}{ logical flag for adaptive lasso variable selection (default=FALSE). }
	\item{ndraw}{ number of MCMC draws. }
	\item{keep}{ thinning parameter, i.e. keep every keepth draw (default=1). }
	\item{prior}{ an S3 object of class \code{"prior"}. If omitted, a diffuse prior will be assumed (see \code{\link{prior}}). }
}
\value{
  An object of class \code{bayesQR}, basically a list of lists.
	For every estimated quantile a list is created containing the following elements:
  \item{method}{ a string containing the method that was used, i.e. indicating whether the dependent variable was continuous or binary and whether adaptive lasso variable selection was used. }
  \item{quantile}{ the quantile that was estimated. }
	\item{names}{ character vector containing the names of the independent variables in the model. }
  \item{betadraw}{ ndraw/keep x nvar(X) matrix of beta draws. }
  \item{sigmadraw}{ ndraw/keep vector of sigma draws (only in case of continuous dependent variable). }
}
\details{
  The function \code{bayesQR} can estimate four types of models, depending on whether the dependent variable is continuous or binary and whether adaptive lasso variable selection is used.
  \itemize{
    \item Continuous dependent variable without adaptive lasso variable selection:
  	\tabular{ll}{
    Model:  \tab y = X\emph{beta} + e \cr
            \tab \emph{e ~ ALD(location=0, scale=sigma, quantile=p)} \cr
    Priors: \tab \emph{beta ~ N(beta0, V0)} \cr
            \tab \emph{sigma ~ invGamma(shape0,scale0)} \cr
  	}
  
    \item Continuous dependent variable with adaptive lasso variable selection:
  	\tabular{ll}{
    Model:  \tab y = X\emph{beta} + e \cr
            \tab \emph{e ~ ALD(location=0, scale=sigma, quantile=p)} \cr
    Priors: \tab \emph{beta ~ ALD(location=0, scale=sigma/lambda, p=0.5)} \cr
            \tab \emph{sigma ~ InvGamma(shape=sigma_shape, scale=sigma_scale)} \cr
            \tab \emph{(lambda/sigma)^2 ~ Gamma(shape=etasq_shape, scale=etasq_scale)} \cr
  	}
  
    \item Binary dependent variable without adaptive lasso variable selection:
    \tabular{ll}{
    Model:  \tab y* = X\emph{beta} + e \cr
            \tab \emph{e ~ ALD(location=0, scale=1, quantile=p)} \cr
            \tab \emph{y = 1, if (y* > 0); y = 0, otherwise} \cr
    Priors: \tab \emph{beta ~ N(beta0, V0)} \cr
    }
  
    \item Binary dependent variable with adaptive lasso variable selection:
    \tabular{ll}{
    Model:  \tab y* = X\emph{beta} + e \cr
            \tab \emph{e ~ ALD(location=0, scale=1, quantile=p)} \cr
            \tab \emph{y = 1, if (y* > 0); y = 0, otherwise} \cr
    Priors: \tab \emph{beta ~ ALD(location=0, scale=1/lambda, p=0.5)} \cr
            \tab \emph{lambda^2 ~ Gamma(shape=lambdasq_shape, scale=lambdasq_scale)} \cr
    }
  }
}
\author{
  Dries F. Benoit, Rahim Al-Hamzawi, Keming Yu and Dirk Van den Poel
}
\references{
  \itemize{  
  \item Continuous dependent variable without adaptive lasso variable selection:

  Yu, K. and Moyeed, R. (2001). Bayesian quantile regression, \emph{Statistics and Probability Letters}, 54, 437-447.

  Also see: \cr
  Kozumi, H. and Kobayashi, G. (2011). Gibbs sampling methods for Bayesian quantile regression, \emph{Journal of Statistical Computation and Simulation}, 81(11), 1565-1578. \cr

  \item Continuous dependent variable with adaptive lasso variable selection:

  Alhamzawi, R., Yu, K. and Benoit, D.F. (2012). Bayesian adaptive LASSO quantile regression, \emph{Statistical Modelling}, 12(3), 279-297.

  Also see: \cr
  Li, Q., Xi, R. and Lin, N. (2010). Bayesian Regularized Quantile Regression. \emph{Bayesian Analysis}, 5, 1-24. \cr

  \item Binary dependent variable without adaptive lasso variable selection:

  Benoit, D.F. and Van den Poel, D. (2012). Binary quantile regression: A Bayesian approach based on the asymmetric Laplace distribution, \emph{Journal of Applied Econometrics}, 27(7), 1174-1188. \cr

  \item Binary dependent variable with adaptive lasso variable selection:

  Al-Hamzawi, R., Benoit, D.F. and Yu, K. (2012). Binary quantile regression with adaptive lasso: A Bayesian approach, working paper.
	}
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Estimate series of quantile regressions with adaptive lasso
out <- bayesQR(y~X, quantile=c(.05,.25,.5,.75,.95), alasso=TRUE, ndraw=5000)

# Initiate plot
## Plot datapoints
plot(X, y, main="", cex=.6, xlab="X")

## Add quantile regression lines to the plot (exclude first 500 burn-in draws)
sum <- summary(out, burnin=500)
for (i in 1:length(sum)){
  abline(a=sum[[i]]$betadraw[1,1],b=sum[[i]]$betadraw[2,1],lty=i,col=i)
}

# Estimate and plot OLS model
outOLS = lm(y~X)
abline(outOLS,lty=1,lwd=2,col=6)

# Add legend to plot
legend(x=0,y=max(y),legend=c(.05,.25,.50,.75,.95,"OLS"),lty=c(1,2,3,4,5,1),
       lwd=c(1,1,1,1,1,2),col=c(1:6),title="Quantile")
}
