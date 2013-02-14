 
\name{prior}
\alias{prior}
\title{ Create prior for Bayesian quantile regression }
\description{
  \code{prior} creates an S3 object of class \code{bayesQR.prior} that contains all necessary prior information to estimate a Bayesian quantile regression model.
}
\usage{
  prior(formula, data, alasso, ...)
}
\arguments{
  \item{formula}{ a symbolic description of the model to be fit. }
  \item{data}{ an optional data frame containing the variables in the model. }
	\item{alasso}{ logical flag for adaptive lasso variable selection (default=FALSE). }
	\item{...}{ the prior parameters that are dependent on the method that is used. If omitted, a standard diffuse prior will be used (see details section). }
}
\value{
  An object of class \code{bayesQR}, basically a list containing the following elements:
  \item{method}{ a string containing the method that was used, i.e. indicating whether the dependent variable was continuous or binary and whether adaptive lasso variable selection was used. }
	\item{...}{ the prior parameters that are dependent on the method that is used (see details section). }
}
\details{
  The function \code{prior} builds the prior for four types of models, depending on whether the dependent variable is continuous or binary and whether adaptive lasso variable selection is used. 
	Every non-specified prior parameter will get the default value.

  Continuous dependent variable without adaptive lasso variable selection:
  \itemize{
    \item \bold{beta0} : nvar(X) x 1 vector of prior means (default: 0)
    \item \bold{V0} : nvar(X) x nvar(X) prior covariance matrix (default: 100*diag(ncol(X)))
    \item \bold{shape0} : shape parameter for inverse Gamma prior for sigma (default: 0.01)
    \item \bold{scale0} : scale parameter for inverse Gamma prior for sigma (default: 0.01)
  }

  Continuous dependent variable with adaptive lasso variable selection:
  \itemize{
    \item \bold{sigma_shape} : shape parameter for the inverse gamma prior on sigma (default: 0.01)
    \item \bold{sigma_scale} : scale parameter for the inverse gamma prior on sigma (default: 0.01)
    \item \bold{etasq_shape} : shape parameter for the gamma prior on etasq (default: 0.01)
    \item \bold{etasq_scale} : scale parameter for the gamma prior on etasq (default: 0.01)
  }

  Binary dependent variable without adaptive lasso variable selection:
  \itemize{
    \item \bold{beta0} : nvar(X) x 1 vector of prior means (default: 0)
    \item \bold{V0} : nvar(X) x nvar(X) prior covariance matrix (default: 100*diag(ncol(X)))
  }

  Binary dependent variable with adaptive lasso variable selection:
  \itemize{
    \item \bold{lambdasq_shape} : shape parameter for the gamma prior on lambdasq (default: 0.01)
    \item \bold{lambdasq_scale} : scale parameter for the gamma prior on lambdasq (default: 0.01)
  }
}
\author{
  Dries F. Benoit
}
\examples{
# Load the Prostate cancer dataset
data(Prostate)

# Create informative prior object
prior <- prior(lpsa~., data=Prostate, beta0=rep(5,9), V0=diag(9)) 

# Investigate structure of bayesQR.prior object
str(prior)

# Estimate the model parameters with informative prior
out <- bayesQR(lpsa~., data=Prostate, prior=prior, ndraw=5000)

# Print results
summary(out)
}
