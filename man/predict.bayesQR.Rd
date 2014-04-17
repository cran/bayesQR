\name{predict.bayesQR}
\alias{predict.bayesQR}
\title{ Calculate predicted probabilities for binary quantile regression }
\description{
  \code{predict.bayesQR} is an S3 method that calculates predicted probabilities for binary quantile regression, both with or without adaptive lasso.
}
\usage{
  \method{predict}{bayesQR}(object, X, burnin=0, ...)
}
\arguments{
  \item{object}{ an output object of the \code{bayesQR} function, S3 class bayesQR, with binary dependent variable, with or without adaptive lasso. At leas 9 estimated quantiles are required. }
  \item{X}{ matrix of predictors (should be of the same type as the variables used to estimate the model). }
  \item{burnin}{ the number of burnin draws that should be discared (default=0, meaning all draws are included). }
  \item{...}{ additional parameters passed to the generic \code{predict} function. }
}
\value{
  A vector containing the predicted probabilities.
}
\author{
  Dries F. Benoit
}
\details{
  \code{predict.bayesQR} is an S3 method that calculates the predicted probabilities based on a matrix of predictors \emph{X} and an object containing the parameter estimates of a sequence of binary quantile regressions.
  The rationale behind the approach is described in Kordas (2006) and applied in Migueis, V.L., Benoit, D.F. and Van den Poel, D. (2012).
  Note that the more quantiles are estimated, the more fine-grained the predicted probabilities will be. 
}
\references{
  Kordas, G. (2006). Binary regression quantiles, \emph{Journal of Applied Econometrics}, \bold{21}(3), 387-407.
  \cr
  \cr
  Migueis, V.L., Benoit, D.F. and Van den Poel, D. (2012). Enhanced decision support in credit scoring using Bayesian binary quantile regression, \emph{Journal of the Operational Research Society}, (\emph{in press}).
}
\examples{
# Simulate data from binary regression model
set.seed(1234)
n <- 200
X <- matrix(runif(n=n*2, min=-5, max=5),ncol=2)
ystar <- cbind(1,X)\%*\% c(1,1.5,-.5) + rnorm(n=n, mean=0, sd=abs(2*X[,1]))
y <- as.numeric(ystar>0)

## Estimate a sequence of binary quantile regression models
out <- bayesQR(y ~ X, quantile=seq(.1,.9,.1), ndraw=4000)

# Calculate predicted probabilities
pred <- predict(object=out, X=cbind(1,X), burnin=2000)

# Make histogram of predicted probabilities 
hist(pred,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(pred>.5))

# Compare with logit model
mylogit <- glm(y ~ X, family=binomial(logit))

# Make histogram of predicted probabilities 
hist(mylogit$fit,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(y==as.numeric(mylogit$fit>.5))
}
