\name{QRb.pred}
\alias{QRb.pred}
\title{ Calculate predicted probabilities with \code{QRb.pred} }
\description{
  \code{QRb.pred} is a function that calculates predicted probabilities for the binary quantile regression
	methods \code{QRb} or \code{QRb.AL}.
}
\usage{
  QRb.pred(QRseq.obj, X, burnin)
}
\arguments{
  \item{QRseq.obj}{ an output object of the \code{QRseq} function with \code{QRb} or \code{QRb.AL} method. At leas 9 estimated quantiles are required. }
	\item{X}{ m x k matrix of predictors. Not that k should equal the number of variables in 'QRb' and 'QRb.AL'. }
	\item{burnin}{ the number of burnin draws that should be discared (default=1, meaning all draws are included). }
}
\value{
  An m x 1 vector containing the predicted probabilities.
}
\author{
  Dries F. Benoit
}
\details{
  \code{QRb.pred} is a function that calculates the predicted probabilities based on a matrix of predictors
	\emph{X} and an object containing the parameter estimates of a sequence of binary quantile regressions
	(preferably estimated with the \code{QRseq} function). The rationale behind the approach is described in 
	Kordas (2006) and applied in Migueis, V.L., Benoit, D.F. and Van den Poel, D. (2012). Note that the more
	quantiles are estimated, the more fine-grained the predicted probabilities will be. 
}
\references{
  Kordas, G. (2006). Binary regression quantiles, \emph{Journal of Applied Econometrics}, \bold{21}(3), 387-407.
	\cr
	\cr
	Migueis, V.L., Benoit, D.F. and Van den Poel, D. (2012). Enhanced decision support in credit scoring using 
	Bayesian binary quantile regression, \emph{Journal of the Operational Research Society}, (\emph{in press}).
}
\examples{
# Simulate data from binary regression model
set.seed(1234)
n <- 200
X <- cbind(1,matrix(runif(n=n*2, min=-5, max=5),ncol=2))
ystar <- X\%*\% c(1,1.5,-.5) + rnorm(n=n, mean=0, sd=abs(2*X[,1]))
y <- as.numeric(ystar>0)

## Set input arguments for a sequence of QRb models
Data <- list(y=y, X=X, p=seq(.1,.9,.1))
Mcmc <- list(R=5000, keep=1)

out <- QRseq(Method="QRb", Data=Data, Mcmc=Mcmc)

# Calculate predicted probabilities
pred <- QRb.pred(QRseq.obj=out, X=Data$X, burnin=2000)

# Make histogram of predicted probabilities 
hist(pred,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(Data$y==as.numeric(pred>.5))

# Compare with logit model
mylogit <- glm(Data$y ~ 0 + Data$X, family=binomial(logit))

# Make histogram of predicted probabilities 
hist(mylogit$fit,breaks=10)

# Calculate Percentage Correclty Classified (PCC)
mean(Data$y==as.numeric(mylogit$fit>.5))
}
