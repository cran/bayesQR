\name{QRplot}
\alias{QRplot}
\title{ Produce quantile plots with \code{QRplot} }
\description{
  \code{QRplot} is a function that produces quantile plots based on the estimates obtained by
	the \code{QRseq} function. Note that the more quantiles \code{p} are estimated with \code{QRseq},
	the more detailed the quantile plot, produced with \code{QRplot}, will be.
}
\usage{
  QRplot(QRseq.obj, var, burnin, credint, main, xlab, ylab, xlim)
}
\arguments{
  \item{QRseq.obj}{ an output object of the \code{QRseq} function }
	\item{var}{ the index of the variable that has to be plotted (default=1) }
	\item{burnin}{ the number of burnin draws that should be discared (default=1, meaning all draws are included) }
	\item{credint}{ the width of the credible interval (default=c(0.025, 0.975)) }
	\item{main}{ Main title of the plot (default="")  }
	\item{xlab}{ Label of the x-axis (default="quantile") }
	\item{ylab}{ Label of the y-axis (default="beta") }
	\item{xlim}{ Plot region of the x-axis (default=c(0,1)) }
}
\value{
  A quantile plot. 
}
\author{
  Dries F. Benoit
}
\details{
  \code{QRplot} is a function that produces quantile plots based on the estimates obtained by
	the \code{QRseq} function. A quantile plot shows how the value of the regression parameter
	changes over a range of quantiles together with the associated credible interval. Note that 
	the more quantiles \code{p} are estimated with \code{QRseq}, the more detailed the quantile 
	plot, produced with \code{QRplot}, will be. The minimum number of quantiles to plot is 2.
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- cbind(1,X)
y <- 1 + 2*X[,2] + rnorm(n=n, mean=0, sd=.6*X[,2])

# Initialize the inputs for QRc
Data = list(y=y, X=X, p=c(.05,.25,.5,.75,.95))
Mcmc = list(R=5000)

# Analyze 5 quantiles using default prior
out = QRseq(Method="QRc", Data=Data, Mcmc=Mcmc)

# Create default quantile plot of first variable
QRplot(out, var=1)

# Create quantile plot of second variable with 90% credible interval
QRplot(out, var=2, credint=c(.05, .95),main="This is an example")
}
