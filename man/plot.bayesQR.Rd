\name{plot.bayesQR}
\alias{plot.bayesQR}
\title{ Produce quantile plots or traceplots with \code{plot.bayesQR} }
\description{
  \code{plot.bayesQR} is an S3 method that produces quantile plots or traceplots based on the estimates obtained by the \code{bayesQR} function. 
  For quantile plots, note that the more quantiles are estimated with \code{bayesQR}, the more detailed the quantile plot will be.
}
\usage{
  \method{plot}{bayesQR}(x, var=NULL, quantile=NULL, burnin=0, credint=c(0.025,0.975), plottype=NULL, main=NULL, xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, ...)
}
\arguments{
  \item{x}{ an output object of the \code{bayesQR} function, i.e. an S3 object of class bayesQR. }
  \item{var}{ vector containing the index or name of the variable(s) that has to be plotted (default=all variables). }
  \item{quantile}{ vector containing the quantile(s) that has to be plotted (default=all quantiles). }
  \item{burnin}{ the number of burnin draws that should be discared (default=0, meaning all draws are included). }
  \item{credint}{ the width of the credible interval (default=c(0.025, 0.975)). }
  \item{plottype}{ should be 'quantile' or 'trace'. }
  \item{main}{ Main title of the plot (default="").  }
  \item{xlab}{ Label of the x-axis; if omitted, the value chosen based on the input data. }
  \item{ylab}{ Label of the y-axis; if omitted, the value chosen based on the input data. }
  \item{xlim}{ Plot region of the x-axis; if omitted, the value chosen based on the input data. }
  \item{ylim}{ Plot region of the y-axis; if omitted, the value chosen based on the input data. }
  \item{...}{ additional arguments that are passed to the generic \code{plot} function }
}
\value{
  A (series of) quantile plot(s) or a (series of) traceplot(s). 
}
\author{
  Dries F. Benoit
}
\details{
  A quantile plot shows how the value of the regression parameter changes over a range of quantiles together with the associated credible interval.
	Note that the more quantiles are estimated, the more detailed the quantile plot will be.
  The minimum number of quantiles to plot is 2.
  A traceplot gives the evolution of the MCMC chain and can be used as graphical check of convergence.
  Note that more formal checks of convergence exist (see, for example, the \pkg{coda} package).
}
\examples{
# Simulate data from heteroskedastic regression
set.seed(66)
n <- 200
X <- runif(n=n,min=0,max=10)
X <- X
y <- 1 + 2*X + rnorm(n=n, mean=0, sd=.6*X)

# Analyze 5 quantiles using default prior
out <- bayesQR(y ~ X, quantile=c(.05,.25,.5,.75,.95), ndraw=5000)

# Check traceplot of first variable of .75 quantile regression 
plot(out, var=1, quantile=.75, plottype="trace")

# Create default quantile plot of first variable
plot(out, var=1, plottype="quantile")

# Create quantile plot of second variable with 90% credible interval
plot(out, var="X", credint=c(.05, .95), plottype="quantile", main="This is an example")
}
