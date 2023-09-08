\name{summary.bayesQR}
\alias{summary.bayesQR}
\title{ Summarize the output of the \code{bayesQR} function }
\description{
  \code{summary.bayesQR} is an S3 method that summarizes the output of the \code{bayesQR} function in an object of class bayesQR.summary.
  For every estimated beta and sigma, this object contains the Bayes estimate and the posterior credible interval is calculated. 
  The object also contains other relevant information about the estimation procedure, such as the quantile, the variable names, etc.
}
\usage{
  \method{summary}{bayesQR}(object, burnin=0, credint=c(.025,.975), quantile=NULL, ...)
}
\arguments{
  \item{object}{ an output object of the \code{bayesQR} function, i.e. an S3 object of class bayesQR. }
  \item{burnin}{ the number of burnin draws that should be discarded (default=0, meaning all draws are included). }
  \item{credint}{ the width of the credible interval (default=c(0.025, 0.975)). }
  \item{quantile}{ the quantile(s) of the quantile regressions that have to be summarized (default: all estimated quantiles in QRobj). }
  \item{...}{ additional arguments passed to the generic \code{summary} function. }
}
\value{
  An object of class \code{bayesQR.summary}, basically a list including elements:
  \item{method}{ a string containing the method that was used, i.e. indicating whether the dependent variable was continuous or binary and whether adaptive lasso variable selection was used. }
  \item{normal.approx}{ logical flag for normal approximation of posterior distribution. }
  \item{quantile}{ the quantile that was estimated. }
  \item{names}{ character vector containing the names of the independent variables in the model. }
  \item{burnin}{ the number of burnin draws that were discarded. }
  \item{retained}{ the number of draws that were retained and used to calculate the summary statistics. }
  \item{credint}{ the width of the credible interval. }
  \item{betadraw}{ the Bayes estimate, credible interval and, if \code{normal.approx=TRUE}, the adjusted credible intervals of the beta draws. }
  \item{sigmadraw}{ the Bayes estimate and credible interval of the sigma draws. }
}
\author{
  Dries F. Benoit
}
\examples{
# Load the Prostate cancer dataset
data(Churn)

# Estimate the model parameters with default prior
out <- bayesQR(churn~gender+recency, data=Churn, ndraw=2000)

# Return Bayes estimates and credible intervals 
sum <- summary(out, burnin=1000)

# Inspect structure of bayesQR.summary object
str(sum)

# Print bayesQR.summary object
sum
}
