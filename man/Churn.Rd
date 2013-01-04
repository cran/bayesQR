\name{Churn}
\docType{data}
\alias{Churn}
\title{Customer Churn Data}
\description{
This dataset is a stratified random sample from all active customers (at the end of June 2006) of a European financial services company. The dependent variable in this dataset is the churn behavior of the customers in the period from July 1st until December 31th 2006. Here a churned customer is defined as someone who closed all his/her bank accounts with the company. Note that all predictor variables are standardized. This dataset is a small subset of the dataset used in Benoit and Van den Poel (2013). The dataset is structured as a dataframe with 400 observations and 5 variables.
}
\usage{data("Churn")}
\format{
  The data frame has the following components:
  \itemize{
    \item \bold{churn} : churn (yes/no) 
    \item \bold{gender} : gender of the customer (male = 1) 
    \item \bold{Social_Class_Score} : social class of the customer 
    \item \bold{lor} : length of relationship with the customer 
    \item \bold{recency} : number of days since last purchase 
  }
}
\source{
Benoit, D.F. and Van den Poel, D. (2013). Quantile regression for database marketing: methods and applications. In: Coussement, K., De Bock, K.W. and Neslin, S.A. (eds.). \emph{Advanced database marketing: Innovative methodologies and applications for managing customer relationships}. Gower Publishing: London (UK), (in press).
}
\references{
Benoit, D.F. and Van den Poel, D. (2013). Quantile regression for database marketing: methods and applications. In: Coussement, K., De Bock, K.W. and Neslin, S.A. (eds.). \emph{Advanced database marketing: Innovative methodologies and applications for managing customer relationships}. Gower Publishing: London (UK), (in press).
}
\keyword{datasets}
