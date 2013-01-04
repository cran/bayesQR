\name{Prostate}
\docType{data}
\alias{Prostate}
\title{Prostate Cancer Data}
\description{
  These data come from a study that examined the correlation between the
  level of prostate specific antigen and a number of clinical measures
  in men who were about to receive a radical prostatectomy.  It is a data
  frame with 97 rows and 9 columns.
}
\usage{data("Prostate")}
\format{
  The data frame has the following components:
  \itemize{
    \item \bold{lcavol} : log(cancer volume)
    \item \bold{lweight} : log(prostate weight)
    \item \bold{age} : age
    \item \bold{lbph} : log(benign prostatic hyperplasia amount)
    \item \bold{svi} : seminal vesicle invasion
    \item \bold{lcp} : log(capsular penetration)
    \item \bold{gleason} : Gleason score
    \item \bold{pgg45} : percentage Gleason scores 4 or 5
    \item \bold{lpsa} : log(prostate specific antigen)
  }
}
\source{
  Stamey, T.A., Kabalin, J.N., McNeal, J.E., Johnstone, I.M., Freiha,
  F., Redwine, E.A. and Yang, N. (1989).
  Prostate specific antigen in the diagnosis and treatment of
  adenocarcinoma of the prostate: II. radical prostatectomy treated
  patients,
  \emph{Journal of Urology} \bold{141}(5), 1076--1083.
}
\references{

  Stamey, T.A., Kabalin, J.N., McNeal, J.E., Johnstone, I.M., Freiha,
  F., Redwine, E.A. and Yang, N. (1989).
  Prostate specific antigen in the diagnosis and treatment of adenocarcinoma 
  of the prostate: II. radical prostatectomy treated patients,
  \emph{Journal of Urology}, \bold{141}(5), 1076--1083.
\cr
\cr
  This dataset also appears in the R-package: \emph{lasso2}. \cr
  Lokhorst, J., Venables, B. and Turlach, B. (2011).
  lasso2: L1 contstrained estimation aka 'lasso',
  published on R CRAN website: 
  \emph{http://cran.r-project.org/web/packages/lasso2/index.html}.
}
\keyword{datasets}
