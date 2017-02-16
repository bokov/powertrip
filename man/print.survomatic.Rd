\name{print.survomatic}
\alias{print.survomatic}
\alias{summary.survomatic}

\title{
Succinctly display the survomatic input most likely to be of interest.
}
\description{
At the moment, this is just a rough placeholder function that mainly
exists to prevent survomatic object from spamming the 
}
\usage{
print.survomatic(xx)
summary.survomatic(xx)
}

\arguments{
  \item{xx}{
A \code{survomatic} object.
}
}
\details{
Not a lot to be said yet; \code{print.survomatic} and
  \code{summary.survomatic} are at the moment the same S3 method and
  return a table with parameter estimates and standard errors copied
  from \code{x.m} and \code{y.m} and hypothesis test information
  from \code{par.differences}. \code{x.m}, \code{y.m}, and
  \code{par.differences} are objects inside a \code{survomatic} object.
}
\value{
    A list containing a data frame with 8 columns of numeric values and
  as many rows as there are parameters in the model plus one more
  representing a test of the null hypothesis that all parameters are the
  same between the two groups. The row names should clarify what
  hypothesis is being tested in each case.
  \describe{
    \item{\code{Xest}}{Parameter estimates for the first group being
  compared.}
    \item{\code{Xse}}{Standard errors for the first group being compared.}
    \item{\code{Yest}}{Parameter estimates for the second group being compared.}
    \item{\code{Yse}}{Standard errors for the first group being compared.}
    \item{\code{LR}}{Likelihood ratio for the hypothesis being tested on
  this row.}
    \item{\code{df}}{Degrees of freedom for the likelihood ratio test
  corresponding to this row.}
    \item{\code{chi2}}{Chi-squared statistic from the likelihood ratio
  test corresponding to this row.}
    \item{\code{p}}{p-value from the likelihood ratio test corresponding
  to this row.}
  }
}
% \references{
% %% ~put references to the literature/web site here ~
% }
% \author{
% %%  ~~who you are~~
% }
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{survwrapper}}
}
\keyword{ htest }
\keyword{ print }
\keyword{ survival }
