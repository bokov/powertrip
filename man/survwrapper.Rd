\name{survwrapper}
\alias{survwrapper}
\title{
Function for fitting and selecting hazard functions to survival data.
}
\description{
Takes one or two vectors of event times (numeric format) and optionally
corresponding vectors of indicator variables to designate right-censored
events. Fits several mortality models, selects the best fitting one/s,
and if two vectors were given, tests hypotheses about the model parameters.
}
\usage{
survwrapper(x, y = NULL, models = c("g", "gm", "l", "lm"), cx = rep(1, length(x)), cy = rep(1, length(y)), ext = F, n = length(c(x, y)), AIC = F, BIC = F, breakties = "AIC", compare.matrix = NULL, constraint.matrix = NULL, thresh = 0.05, smooth = 7)
}
\arguments{
  \item{x}{
A numeric vector of event times. For example, number of days an
individual has survived. 
}
  \item{y}{
An optional second numeric vector of event times, in the same units as
\code{x}. Specifying this vector enables fitting joint models and
hypothesis testing.
}
  \item{models}{
A character vector of model names: \code{g} = Gompertz, \code{gm} =
Gompertz-Makeham, \code{l} = logistic, \code{lm} = logistic-Makeham, and
\code{w} = Weibull. The default is to fit all the models except Weibull.
}
  \item{cx}{
A vector of 0 and 1 the same length as \code{x}. Ordinary events
(e.g. death) are indicated by 1 and censored events (e.g. withdrawal
from the study) are indicated by 0.
}
  \item{cy}{
A vector of 0 and 1 the same length as \code{y}. Ordinary events
(e.g. death) are indicated by 1 and censored events (e.g. withdrawal
from the study) are indicated by 0.
}
  \item{ext}{
Not implemented.
}
  \item{n}{
Total sample size. Should normally be left for the script to
automatically calculate, but can be specified when survwrapper is called
from another script repeatedly in order to speed up runtimes.
}
  \item{AIC}{
Whether to calculate the AIC (Akaike Information Criterion) for each
candidate model.
}
  \item{BIC}{
Whether to calculate the BIC (Bayes Information Criterion) for each
candidate model.
}
  \item{breakties}{
What criterion to use for choosing a model if more than one is justified
by the comparisons.
}
  \item{compare.matrix}{
A matrix for specifying a customized comparison algorithm.
}
  \item{constraint.matrix}{
A matrix of 1's and 0's for specifying a customized set of parameter
constraints to test.
}
  \item{thresh}{
Significance cutoff.
}
  \item{smooth}{
Not yet supported.
}
}
\details{
This function takes vectors (assumed to be times-to-event) and uses
numerical methods to find maximum likelihood estimates of model
parameters for a one or more models (by default, these are exponential,
Gompertz, Gompertz-Makeham, logistic, and Logistic-Makeham). Censored
events can be specified with the \code{cx} and \code{cy} arguments,
which should be vectors of 1's and 0's the same length as \code{x} and
\code{y}, respectively. If both \code{x} and \code{y} are specified, the
best joint model/s are chosen (such that the same type of model is
fitted to both populations, and the likelihood ratio between this joint model
and the corresponding joint model but one fewer parameter is
significantly different from 1 according to the chi-squared
distribution). Then, for each parameter in the model/s so chosen, a test
is performed on the null hypothesis that constraining the parameter
such that it is identical between the two populations will result in a
joint model that does not fit significantly worse than the full
model. If this null hypothesis is rejected, the interpretation is that
the corresponding parameter significantly differs between
populations. The null hypothesis of all parameters not being significantly
different is also tested by default. The user can also specify which
hypotheses to test.
}
\value{
  \item{x.m}{
    A data frame containing a column for the log-likelihood and each of
    the model parameter estimates. Each model fitted is represented by
    its own row.
  }
\item{y.m}{
  If \code{y} was specified, a data frame like \code{x.m} but for the
  \code{y} population.
}
\item{xy.sm}{
  If \code{y} was specified, a data frame showing the results of a
  series of full-reduced tests where the reduced (\code{simple}) joint
  model is compared to the full (\code{complex}). \code{df} indicates
  the different in degrees between the complex and simple
  models. \code{LLc} and \code{LLs} are the log-likelihoods for the
  complex and simple models, respectively. \code{LR} is the likelihood
  ratio between the complex and simple models. \code{chi2} is the chi
  squared statistic calculated from this ratio. \code{npar} is the
  number of parameters in the complex model. \code{AIC} and \code{BIC}
  are respectively the Akaike and Bayes information criteria for the
  complex model. \code{p} is the P-value for the full-reduced test and
  \code{sig} is whether this was deemed significant, and \code{chosen}
  is whether, based on the overall pattern of significant and
  non-significant full-reduced tests, the complex joint model for that
  row was the one selected for testing hypotheses about parameter
  differences.
}
\item{par.differences}{
  IN PROGRESS
}
\item{x}{
  The original values in the \code{x} argument.
}
\item{y}{
  The original values in the \code{y} argument, if any.
}
\item{cx}{
  The original values in the \code{cx} argument.
}
\item{cy}{
  The original values in the \code{cy} argument, if any.
}
\item{x.d}{
  IN PROGRESS
}
\item{y.d}{
  IN PROGRESS
}
\item{suggested.models}{
  A string with the abbreviation of the joint model/s chosen for hypothesis testing.
}
\item{nx}{
  An integer representing the sample size of the \code{x} group.
}
\item{ny}{
  An integer representing the sample size of the \code{y} group, if any.
}
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\author{
Alex F. Bokov (bokov@uthscsa.edu), Jon A. Gelfond
}
\note{
Uses Nelder-Mead algorithm to find maximum likelihood estimates of model
parameters.
}

\seealso{
  \code{\link{optim}}
}
\examples{
## Generate two sets of survival times.
population1 <- simsurv(629,type='g',p=c(7.33e-4,0.1227,0,0));
population2 <-
  simsurv(574,type='lm',p=c(5.4818e-5,0.1543,0.0023,0.6018));

## Fit models to the populations and compare the parameters.
models1vs2 <- survwrapper(population1,population2);
}
\keyword{ survival }
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
