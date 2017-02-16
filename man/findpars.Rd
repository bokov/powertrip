\name{findpars}
\Rdversion{1.1}
\alias{findpars}
\title{
Fit and compare mortality models.
}
\description{
Fit and compare all supported models to a dataset evaluate the likelihood of a parameter differing between two datasets. Instead of using 'reasonable' starting values, this function takes advantage of the fact that many mortality models are nested in each other and all are ultimately nested in the exponential model whose sole parameter can be calculated directly. Each successively more complex model uses the parameters fitted to the model it is nested in as starting values.
}
\usage{
findpars(x, y = NULL, cx = NULL, cy = NULL, nil = 0, bnil = 0, wbnil = 1, pf = mean, label = NULL, summary = F, id = 0, tlog = F, digits = 22, sig = 0.05, models = NULL)
}
\arguments{
  \item{x}{
A vector of ages at death (integers, usually representing days).
}
  \item{y}{
An optional second vector of ages at death (integers, usually representing days).
}
  \item{cx}{
A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the control group.
}
  \item{cy}{
A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the experimental group.
}
  \item{nil}{
A starting numeric value for a newly added parameter. For example, the Logistic model uses the best-fit parameters for the Gompertz model as its starting values, but the Gompertz model has two parameters while the Logistic model has three. The third parameter takes on the value of the nil argument.
}
  \item{bnil}{
Like nil, but bnil only gets substituted into the b parameter.
}
  \item{wbnil}{
The numeric value to substitute for missing values of the 'b' parameter when finding parameters for the Weibull model. Defaults to 1.
}
  \item{pf}{
A function to call when two different parameters are constrained and need to produce a single starting value. In addition to \code{mean} (default), other valid functions include \code{median}, \code{gmean}, \code{max}, and \code{min}.
}
  \item{label}{
A character prefix based on which the output files are named. You may specify a file-path in the prefix, but do not specify a file extension as the extension '.txt' will automatically be appended to the filenames. If the label is not specified, findpars will not save an output file.
}
  \item{summary}{
Logical value telling findpars whether or not to echo a summary table to the console.
}
  \item{id}{
A numeric value, copied to the 'id' column of the output table. This is used to distinguish different iterations when findpars is called by another function or script.
}
  \item{tlog}{
Logical value passed as the tlog argument to the opsurv function, instructing it to find best fits for the natural log of the starting values and then convert them back before outputting them. Maybe result in improved speed but decreased accuracy and is still being tested.
}
  \item{digits}{
Integer specifying many digits to show in the summary output to the console. Only useful if summary is set to TRUE.
}
  \item{sig}{
Significance cutoff used by the chi-squared test.
}
  \item{models}{
Optional character vector containing the names of models to test. If you specify only the models you're interested in, findpars is smart enough to automatically fill in the prerequisite models and fit them as well.
}
}
\details{
This function takes one or two vectors of ages at death, fits various mortality models to them, and compares the fits of each mortality model to the next simpler mortality model that it's based on to find the best models. If two vectors (x and y) have been given, for each parameter of each joint model, this function will compare the unconstrained joint model to a joint model that is constrained such that that parameter is the same between the two groups. A significantly worse fit than the unconstrained model indicates the parameter being tested is in fact different. While running the function draws a progress bar on the console and prints the name of a model (the unconstrained version of it) each time a fit is complete. This is normal behavior whose purpose is to reassure the user that the function hasn't hung and is indeed working, because for some datasets it may take a long time to converge on the best parameters.
}
\value{
A data.frame, invisibly returned. If summary=T, the table is also printed to the console. If a label argument is given, the data is also saved to a tab-delimited file of that name. Cells in the table for which there are no applicable values are filled with NAs. The output contains the following columns:
\item{MLE}{Maximum likelihood estimate for the model described in that row.}
\item{# cycles}{The number of iterations the optimization function had to go through before it converged.}
\item{OK?}{If this model has a better fit than the next simplest model, this value will be '1'. If the fit is the same, this value will be 0. If the fit is worse, this value will be -1 and in such cases you are encouraged to send your data and your command history to the author because you have most likely uncovered a bug in the software.}
\item{a1}{The lambda parameter of the first group, for all models.}
\item{b1}{The beta parameter of the first group for a Weibull model and the gamma parameter of the first group for a Gompertz, Logistic, Gompertz-Makeham, or Logistic-Makeham model.}
\item{c1}{The c parameter of the first group for a Gompertz-Makeham or Logistic-Makeham model.}
\item{s1}{The s parameter of the first group for a Logistic or Logistic-Makeham model.}
\item{a1}{The lambda parameter of the second group, for all models.}
\item{b1}{The beta parameter of the second group for a Weibull model and the gamma parameter of the second group for a Gompertz, Logistic, Gompertz-Makeham, or Logistic-Makeham model.}
\item{c1}{The c parameter of the second group for a Gompertz-Makeham or Logistic-Makeham model.}
\item{s1}{The s parameter of the second group for a Logistic or Logistic-Makeham model.}
\item{model}{The model being fitted in this row.}
\item{id}{A copy of the 'id' argument.}
\item{nil}{A copy of the 'nil' agrument.}
\item{bnil}{A copy of the 'bnil' argument.}
\item{LR}{Log ratio for this fit-- i.e. log(MLE this model) - log(MLE next simpler model).}
\item{AIC}{Akaike's Information Criterion for this model fit. This can be used for comparing the fit of models that are not nested in each other.}
\item{BIC}{Bayesian Information Criterion for this model fit, which is similar to the AIC but with the additional feature of penalizing for the number of parameters.}
\item{p (chi squared)}{The significance level of the observed LR assuming a chi-squared null distribution of LRs with one degree of freedom.}
\item{sig?}{Whether or not the significance level is below the threshold set by the 'sig' argument.}
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\author{
Alex F. Bokov
}

\seealso{
\code{\link{opsurv}}, \code{\link{simdist}}, \code{\link{empdist}}
}

\keyword{ survival }