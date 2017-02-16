\name{opsurv}
\Rdversion{1.1}
\alias{opsurv}
\title{
Optimization function for finding the set of parameters that return the largest MLE (maximum likelihood estimate) given a mortality model (that can be constrained) and a set of observed data.
}
\description{
This function takes one or two datasets, starting parameters, constraints, and various other user settings and formats them in a manner that can be repeatedly sent to the \code{\link{optim}} function until the MLE stops improving.
}
\usage{
opsurv(x, y = NULL, model = "g", par = c(2.6e-06, 0.004, 1e-07, 0.1), cons = 1, usegr = F, usehs = F, debug = F, lb = c(1e-14, 1e-04, 0, 0), ub = c(0.5, 0.5, 0.5, 2), cx = NULL, cy = NULL, mvers = "", method = "Nelder-Mead", tlog = F)
}
\arguments{
  \item{x}{
A vector of ages at death (integers, usually representing days).
}
  \item{y}{
An optional second vector of ages at death (integers, usually representing days). Null by default.
}
  \item{model}{
A character string that can be one of 'w','g','gm','l', or 'lm' (Weibull, Gompertz, Gompertz-Makeham, Logistic, and Logistic-Makeham, respectively) and specifies the model to fit parameters for.
}
  \item{par}{
Starting parameters, which are a vector of four values for the a, b, c, and s parameters respectively. If not left at its default value, all four parameters must be specified but the values of anything other than the first parameter can be 0.
}
  \item{cons}{
Constraints-- i.e., which of the four parameters should be held equal between the two groups. A vector of four logical or numeric values (which will be coerced to logical) but if a single value is given, it gets copied to a vector of length four. A value of FALSE means that parameter will be fitted such that it is the same between the 'x' and 'y' group. More than one parameter can be constrained at the same time. The parameters are specified in the following order: a, b, c, s, therefore cons=c(1,0,1,1) would be interpreted as the b parameter being the same between the two populations. Parameters that are not used by a given model (for example, c or s in the case of the Gompertz model) are ignored. It is meaningless to specify constraints if only the 'x' group is specified in the arguments.
}
  \item{usegr}{
Logical value specifying whether or not to calculate a gradient.
}
  \item{usehs}{
Logical value specifying whether or not to calculate a hessian.
}
  \item{debug}{
Logical value specifying whether or not to turn on debug options. Not supported, and subject to change without notice.
}
  \item{lb}{
A vector of four numeric values, specifying the lower feasibility bound for the fitted parameters.
}
  \item{ub}{
A vector of four numeric values, specifying the upper feasibility bound for the fitted parameters.
}
  \item{cx}{
A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the control group.
}
  \item{cy}{
A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the experimental group.
}
  \item{mvers}{
String representing the version of a function to call. Used for debugging and not currently applicable.
}
  \item{method}{
String representing the optimization method to use in finding the best-fit parameters. Currently only Nelder-Mead is supported and is the default.
}
  \item{tlog}{
Logical value indicating whether or not the starting parameters should be converted to log form before the optimization step. Defaults to FALSE.
}
}
\details{
In principle, the only argument that needs to be specified is 'x', and if you wish to fit a model other than Gompertz, then also the 'model' argument. If both 'x' and 'y' are specified, then constraints can also be specified. It is easy to make this function error out by giving it bad starting values. Furthermore, it is tedious to call \code{\link{opsurv}} for each possible model and set of constraints. For these reasons it is recommended that \code{opsurv} be invoked via \code{\link{findpars}} rather than used directly.
}
\value{
A list, with the following items:
\item{estimate}{The set of parameters that gives the best MLE. Only the parameters used by the model you chose are returned. If a parameter is constrained, it is omitted from the first group, so where an unconstrained model would have N parameters shown here, a model with one constraint will have N-1. }
\item{maximum}{Maximum likelihood estimate (MLE).}
\item{iterations}{The number of times the Nelder-Mead optimimization algorithm had to run on the final iteration before it converged.}
\item{code}{Integer code, where '0' indicates successful convergence and anything else is an error state. For more information, see the documentation for \code{\link{optim}}.}
\item{message}{Message returned by \code{\link{optim}}. In the absence of errors should be NULL.}
\item{gradient}{If usegr is TRUE, the gradient for the selected model, under the selected constraints, using the fitted parameters. }
\item{titer}{The total number of iterations for all the times the Nelder-Mead optimization had to run (Nelder-Mead is prone to stopping prematurely, so it is called repeatedly, using the previous estimate as the starting value for the new run, until two successive estimates are the same to the limit of machine precision).}
\item{runtime}{Total amount of time the process took to arrive at the estimate.}
\item{group0}{If both 'x' and 'y' are specified but there are no constraints, they are run singly as unconstrained models. This item is a copy of the usual output list for 'x' only.}
\item{group1}{If both 'x' and 'y' are specified but there are no constraints, they are run singly as unconstrained models. This item is a copy of the usual output list for 'y' only.}
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\author{
Alex F. Bokov
}
\seealso{
\code{\link{findpars}}, and \code{\link{objf}}, which is the function that actually evaluates the MLE at each step of the optimization process.
}
\keyword{ survival }
