\name{simdist}
\Rdversion{1.1}
\alias{simdist}
\title{
Function that uses models fitted to the original data to obtain simulated null distributions that can in turn be used to select the most appropriate model for a population.
}
\description{
This function identifies the best-fit parameters for a sample, then uses these parameters to obtain null distributions of MLE log-ratios by sampling from the population that this sample would represent if each candidate model were, in turn, correct.
}
\usage{
simdist(x, label, cx = NULL, nil = 0, bnil = 0, wbnil = 1, pf = mean, rounds = 5000, models = c("w", "g", "gm", "l", "lm"), dropcols = c("a2", "b2", "c2", "s2"), sims = NULL, pars = NULL)
}
\arguments{
  \item{x}{
An integer vector of ages at death.
}
  \item{label}{
A character string representing the base name for the output files. Each output file has the same format as the output from \code{\link{findpars}}, stacked together by rows, representing every simulation. There is a separate such output file for each model selected (in the \code{models} argument) and its name is 'LABEL_M' where 'LABEL' is the value of this argument and 'M' is the abbreviation of the model name ('w','g','gm','l', or 'lm').
}
  \item{cx}{
A vector of ones and zeros, with ones representing a natural death event, and zeros representing a censored event. If omitted, it is assumed that all the deaths are natural in the control group.
}
  \item{nil}{
The numeric value to substitute missing values of fitted parameters in order to avoid rounding errors caused by machine precision limits. That problem has since been solved and this argument may be removed in future versions.
}
  \item{bnil}{
The numeric value to substitute for missing values of the 'b' parameter in order to avoid rounding errors caused by machine precision limits. That problem has since been solved and this argument may be removed in future versions.
}
  \item{wbnil}{
The numeric value to substitute for missing values of the 'b' parameter when finding parameters for the Weibull model. Defaults to 1.
}
  \item{pf}{
A function to call when two different parameters are constrained and need to produce a single starting value. In addition to \code{mean} (default), other valid functions include \code{median}, \code{gmean}, \code{max}, and \code{min}.
}
  \item{rounds}{
How many samples to simulate.
}
  \item{models}{
A character vector of models to try which can be any combination of: 'w','g','gm','l', or 'lm'. In some cases, models that are not listed in this argument are fitted anyway because they are needed to obtain starting parameters for the models that are listed.
}
  \item{dropcols}{
Columns to drop from the output. By default this is \code{c('a2','b2','c2','s2')} because \code{simdist} is designed for just one sample rather than a comparison of two samples, and therefore those columns will never be used by this function. The name of any unwanted output columns can be included in this argument however.
}
  \item{sims}{
Simulating samples is a time consuming process. If the simulations from a previous run exist as an R object in the current environment, that object can be specified in the \code{sims} argument, causing the simulation step to be bypassed and all the simulated data to instead be taken from the object.
}
  \item{pars}{
Similarly to the \code{sims} argument, \code{pars} allows you to specify an already calculated collection of parameters. This feature is experimental.
}
}
\details{
\code{simdist} calls \code{findpars} on a sample in order to obtain estimates of each parameter for each model of interest as well as log-ratios of MLEs for the respective model comparisons. This function then simulates samples of the same size based on the fitted parameters, and does this \code{rounds} times. Then, findpars is run again on each model of interest, to obtain a distribution of MLE log-ratios. The original MLE log-ratio is compared against this distribution, in order to obtain a more robust estimate of significance level.
}
\value{
This function does not return anything to the console, but instead saves its output as files. The following files are generated (replacing LABEL with the value of the \code{label} argument):
\item{LABEL.rdata}{An R data file containing the output from the initial (non-simulated) \code{findpars} in a data.frame called \code{ihaz}.}
\item{LABEL.sims.rdata}{An R data file with all the simulated samples in a list object called \code{sims} containing one N x rounds matrix of survival times for each model used as the null hypothesis. It also contains a list object called \code{pars} that contains the parameter estimates from all the fitted models labeled AB where A is the null model and B is the alternative model. This file can be opened from R using the \code{load()} command but is not usually necessary unless one is trying to recreate the results of a previous simulation.}
\item{LABEL_M}{A collection of tab delimited text files is generated, one for each model tested against its respective null hypothesis models. The 'M' in the label name is replaced with the abbreviation of the model name (by default, 'g', 'gm', 'l', and 'lm'). These files are like those produced by \code{findpars()}} except they lack any of the columns named in the \code{dropcols} argument and they have the following additional columns:
\item{null_pars.a1}{The best fit for the \code{a} parameter of the null model.}
\item{null_pars.b1}{The best fit for the \code{b} parameter (if applicable) of the null model.}
\item{null_pars.c1}{The best fit for the \code{c} parameter (if applicable) of the null model.}
\item{null_pars.s1}{The best fit for the \code{s} parameter (if applicable) of the null model.}
\item{null_model}{An abbreviation identifying the model representing the null hypothesis (i.e. a model that has one less parameter than the target model, which is the one that is being compared to it).}
\item{target_model}{An abbreviation identifying the model representing the alternative hypothesis (i.e. a model that has one more parameter than the null model, which it is being compared to).

The values in the 'id' column correspond to which column (from the left) in the corresponding 'sims' object for that model comparison contains the data that produced that given row.
}
}
\references{
Pletcher,S.D., Khazaeli,A.A., and Curtsinger,J.W. (2000). Why do life spans differ? Partitioning mean longevity differences in terms of age-specific mortality parameters. Journals of Gerontology Series A-Biological Sciences and Medical Sciences 55, B381-B389
}
\note{
Even though the examples use 'rounds=100' for demonstration purposes, for publication it is recommended to leave the 'rounds' argument at its default of 5000 or give it a larger value. This may take several hours or even days depending on the dataset. This is why the output is saved to a text file rather than an object within R-- the user can monitor the progress of this function by accessing the data file with \code{\link{constrshow}} from another R session.

Currently this function spams the console with periods and abbreviated model names. This behavior is intentional, for troubleshooting purposes and to indicate that R has not crashed or hung.
}
\author{
Alex F. Bokov
}
\seealso{
\code{\link{findpars}}, \code{\link{modelshow}}
}
\section{Warning}{
These examples may take 10 minutes or more to finish running.
}
\keyword{ survival }
