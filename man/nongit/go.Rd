\name{go}
\alias{go}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Automatic analysis of survival/longevity data. }
\description{
  The main function of Survomatic, which takes the raw survival data and performs log-rank, quantile regression and mortality modeling analysis.
}
\usage{
go(x, y, xynames = c(), save = T, path, prompt = 2, xlim = 1350, slwd = 4, scol = c("darkred", 1), spch = 24:25, spcex = 0, splwd = slwd, spbg = scol, sleg = T, sxcex = 1.5, sxlwd = 3, qstart = 1, demint = 30, qrse = "boot", dohz = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Vector of ages at death (in days) for control group.}
  \item{y}{Vector of ages at death (in days) for experimental group.}
  \item{xynames}{Character vector with labels for the two groups.}
  \item{save}{Boolean value of whether to save the various output files or not.}
  \item{path}{Path for where to save the output files.}
  \item{prompt}{How aggressively to prompt you to make choices, versus deciding on your behalf.}
  \item{xlim}{The length of the x axis for the Kaplan-Meier curves.}
  \item{slwd}{Line width for the survival curve plots.}
  \item{scol}{Line colors for the survival surves, first the experimental group and second the control.}
  \item{spch}{Shapes to represent data points on the survival curves.}
  \item{spcex}{Size of the data points on the survival curves (set to 0 by default and therefore the survival curves are by default plotted as lines only).}
  \item{splwd}{Thickness of the lines used to plot the data points on the survival curves.}
  \item{spbg}{Fill color for the data points on the survival curves.}
  \item{sleg}{Whether or not to plot a legend on the survival curves.}
  \item{sxcex}{Font size for the axes in the survival curves.}
  \item{sxlwd}{Line thickness of the axes in the survival curves.}
  \item{qstart}{Which collection of quantiles to attempt to use first in quantile regression-- an integer between 1 and 4 inclusive, and should not ordinarily be changed. This argument will probably removed in future version.}
  \item{demint}{By how many days to smooth the observed mortality rate. Should not be changed, and will probably be removed in future versions.}
  \item{qrse}{Algorithm for obtaining confidence intervals for quantile regression. See the help file for \code{summary.rq} for more information.}
  \item{dohz}{Boolean value controling whether or not to fit mortality models to the data.}
}
\details{
  This is intended to be a wrapper that automagically calls various other 
  functions on your survival data. At the moment these are log-rank, quantile regression, quantile score test, and fitting mortality models to the data.
}
\value{
  A list object with all potentially useful output about the comparison of the two populations. This list also contains a copy of go() which can be invoked as listname$go(), in which case it will print out and plot the previous output without recalculating it.
  \item{smry}{Quick summary of the quantiles of the two groups being compared.}
  \item{zsc}{A table with all the output of the quantile score test.}
  \item{x}{The ages at death of the control group.}
  \item{y}{The ages at death of the experimental group.}
  \item{xy}{The ages at death of the two groups together (first x, then y).}
  \item{xynames}{The names of the two groups.}
  \item{path}{The path where to save files by default.}
  \item{group}{A vector of group assignments the same length as xy}
  \item{lr}{Output from the log-rank test on the two groups.}
  \item{tt}{Output from the Student's t-test on the logs of the two groups}
  \item{tests}{The names of all the tests performed and the names of the items that contain the results.}
  \item{sig.tests}{The subset of tests which returned significant results.}
  \item{qreg}{Output from quantile regression.}
  \item{qreg.tab}{Table with quantile regression results with standard errors, test statistics, and p-values added.}
  \item{qreg.sum}{Output from the summary.rq functon.}
  \item{demint}{}
  \item{xd}{Demography table for x, including the observed hazard (ux).}
  \item{yd}{Demography table for x, including the observed hazard (ux).}
  \item{mod}{Output from mortality model fitting.}
  \item{report}{A table listing some summary statistics and the quantiles flagged by quantile regression and the score test as potentially significant.}
  \item{sigqreg}{Table with only the significantly different quantiles, values at that those quantiles, and p-values from qreg.tab.}
  \item{sigzsc}{Table with only the fraction surviving from each group at quantiles where a significant difference was observed and the p-value of that difference.}
  \item{xysurvfit}{Output from survfit of the two groups.}
  \item{slwd}{Line width for the survival curve plots.}
  \item{scol}{Line colors for the survival surves, first the experimental group and second the control.}
  \item{spch}{Shapes to represent data points on the survival curves.}
  \item{spcex}{Size of the data points on the survival curves (set to 0 by default and therefore the survival curves are by default plotted as lines only).}
  \item{splwd}{Thickness of the lines used to plot the data points on the survival curves.}
  \item{spbg}{Fill color for the data points on the survival curves.}
  \item{sleg}{Whether or not to plot a legend on the survival curves.}
  \item{sxcex}{Font size for the axes in the survival curves.}
  \item{sxlwd}{Line thickness of the axes in the survival curves.}
  \item{dohz}{Boolean value controling whether or not to fit mortality models to the data.}
  \item{go}{A copy of the function that originally generated the output, so the output can be used to plot/print itself.}
  \item{sys}{The exact phrasing of how this function was invoked in order to produce the current output.}
}


\references{ 
Pletcher, S. D. et al., (2000) Why Do Lifespans Differ? Partitioning Mean Longevity Differences in Terms of Age-Specific Mortality Parameters, _Journal of Gerontology: Biological Sciences_, *55A*, B381-B389

Wang, C. et al., (2004) Statistical Methods for Testing Effects on "Maximum Lifespan",_Mechanisms of Ageing and Development_, *125*, 629-632
}
\author{ Alex F. Bokov }
\note{Some older version of MacOS X do not load tcltk properly, and therefore will have problems running this package. Please make sure you are using the latest revision of Survomatic.}
\seealso{\code{\link{rq}},\code{\link{summary.rq}},\code{\link{surv2.logrank}}}
\examples{
# Can be directly called with no arguments, and it will prompt you for everything.
# Can also be called with arguments.
# a<-simsurv(20);
# b<-simsurv(20,p=c(1e-6,0.003,0,0));
# You can capture the output within R.
# output<-go(a,b,xynames=c('a','b'));
# output$go();
}
\keyword{ survival }