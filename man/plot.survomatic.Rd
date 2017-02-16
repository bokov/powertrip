\name{plot.survomatic}
\alias{plot.survomatic}
\title{Plot Survomatic Object}
\description{
Plot a \code{survomatic} object.
}
\usage{
plot.survomatic(xx, what = c("srv", "haz", "den"), col = c("black", "black", "darkred", "darkred"), lwd = c(2, 2, 2, 2), lty = c(1, 2, 1, 2), xlim = NULL)
}
\arguments{
  \item{xx}{
A \code{survomatic} object (as produced by the \code{\link{survwrapper}} function).
}
  \item{what}{
\code{"srv"} will plot a survival curve, \code{"haz"} will plot a hazard
function, and \code{"den"} will plot a density function.
}
  \item{col}{
Line colors for the first group data, the first group predicted values,
the second group data, and the second group predicted values, respectively.
}
  \item{lwd}{
Line widths for the first group data, the first group predicted values,
the second group data, and the second group predicted values, respectively.
}
  \item{lty}{
Line types for the first group data, the first group predicted values,
the second group data, and the second group predicted values, respectively.
}
  \item{xlim}{
Range for x-axis, passed to \code{\link{plot}}.
}
}
\details{
For \code{what="surv"} plots ordinary survival curves. For
\code{what="haz"} and \code{what="den"} also plots predicted values. To
suppress one or the other, set the first and third or the second and
fourth values of the \code{col} argument to \code{NA}.
}
\value{
\code{NULL}
}
%% \references{	    
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{survwrapper}}
}
\examples{
## Generate two sets of survival times.
population1 <- simsurv(250,type='g',p=c(7.33e-4,0.1227,0,0));
population2 <-
simsurv(200,type='lm',p=c(5.4818e-5,0.1543,0.0023,0.6018));
     
## Fit models to the populations and compare the parameters.
models1vs2 <- survwrapper(population1,population2);

## Plot the resulting survomatic object
plot(models1vs2)
}
\keyword{ survival }
\keyword{ hplot }

