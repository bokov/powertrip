\name{simsurv}
\alias{simsurv}
\title{
Simulate Survival Data
}
\description{
Generate a sample of \code{n} survival times from a distribution
specified by \code{type}, and parameters specified \code{p}.
}
\usage{
simsurv(n, type = "g", p = c(2.433083e-05, 0.005, 3e-11, 0.0015))
}
\arguments{
  \item{n}{
Number of observations to simulate.
}
  \item{type}{
One of \code{"e"} (Exponential), \code{"w"} (Weibull), \code{"g"} (Gompertz),
\code{"gm"} (Gompertz-Makeham), \code{"l"} (Logistic), \code{"lm"} (Logistic-Makeham)
}
  \item{p}{
A numeric vector of length 4 specifying the model parameters. Some
models ignore some parameters: \code{"e"} only uses the first parameter,
\code{"w"} and \code{"g"} use the first two, \code{"gm"} uses the first three,
\code{"l"} uses the first, second, \emph{and} fourth, and \code{"lm"} uses all
four. The model parameters are in the same order as they would appear in
\code{\link{survwrapper}} output.
}
}
\details{
If a parameter that is used by a model
is set to zero, the \code{type} argument is changed to that of a simpler
model that does not use that parameter.  \code{"l"} and \code{"lm"} use
the \pkg{Survomatic} functions \code{\link{rlgst}} and
\code{\link{rlogmake}}, respectively.  \code{"g"} and \code{"gm"} use the
\pkg{eha} functions \code{\link[eha]{rgompertz}} and
\code{\link[eha]{rmakeham}}, respectively.  \code{"e"} and
\code{"w"} use the \pkg{stats} functions \code{\link[stats]{rexp}} and
\code{\link[stats]{rweibull}} functions, respectively. Validation is not
enforced, but values outside [0,1] will probably give nonsensical
results or errors.
}
\value{
A numeric vector of length \code{n}.  
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

\seealso{
\code{\link{rlgst}}, \code{\link{rlogmake}},
\code{\link[eha]{rgompertz}}, \code{\link[eha]{rmakeham}},
\code{\link[stats]{rexp}}, 
\code{\link[stats]{rweibull}}
}
\examples{
## The following two are from the same Gompertz distribution because the
## zero parameter values cause the type to get changed from lm to g
s0 <- simsurv(1e4, type = "g", p=c(2e-5,5e-3,0,0))
s1 <- simsurv(1e4, type = "lm", p=c(2e-5,5e-3,0,0))
plot(survfit(Surv(s0)~1),col='red',conf.int=FALSE)
lines(survfit(Surv(s1)~1),col='blue',lty=2)
## ...but this one really is from a Logistic-Makeham distribution:
s2 <- simsurv(1e4, type = "lm", p=c(2e-5,5e-3,1e-4,0.4))
lines(survfit(Surv(s2)~1),col='green')
legend('topright',bty='n',col=c('red','blue','green'),lty=c(1,2,1),legend=c('s0','s1','s2'))
}
\keyword{ datagen }
\keyword{ distribution }
\keyword{ survival }
