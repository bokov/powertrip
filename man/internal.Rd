\name{survomatic-internal}
\alias{weibhaz}
\alias{rlgst}
\alias{rlogmake}
\alias{srvhaz}
\alias{srvshp}
\alias{cpars}
\alias{smoother}
\alias{modelinfo}
\alias{modelpars}
\alias{gmean}
\alias{checkcons}
\alias{fixpar}
\alias{objf}

\title{
Functions that are currently internal to Survomatic.
}

\description{
  Various functions that are internally used by
  \code{\link{survwrapper}} and other front-end functions.
}
\usage{

weibhaz(x, a, b, c = 0, s = 0, i = 1) 
rlgst(n, a, b, s) 
rlogmake(n, a, b, c, s) 
srvhaz(x, a, b, c = 0, s = 0, i = 1) 
srvshp(x, a, b, c = 0, s = 0, i = 1, model = "g") 
cpars(x, y, xpars, ypars, model, LL, const = NULL, cx = rep(1, 
    length(x)), cy = rep(1, length(y))) 
smoother(x, s = smooth) 
modelinfo(x, what = c("comp", "rcomp", "dfcomp", "dfrcomp", "tblcomp", 
    "tblrcomp", "tblcomp2", "deps", "rdeps", "tbldeps", "tblrdeps", 
    "vars", "keep", "fullname"), extended = T) 
modpars(x, modeli, modelo = NULL, cni = rep(T, 4), cno = NULL, 
    nil = 1e-07, wbnil = 1, pf = mean, trim = F, onegrp = F) 
gmean(x) 
checkcons(cons, np, keep) 
fixpar(par, np, keep) 
objf(par, lb = 0, ub = 10, cons = rep.int(1, 4), x, y = NULL, 
    keep = 1:2, np = 2, ex = "ofg", nx, ny, cx, cy, tlog = F) 
}

\arguments{
  \item{a}{The a parameter, in literature sometimes called lambda.}
  \item{b}{The b parameter, in literature sometimes called gamma.}
  \item{c}{The c parameter (exogenous mortality).}
  \item{s}{The s parameter (late mortality deceleration).}
  \item{i}{Unit conversion factor by which the raw observations are multiplied} 
  \item{n}{How many values to return}
  \item{model}{Model to use, can be \code{'e'}, \code{'w'}, \code{'g'}, \code{'gm'},
  \code{'l'}, or \code{'lm'}} 
  \item{cni}{Vector of 1s and 0s. The 0s are parameters were constrained to be equal between
  the groups in the input model.}
  \item{cno}{Optional vector of 1s and 0s. The 0s are the parameters
  that should be constrained to be equal between the groups in the
  output model}
  \item{cons}{Parameters that should be held constant}
  \item{const}{Parameters that should be held constant}
  \item{cx}{Censored observations from the first group}
  \item{cy}{Censored observations from the second group}
  \item{ex}{Name of the objective function to try}
  \item{extended}{Whether to return extended model info}
  \item{keep}{Which parameters to keep}
  \item{lb}{Lower bound}
  \item{ub}{Upper bound}
  \item{modeli}{Model converting from}
  \item{modelo}{Model converting to}
  \item{nil}{A number with which to fill formerly variable positions
  that were \code{NA} in the input model but must exist in the output model.}
  \item{wbnil}{A weibull-specific argument that, if it exists and the
  model is weibull, will overwrite \code{nil}.}
  \item{np}{Number of parameters.}
  \item{nx}{Number of observations in first group}
  \item{ny}{Number of observations in second group}
  \item{onegrp}{If \code{TRUE}, then parameters from just one
  group. Otherwise (default) parameters from a two-group joint model.}
  \item{par}{A vector of model parameters}
  \item{pf}{Parameter function. If given a pair of parameters an unconstrained model, how to
  find a starting value for those parameters in a model where they are
  constrained to equal each other.}
  \item{tlog}{If true in \code{objf}, the parameters are exponentiated first.}
  \item{trim}{Whether to strip \code{NA}s from the \code{modpars} output.}
  \item{what}{What type of information is needed from \code{modelinfo}.}
  \item{x}{Event times from first group.}
  \item{y}{Event times from second group.}
  \item{xpars}{To be added.}
  \item{ypars}{To be added.}  
  \item{LL}{To be added.}
}

\details{
Internal functions.
}

\keyword{ internal }
