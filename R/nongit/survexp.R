survexp <- function(pars,model='lm',subdivisions=1e6) integrate(srvshp,lower=0,upper=Inf,a=pars[1],b=pars[2],c=pars[3],s=pars[4],model=model,subdivisions=subdivisions)
