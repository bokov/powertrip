## Concise and hopefully more efficient replacement for modpars
## this one only does single groups
modparsng<-function(pars,min=c('e','g','gm','l','lm'),mout=c('e','g','gm','l','lm'),pad){
  ## pars = input params, AS USED BY optim; if you have them padded out to 4, strip them or modparsng them down first
  ## min, mout = the names of the in and out models
  ## pad = optional, what to fill the empty spaces with, if any
  min<-match.arg(min);mout<-match.arg(mout); pars<-c(pars);
  if(min==mout) return(pars);
  if(mout=='e') return(pars[1]);
  if(missing(pad)) nullvals<-c(NA,1e-4,1e-11,1e-9) else nullvals <-rep(pad,len=4);
  names(nullvals)<-c('a','b','c','s');
  if(min=='e') pars<-c(pars,nullvals[2]);
  if(mout=='lm') {
    if(min=='gm') return(c(pars,nullvals[4]));
    return(c(pars[1:2],switch(min,e=,g=nullvals[3:4],l=c(nullvals[3],pars[3]))));
  } else if(min!='lm') {
    return(c(pars[1:2],switch(mout,g=c(),gm=nullvals[3],l=nullvals[4])));
  } else switch(mout,g=pars[1:2],gm=pars[1:3],l=pars[c(1,2,4)]);
}
