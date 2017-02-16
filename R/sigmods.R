sigmods <- function(x,y=NULL,n=NULL,comp='tblcomp2',breakties='AIC',AIC=T,BIC=T,cmat=NULL,thresh=.05){
  # x is a data.frame with the columns 'model' and 'LL'
  # the optional y is in the same format if it exists
  # we assume that x and y have the same number of rows
  # and in fact, identical 'model' columns
  # departures from this assumptions might still work, but
  # are not supported
  if(!is.null(y)) {
    x$LL <- x$LL+y$LL[match(x$model,y$model)];
    mult <- 2;
  } else mult <- 1;
  if(is.null(cmat)){
    cmat <- matrix(c(
                     1,0,NA,NA,NA,NA,NA,     # w ???
                     NA,1,0,0,0,NA,NA,       # g
                     NA,1,1,0,0,0,NA,        # gm (better than g)
                     NA,1,1,0,NA,0,1,        # gm (gm and lm better than g and l)
                     NA,1,0,1,0,NA,0,        # l (better than g)
                     NA,1,0,1,NA,1,0,        # l (lm and l better than gm and g)
                     NA,1,NA,NA,1,NA,NA,     # lm (better than g)
                     NA,1,1,NA,NA,1,NA,      # lm (better than gm)
                     NA,1,NA,1,NA,NA,1,      # lm (better than l)
                     NA,0,NA,NA,NA,NA,NA,    # e ???
                     NA,1,1,1,0,0,0),nrow=7, # true tie
                   dimnames=list(NULL,c('w','g','gm','gm','l','l','lm','lm','lm','e','tie')));
  }
  rownames(cmat) <- c('e.e','e.g','g.gm','g.l','g.lm','gm.lm','l.lm');
  m <- modelinfo(x$model,comp);
  # The 'tblcomp2' default argument to comp is magical, and m is a data.frame with simple, complex, and df
  cmat <- cmat[paste(m$simple,m$complex,sep='.'),c(grep(paste("^",x$model,"$",sep="",collapse="|"),colnames(cmat),val=T),'tie')];
  # The above is supposed to fix an error when fewer than a list of all the gompertz family models is passed
  # but now it causes incorrect model choices; for now can rely on AIC to resolve the ties
  df <- mult*m$df;
  LLc <- x$LL[match(m$complex,x$model)]; LLs <- x$LL[match(m$simple,x$model)];
  LR <- LLc - LLs;
  chi2 <- 2*LR; p <- pchisq(chi2,df,lower.tail=F);
  sig <- p<thresh;
  chosen <- unique(colnames(cmat)[apply(cmat==sig,2,all,na.rm=T)]);
  npar <- mult*sapply(m$complex,function(x) length(modelinfo(x,'var')));
  if(AIC|breakties=='AIC') AIC <- 2*npar-2*LLc;
  if(BIC|breakties=='BIC') BIC <- 2*log(n)*npar-2*LLc;
  if(!is.null(breakties)&(length(chosen)>1|chosen=='tie')){
    if(length(chosen)>1){
      chosen <- chosen[which.min(environment()[[breakties]][chosen])];
    } else if(chosen=='tie') chosen <- names(which.min(environment()[[breakties]]));
  }
  #browser();
  out <- cbind(m[,1:2],df,LLc,LLs,LR,chi2,npar,AIC={if(exists('AIC')) AIC else NULL},BIC={if(exists('BIC')) BIC else NULL},p,sig,chosen=m[,2]==chosen);
  return(out);
}
