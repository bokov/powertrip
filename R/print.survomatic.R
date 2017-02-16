
print.survomatic <- function(xx,...){
  # get the available models
  models <- names(xx$par.differences);
  out <- list();
  for(ii in models){
    ivars <-  modelinfo(ii,'vars');
    iout <- xx$par.differences[[ii]][-1,c('a','b','c','s','LR','df','chi2','p')];
    itested <- apply(iout[,c('a','b','c','s')],1,function(jj) paste(c('a','b','c','s')[!jj],collapse=""));
    iout$Yse <- iout$Yest <- iout$Xse <- iout$Xest <- NA;
    for(jj in ivars) {
      iout[itested==jj,c('Xest','Yest')] <- xx$par.differences[[ii]][1,paste(jj,1:2,sep="")];
      iout[itested==jj,'Xse'] <- xx$x.m[xx$x.m$model==ii,paste(jj,"se",sep=".")];
      iout[itested==jj,'Yse'] <- xx$y.m[xx$y.m$model==ii,paste(jj,"se",sep=".")];
    }
    out[[ii]] <- iout[,c('Xest','Xse','Yest','Yse','LR','df','chi2','p')];
  }
  print(out);
}

summary.survomatic <- print.survomatic;
