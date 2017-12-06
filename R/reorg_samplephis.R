#' Trying again at samplephis from the  inside out
#' 
#' To work on these individually, load `data/testdata.rda` the following hopefully
#' self-explanatory objects currently exist:
#' `test_iisims`,`test_phis`,`test_pnlst`,`test_radii`,`test_tfresp`
#' 
resp_preds <- function(tfresp,radii,glmfit,saveglm=F,power=0.8){
  if(missing(glmfit)) {glmfit <- glm(tfresp~radii,family='binomial') } else {
    glmfit <- update(glmfit) };
  out <- MASS::dose.p(glmfit,p=power);
  out <- setNames(c(out[1],attr(out,'SE')
           ,with(predict(glmfit,newdata=data.frame(radii=out[1]),type='resp',se.fit=T)
                 ,c(fit[1],se.fit[1]))
           ,glmfit$converged)
           ,c('radest','radse','respest','respse','conv'));
  if(saveglm) attr(out,'glmfit') <- glmfit;
  out;
  # fits glmfit object to tfresp ~ radius or updates optionally passed glmfit
  # returns vector that contains reverse predictor for radius, se for reverse
  # predictor, and se of prediction of glmfit when evaluated at reverse predictor
}

#' given a vector of radii and data.frame or matrix of responses with the same 
#' number of columns as the length of radii, res_preds() should iterate over each
#' column

preds_lims <- function(preds,tolse=1,tol=0.01,limit=Inf,radci=2,...){
  # takes the result of iterating res_preds() over each column of tfresps returned
  # by the panel
  # determines which panels failed on latest round 
  notfailed <- preds['conv',]==1 & preds['radest',]>0 & preds['radest',] < limit;
  # determines which panels not yet finished::
  #     sepred*setol > tol  | min < 0 | max > limits
  notdone <- preds['respse',]*tolse > tol;
  # if all panels failed, catch failure and signal accordingly
  if(!any(notfailed)) return(c(min=NA,max=NA,status=-1)); #TODO: log failure
  # if all non-failed panels converged, return final estimates with failed ones flagged
  # if there are any non-failed non-converged panels
  # returns a single max and min for the next round of radii based on those
  if(status<-!any(select <- notfailed & notdone)) select <- notfailed;
  return(c(
     min=max(min(preds['radest',select] - preds['radse',select]*radci),0)
    ,max=min(max(preds['radest',select] + preds['radse',select]*radci),limit)
    ,status=status));
}

#' Not needed, one liner: runif(nn,lims[1],lims[2]);
#lims_radii <- function(lims,phis,opts,...){
  # using limits as provided by preds_lims() and a set of phis, uniformly sample
  # radii within those limits
  # returns a vector of radii
#}

radii_res <- function(radii,phis,opts,ptsim,pnlst,...){
  # using radii as provided by lims_radii() and phis, a vector of angles
  # simulates data at each value of radii and the entire vector of phis using
  # the ptsim() function and iterates over all the functions in pnlst()
  # returns a data structure which for each value of radii has the verdict returned
  # by each function of pnlst() (T/F)
}

testenv <- new.env();
load('data/testdata.rda',envir = testenv);

phi_radius <- function(phi=c(2.3,5.12)
                         ,philabel=paste0('phi_',paste0(round(phi,3),collapse='_'))
                         ,maxs=c(2,4.5,3),mins=c(-3.1,-1.3,-0.5)
                         ,nrads=20
                         ,pnlst=list(lm=ptpnl_lm,lm2=update(ptpnl_lm,fname="lm2",frm=yy~(.)*group),summ=ptpnl_summary)
                         ,ptsim=ptsim_nlin
                         ,...){
  # The function which will plug the above modules into each other and test them
  # jointly
  phi_radius_env <- environment();
  on.exit(.GlobalEnv$phi_radius_env <- phi_radius_env);
  # First determine which functions in pnlst are evaluable (i.e retrun verdicts 
  # rather than just summary statistics)
  pneval <- sapply(pnlst,attr,'eval');
  pnfit <- names(pneval)[pneval];
  pninfo <- names(pneval)[!pneval];
  # obtain the maximum allowed radius for the current phis
  # (derived from the cartesian limits maxs and mins)
  maxrad<-pollim(phi,maxs=maxs,mins=mins);
  # lenght of current verdicts
  #nntf <- length(list_tfresp);
  cycle <- 1;
  lims <- c(min=0,max=maxrad,status=0);
  list_tfresp <- list_radii <- list();
  tfoffset <- 0;
  while(lims['status']==0){
    list_radii[[cycle]] <- runif(nrads,lims['min'],lims['max']);
    for(ii in 1:nrads){
      iicoords <- c(list_radii[[cycle]][ii],phi);
      iidat <- ptsim(pol2crt(iicoords));
      list_tfresp[[tfoffset+ii]] <- sapply(pnlst[pneval],function(xx) any(xx(iidat,iicoords)));
      if(any(is.na(list_tfresp[[tfoffset+ii]]))) list_radii[[cycle]][ii] <- NA;
    }
    # then fit models on the panel verdicts (T/F), tfresp
    testtf<-na.omit(data.frame(do.call(rbind,list_tfresp)));
    testrd<-na.omit(unlist(list_radii));
    if(length(testrd)!=nrow(testtf)) browser();
    preds <- sapply(testtf,resp_preds,radii=testrd);
    #preds<-sapply(na.omit(data.frame(do.call(rbind,list_tfresp))),resp_preds,radii=na.omit(unlist(list_radii)));
    lims <- preds_lims(preds,limit=maxrad);
    # Currently, we give up on this entire set of phis and exit from phi_radius()
    # the first cycle when glm fails for all panel functions regardless of why.
    # TODO: this doesn't go far enough-- find a way to detect > X% NAs in the 
    #       verdicts and also disqualify. Perhaps within the panel function
    # TODO: this may run too long-- set a cycle-limit exceeding which would also 
    #       end the attempt to estimate the radius for these phis
    # TODO: implement logging at the phi level-- don't retain model results or 
    #       stats about the current simulated dataset, just the global part of 
    #       ptpnl_summary especially the phi
    # TODO: consider not bothering with na.omit since glm probably does it
    #       internally anyway. Makes it easier to find fraction NA too...
    cycle <- cycle+1; tfoffset <- length(list_tfresp);
  };
  cat(cycle,'|');
  if(lims['status']==1){
    cat('Success: ');
    for(pp in pnfit) {
      ppdat <-ptsim(pol2crt(c(preds['radest',pp],phi))); 
      pnlst[[pp]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,pp));
      # for each pp (verdict-returning panel function) we generate a separate dataset therefore
      # we need to iterate over all the summary-only non-verdict functions for each of these datasets
      for(qq in pninfo) pnlst[[qq]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,qq));
    };
  } else cat('Failure: ');
  cat('radii= ',try(preds['radest',]),'\tphi= ',try(phi),'\n');
  # TODO: add a phi -> radius prediction step (multivariate, whole parameter space)
  # TODO: run all the way through
  # TODO: finalize outer function
  # TODO: launch the linear model version
  # TODO: try to resurrect simsurve and survwrapper
}

test_harness<-function(logenv=logenv,maxs=c(2,4.5,6),mins=c(-3.1,-1.3,-6)
                       ,npoints=50,nphi=2,nrads=20
                       ,pnlst=list(lm=ptpnl_lm,lm2=update(ptpnl_lm,fname='lm2',frm=yy~(.)*group),summ=ptpnl_summary)
                       ,ptsim=ptsim_nlin,...){
  phis <- matrix(runif(npoints*nphi,0,2*pi),nrow=npoints,ncol=nphi);
  for(ii in 1:npoints){
    phi_radius(phi=phis[ii,],maxs=maxs,mins=mins,nrads=nrads,pnlist=pnlist,ptsim=ptsim);
  }
  browser();
  data.frame(t(sapply(logenv$coords,function(xx) with(xx$summ[[1]],c(
    phi=phi,radii=ifelse(preds['conv',]==1,preds['radest',],NA)
    ,time=time,nsims=nsims,cycle=cycle)))));
};
