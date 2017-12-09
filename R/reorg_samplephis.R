#' Trying again at samplephis from the  inside out
#' 
#' To work on these individually, load `data/testdata.rda` the following hopefully
#' self-explanatory objects currently exist:
#' `test_iisims`,`test_phis`,`test_pnlst`,`test_radii`,`test_tfresp`
#' 
#' # Helper functions
#' 
#' Lazy way to trigger a debug or save event in another session that's running
#' the main powertrip 
ptmg <- function(action=c('save','debug')
                 ,wd=paste0(getwd(),'/')
                 ,savetrigger=paste0(wd,'pt_savedata')
                 ,debugtrigger=paste0(wd,'pt_debug')){
  action <- match.arg(action,several.ok = T);
  if('save' %in% action) file.create(savetrigger);
  if('debug' %in% action) file.create(debugtrigger);
}

#' Shortcut for loading a saved powertrip environment
load.ptenv <- function(file='pt_result.rdata',env=new.env(),logenvonly=T,savewait=0.5,...) {
  if(savewait>0) ptmg(action='save',...);
  Sys.sleep(savewait);
  load(file,envir = env); class(env)<-c('ptenv',class(env)); 
  if(logenvonly) env$logenv else env;
}

#' Extract from a pt-derived dataframe a cartesian one with the selected 
#' columns as radius and phis
dfcrt <- function(data,radius,phis=c('phi.Var1','phi.Var2'),subset=T){
  data.frame(pol2crt(subset(data,subset=subset)[,c(radius,phis)]));
};

#' Set up lm fits and predictions in logenv to subsequently update
#' Note: this function intentionally blows away anything already in
#' logenv[[pathtop]] that may have the same name. For updates, use
#' env_fitupdt()
#

env_fitinit <- function(logenv
                        ,fields=alist(rad=ifelse(preds['conv',]==1,preds['radest',],NA)
                                      ,phi=phi
                                      #,convs=sum(preds['conv',])
                                      ,nsims=nsims)
                        #,pathtop='fits'
                        ,...){
  if(length(intersect(c('rad','phi','nsims'),names(fields)))<3){
    stop("The fields argument must be an alist that includes 'rad','phi',and 'nsims' among its elements, see 'args(env_fitinit)'");
  };
  logenv$fits$fields <- fields;
  #if(!pathtop %in% names(logenv)) logenv[[pathtop]] <- list();
  # hardcoding pathtop in hopes of making functions less brittle
  if(!'fits' %in% names(logenv)) logenv$fits <- list();
  # obtain the updated data we will need
  logenv$fits$radsphis <- pt2df(ptenv=logenv,fields=fields,...);
  # here is why we needed rads, phis, and nsims...
  with(logenv$fits,{
    logenv$fits$radnames <- radnames <- grep('^rad\\.',names(radsphis),val=T);
    logenv$fits$phinames <- phinames <- grep('^phi\\.',names(radsphis),val=T);
    # formula for modeling number of simulations needed to converge as a function 
    # of angle. This will get updated to create all the formulas for predicting
    # the radii for the respective ptpnl_ functions
    nsimsfrm<-update(as.formula(paste('nsims~('
                                      ,paste0(sprintf('%1$s+cos(%1$s)+sin(%1$s)',phinames)
                                              ,collapse='+'),')^2')),.~.);
    logenv$fits$frms <- frms <- c(nsims=nsimsfrm
                                ,sapply(radnames
                                        ,function(xx) update(nsimsfrm
                                                             ,as.formula(paste0(xx,'~.')))));
  logenv$fits$models <- sapply(frms,function(xx) {
      oo<-lm(formula=xx,logenv$fits$radsphis); oo$call$formula <- xx; oo;
      },simplify=F);
  });
};

#' To run each time a new set of phis has been completed
env_fitupdt <- function(logenv,...){
  logenv$fits$radsphis <- pt2df(ptenv = logenv,fields=logenv$fits$fields,...);
  logenv$fits$models <- sapply(logenv$fits$models,update,data=logenv$fits$radsphis);
};

#' To get predictions
env_fitpred <- function(logenv,newdata
                        # someday you'll forget what the hell was the structure
                        # of your input data for this particular run and when you
                        # do, set example=T and run, with logenv as the only other
                        # required argument
                        ,example=F
                        ,maxrad=NULL
                        ,...){
  if(example) return(summary(logenv$fits$radsphis[,logenv$fits$phinames]));
  if(missing(newdata)) newdata <- logenv$fits$radsphis[,logenv$fits$phinames];
  oo<-do.call(data.frame,sapply(logenv$fits$models
                            ,function(xx) with(predict(xx,newdata=newdata,se=T)
                                               ,cbind(fit=fit,se=se.fit)),simplify=F));
  if(!is.null(maxrad)){
    fnames<-paste0(logenv$fits$radnames,'.fit');
    snames<-paste0(logenv$fits$radnames,'.se');
    ex <- apply(oo[,fnames<-paste0(logenv$fits$radnames,'.fit'),drop=F],2,function(xx) xx<0|xx>maxrads);
    oo[fnames][ex] <- NA;
    oo[snames][ex] <- NA;
  }
  oo;
};

env_state <- function(logenv,coords='coords',summ='summ',fits='fits'
                      ,radsphis = 'radsphis'
                      ,fitslist = c('fields','frms','models','phinames','radnames')
                      ,...){
  if(!coords %in% names(logenv) ||
     length(logenv[[coords]]) <= 1 ||
     (ndata <- sum(sapply(logenv[[coords]],function(xx) summ %in% names(xx)))) <= 1
     ) return('needsdata') else {
       if(length(fnames<-names(logenv[[fits]]))==0 || 
          length(setdiff(fnames,c(fitslist,radsphis)))>0 ||
          is.null(rprows <- nrow(logenv[[fits]][[radsphis]]))
          ) return('needsinit') else {
            if(rprows != ndata) return('needsupdate') else {
              return('uptodate');
            }
          }
     }
}

#' Only the first argument is required if logenv is not initialized
#' 
make_phis <- function(logenv,npoints,maxs,mins,nphis,phiprefix='phi',bestfrac=0.5,numse=2,fresh=F,...){
  # Note: do not change phi prefix without also changing the fields argument of
  # pt2df() or the generation of subsequent rounds of phis will break
  # if no data obtained yet or fresh manually set to T then phiprefix and nphis
  # are required arguments' logenv and npoints are always required
  # so are maxs and mins until/unless they get saved in logenv as well...
  # bestfrac is the quantile (of the filtering criterion) above which to keep 
  # the candidate phis... in order to target the most informative and least
  # computationally expensive parts of the current parameter space
  switch(env_state(logenv,...)
         ,needsdata={fresh<-T}
         ,needsinit={env_fitinit(logenv)}
         ,needsupdate={env_fitupdt(logenv)});
  if(!fresh) {
    phinames<-logenv$fits$phinames;
    nphis<-length(phinames);
  } else {
    phinames <- paste0(phiprefix,seq_len(nphis));
  }
  oo<-data.frame(matrix(runif((nphis-1)*npoints,0,pi),nrow=npoints),runif(npoints,0,2*pi));
  colnames(oo) <- paste0(phiprefix,seq_len(nphis));
  maxrad<-apply(oo,1,pollim,maxs=maxs,mins=mins); 
  if(!fresh){
    fp <- env_fitpred(logenv,newdata = oo,maxrad = maxrad);
    snames <- paste0(logenv$fits$radnames,'.se');
    fnames <- paste0(logenv$fits$radnames,'.fit');
    if(bestfrac<1){
      filter <- rowMeans(apply(fp[,snames],2,rank,na.last = 'keep'),na.rm = T)/fp$nsims.fit;
      filterkeep <- filter>quantile(filter,bestfrac,na.rm=T);
      oo <- oo[filterkeep,];
      maxrad <- maxrad[filterkeep];
      fp <- fp[filterkeep,];
    }
    oo$mins <- pmax(apply(fp[fnames]-numse*fp[snames],1,min,na.rm=T),0);
    oo$maxs <- pmin(apply(fp[fnames]+numse*fp[snames],1,max,na.rm=T),maxrad);
  } else {oo$mins <- 0; oo$maxs <- maxrad};
  cbind(oo,maxrad);
}

#' Any arguments in ... are treated as expressions to evaluate in the context of
#' ptenv$coords[[XXX]][[summname]]
pt2df <- function(ptenv,summname='summ'
                  ,fields=alist(rad=preds['radest',],phi=unname(phi),conv=preds['conv',],nsims=nsims,time=time,maxrad=maxrad,lims=lims)
                  ,...){
  dots <- as.list(substitute(list(...))[-1]);
  # to the default fields below add in maxrad=maxrad next time logenv is rebuilt
  fields <- c(fields,dots);
  oo<-data.frame(t(sapply(ptenv$coords,function(xx) with(xx[[summname]][[1]],do.call('c',fields)))));
  cat('Read',nrow(oo),'rows,',ncol(oo),'columns.\n');
  oo;
}


#' # Core functions
#' 
resp_preds <- function(tfresp,radii,glmfit,saveglm=F,power=0.8){
  # think about what to do when this function fails
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

preds_lims <- function(preds,tol=0.01,limit=1e6,numse=2,...){
  # takes the result of iterating res_preds() over each column of tfresps returned
  # by the panel
  # determines which panels failed on latest round 
  notfailed <- preds['conv',]==1 & preds['radest',]>0 & preds['radest',] < limit;
  # determines which panels not yet finished::
  #     sepred*setol > tol  | min < 0 | max > limits
  notdone <- preds['respse',] > tol;
  # if all panels failed, catch failure and signal accordingly
  if(!any(notfailed)) return(c(min=NA,max=NA,status=-1,notfailed=notfailed,notdone=notdone)); #TODO: log failure
  # if all non-failed panels converged, return final estimates with failed ones flagged
  # if there are any non-failed non-converged panels
  # returns a single max and min for the next round of radii based on those
  pnstatus <- 
  if(status<-!any(select <- notfailed & notdone)) select <- notfailed;
  return(c(
     min=max(min(preds['radest',select] - preds['radse',select]*numse),0)
    ,max=min(max(preds['radest',select] + preds['radse',select]*numse),limit)
    ,status=status,notfailed=notfailed,notdone=notdone));
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

phi_radius <- function(phi,maxrad,pnlst
                       ,pneval=sapply(pnlst,attr,'eval')
                       ,pnfit=names(pneval)[pneval]
                       ,pninfo=names(pneval)[!pneval]
                       ,philabel=paste0('phi_',paste0(round(phi,3),collapse='_'))
                       ,logenv=logenv,nrads=20
                       ,max=maxrad,min=0,numse=1.5
                       ,ptsim=ptsim_nlin
                       ,wd=paste0(getwd(),'/')
                       # when a file having the name specified by this variable is found, 
                       # the dataenv,logenv, and errenv objects are saved
                       ,savetrigger=paste0(wd,'pt_savedata')
                       ,debugtrigger=paste0(wd,'pt_debug')
                       # name of file to which to save when savetrigger encountered
                       ,savefile=paste0(wd,'pt_result.rdata')
                       # name of file that will be sourced (once and then moved) if found
                       ,sourcepatch=paste0(wd,'pt_sourcepatch.R')
                       ,phicycle=0
                       ,...){
  # The function which will plug the above modules into each other and test them
  # jointly
  
  # for debugging
  phi_radius_env <- environment();
  on.exit(.GlobalEnv$phi_radius_env <- phi_radius_env);
  # end debugging 

  # First determine which functions in pnlst are evaluable (i.e retrun verdicts 
  # rather than just summary statistics)
  # obtain the maximum allowed radius for the current phis
  # (derived from the cartesian limits maxs and mins)
  #if(missing(maxrad)) maxrad<-pollim(phi,maxs=maxs,mins=mins);
  # lenght of current verdicts
  #nntf <- length(list_tfresp);
  cycle <- 1;
  lims <- c(min=min,max=max,status=0);
  list_tfresp <- list_radii <- list();
  tfoffset <- 0;
  t0 <- Sys.time();
  while(lims['status']==0){
    list_radii[[cycle]] <- runif(nrads,lims['min'],lims['max']);
    for(ii in 1:nrads){
      iicoords <- c(list_radii[[cycle]][ii],phi);
      iidat <- ptsim(pol2crt(iicoords));
      list_tfresp[[tfoffset+ii]] <- sapply(pnlst[pneval],function(xx) any(xx(iidat,iicoords)));
      if(any(is.na(list_tfresp[[tfoffset+ii]]))) list_radii[[cycle]][ii] <- NA;
    }
    # then fit models on the panel verdicts (T/F), tfresp
    # na.omits might be unnecessary
    testtf<-data.frame(do.call(rbind,list_tfresp));
    testrd<-unlist(list_radii);
    #if(length(testrd)!=nrow(testtf)) browser();
    preds <- sapply(testtf,resp_preds,radii=testrd);
    #preds<-sapply(na.omit(data.frame(do.call(rbind,list_tfresp))),resp_preds,radii=na.omit(unlist(list_radii)));
    lims <- preds_lims(preds,limit=maxrad,numse = numse);
    #if(any(na.omit(lims[c('min','max')])<0|na.omit(lims[c('min','max')])>maxrad)) browser(text='Invalid limits!');
    # Currently, we give up on this entire set of phis and exit from phi_radius()
    # the first cycle when glm fails for all panel functions regardless of why.
    # TODO: this doesn't go far enough-- find a way to detect > X% NAs in the 
    #       verdicts and also disqualify. Perhaps within the panel function
    # TODO: this may run too long-- set a cycle-limit exceeding which would also 
    #       end the attempt to estimate the radius for these phis
    # TODO: implement logging at the phi level-- don't retain model results or 
    #       stats about the current simulated dataset, just the global part of 
    #       ptpnl_summary especially the phi
    # TODO: implement error logging for ptpnl_ functions
    # TODO: in normal result logging for ptpnl_ functions retain the verdict as well
    # DONE: consider not bothering with na.omit since glm probably does it
    #       internally anyway. Makes it easier to find fraction NA too (above TODO)
    cycle <- cycle+1; tfoffset <- length(list_tfresp);
    if(file.exists(savetrigger)) {
      save(phi_radius_env,logenv,file=savefile);
      file.remove(savetrigger);
    }
    if(file.exists(debugtrigger)){
      browser();
    }
    if(file.exists(sourcepatch)) {
      source(sourcepatch,local = T);
      file.rename(sourcepatch,paste0(sourcepatch,'.bak'));
    }
  };
  if(lims['status']==1){
    cat('Success: ');
    for(pp in pnfit) {
      ppdat <-ptsim(pol2crt(c(preds['radest',pp],phi))); 
      pnlst[[pp]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,pp));
      # for each pp (verdict-returning panel function) we generate a separate dataset therefore
      # we need to iterate over all the summary-only non-verdict functions for each of these datasets
      for(qq in pninfo) pnlst[[qq]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,qq),time=Sys.time()-t0);
    };
  } else cat('Failure: ');
  cat('radii= ',try(preds['radest',]),'\tphi= ',try(phi),'\tlims= ',c(maxrad,lims[-3]),'\n');
  # DONE: add back in the dynamic script execution and the external exit directive
  # DONE: add a phi -> radius prediction step (multivariate, whole parameter space)
  # DONE: for prioritizing phis, also model the runtime to get most uncertainty
  #       per second of runtime or per simulation
  # TODO: also, model the 'dead-zones' -- places where we had to give up --
  #       and exclude them
  # DONE?: figure out why negative radii are being allowed and fix
  # DONE?: figure out why plots look so wierd... is pol2crt wrong?
  # TODO: finalize outer function
  # DONE?: launch the linear model version
  # TODO: try to resurrect simsurve and survwrapper
}

test_harness<-function(logenv=logenv
                       #,maxs=c(2,4.5,6),mins=c(-3.1,-1.3,-6)
                       ,maxs=c(20,20,20),mins=c(-20,-20,-20)
                       ,npoints=50,nphis=2,nrads=20,numse=2
                       # which fraction of the most impactful phis should we 
                       # model each time?
                       ,wd=paste0(getwd(),'/'),savetrigger=paste0(wd,'pt_savedata')
                       ,debugtrigger=paste0(wd,'pt_debug')
                       ,savefile=paste0(wd,'pt_result.rdata')
                       ,sourcepatch=paste0(wd,'pt_sourcepatch.R')
                       ,bestfrac = 0.5
                       ,pnlst=list(#lm=ptpnl_lm
                                   lm2=update(ptpnl_lm,fname='lm2',frm=yy~(.)*group)
                                   ,tt=ptpnl_tt
                                   ,wx=ptpnl_wx
                                   ,summ=ptpnl_summary)
                       ,ptsim=ptsim_nlin
                       # technically this can run forever continuously generating
                       # better predictions and get stopped from an external R
                       # session via ptmg()... but it feels wrong to not install
                       # an automatic off-switch, so we will have maxphicycle
                       # if we do need literally indefinite runtime just set 
                       # this to Inf
                       ,maxphicycle=1e6,...){
  #phis <- cbind(matrix(runif(npoints*(nphi-1),0,pi),nrow=npoints,ncol=nphi-1),runif(npoints,0,2*pi));
  # trying more regularly spaced phis, to get a handle on the weird output
  # 10 intervals per phi results in 10,000 distinct points
  #phiseq <- sample(seq(0,pi,length.out=npoints),npoints,rep=F);
  #browser();
  #phis <- as.matrix(do.call(expand.grid,c(replicate(nphi-1,phiseq,simplify=F),list(sample(seq(0,2*pi,length.out = npoints),npoints,rep=F)))));
  pneval_ <- sapply(pnlst,attr,'eval');
  pnfit <- names(pneval_)[pneval_];
  pninfo <- names(pneval_)[!pneval_];
  phicycle <- 1; 
  while(phicycle < maxphicycle){
    phis <- make_phis(logenv=logenv,npoints = npoints,maxs=maxs,mins=mins
                      ,nphis = nphis,numse = numse);
    #logenv$temp$maxrads <- maxrads <- apply(phis,1,pollim,maxs=maxs,mins=mins);
    actualpoints <- nrow(phis);
    for(ii in seq_len(actualpoints)){
      cat(ii,'\t|');
      phi_radius(phi=unlist(phis[ii,seq_len(nphis)]),maxrad=phis[ii,'maxrad'],pnlst=pnlst
                 ,pneval=pneval_,pnfit=pnfit,pninfo=pninfo
                 ,logenv=logenv,max=phis[ii,'maxs'],min=phis[ii,'mins'],nrads=nrads
                 ,numse=numse
                 ,ptsim=ptsim
                 ,wd=wd
                 ,savetrigger=savetrigger
                 ,debugtrigger=debugtrigger
                 ,savefile=savefile
                 ,sourcepatch=sourcepatch
                 ,phicycle=phicycle);
    }
    phicycle<-phicycle+1
  }
  data.frame(t(sapply(logenv$coords,function(xx) with(xx$summ[[1]],c(
    phi=phi,radii=ifelse(preds['conv',]==1,preds['radest',],NA)
    ,time=time,nsims=nsims,cycle=cycle)))));
  browser();
};
# Timing stopped at: 6573 5411 6740 ... about 1.87 hours to try 1000 pairs of phis
