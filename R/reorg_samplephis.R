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
dfcrt <- function(data,radius,phis=c('phi1','phi2')
                  ,refcoords=0,transform=identity,subset=T){
  subset<-substitute(subset);
  oo<-data.frame(pol2crt(subset(data,subset=eval(subset))[,c(radius,phis)]));
  if(any(refcoords!=0)){
    if(length(refcoords)==1) refcoords <- rep_len(refcoords,length(phis)+1);
    oo <- oo + rbind(refcoords)[rep_len(1L,nrow(oo)),]
  }
  oo[is.infinite(as.matrix(oo))]<-NA;
  invisible(transform(oo));
};

# required: maxs, mins, if phis missing then 
make_polgrid <- function(maxs,mins,phis,nticks=5,refcoords=0*maxs
                         ,cartesian=T,radticks=0
                         # maxes and mins relative to refcoords
                         ,relmaxs=maxs-refcoords,relmins=mins-refcoords,...){
  # TODO: more error checking
  if(length(relmaxs)!=length(relmins)) stop('The maxs and mins arguments should be the same length');
  nphis<-length(relmaxs)-1;
  if(missing(phis)) {
    phiseq <- seq(0,pi,length.out = nticks);
    phis<-do.call(expand.grid,c(replicate(nphis-1,phiseq,simplify=F),list(seq(0,2*pi,length.out = 40*nticks))));
  } else if(nphis != ncol(phis)) stop('The length of both maxs and mins arguments should equal the number of columns in the phis argument, if used.');
  phinames <- names(phis);
  phis$rad <- apply(phis,1,pollim,maxs=relmaxs,mins=relmins);
  phis$ismax<-1;
  # for each point on the outer boundaries draw a line from the center that goes to it
  # and includes radticks points including the center and outer point
  if(radticks>0){
    phis<- data.frame(do.call(rbind,lapply(split(phis,phis[,phinames]),function(xx) {
      oo<-v2mat(xx,radticks);
      # we oimit the zero point to avoid redundancy
      oo[,'rad']<-seq(0,oo[1,'rad'],len=radticks+1)[-1];
      oo[-radticks,'ismax']<-0;
      oo})));
    # this gives us the one zero point we need
    phis <- rbind(phis,0);
    }
  if(cartesian) phis <- cbind(phis,dfcrt(phis,'rad',phinames,refcoords,...));
  #dfcrt(data=phis,radius = 'maxrad',phis=colnames(phis),refcoords=refcoords,...)
}

#' makena takes a prefectly good object and ruins it by recursively turning everything
#' in it into an NA! Why would anybody want to do that? Well one reason might be
#' to create placeholder objects to return when an error is encountered that prevents
#' the correct response from being generated-- placeholders with, as far as I can
#' tell at this time, identical dimensions, sub-objects, names, everything.
makena <- function(xx,fn){if(missing(fn)) fn<-sys.function(); if(is.list(xx)) xx[]<-lapply(xx,fn,fn=fn) else is.na(xx)<-T; xx}

#' v2mat() stacks a vector on itself to make a matrix with the same number of 
#' columns as the length of the vector and nr rows. Because for some reason I
#' keep having to create matrices like that for this project
#' The extra `as.matrix()` is needed here to properly handle data.frames
v2mat <- function(vv,nr) matrix(as.matrix(vv),nrow=nr,ncol=length(vv),byrow=T,dimnames=list(NULL,colnames(vv)));

#' Set up lm fits and predictions in logenv to subsequently update
#' Note: this function intentionally blows away anything already in
#' logenv[[pathtop]] that may have the same name. For updates, use
#' env_fitupdt()
#

env_fitinit <- function(logenv,...){
  #if(length(intersect(c('rad','phi','nsims'),names(fields)))<3){
  #  stop("The fields argument must be an alist that includes 'rad','phi',and 'nsims' among its elements, see 'args(env_fitinit)'");
  #};
  #logenv$fits$fields <- fields <- alist(rad=ifelse());
  #if(!pathtop %in% names(logenv)) logenv[[pathtop]] <- list();
  # hardcoding pathtop in hopes of making functions less brittle
  if(!all(c('phinames','notfailed','fields','radnames') %in% names(logenv$names))) stop('The environment passed through the logenv variable needs to have its names initialized');
  if(!'fits' %in% names(logenv)) logenv$fits <- list();
  # obtain the updated data we will need
  logenv$fits$radsphis <- pt2df(ptenv=logenv,fields=logenv$names$fields);
  # here is why we needed rads, phis, and nsims...
  with(logenv$fits,{
    #logenv$fits$radnames <- radnames <- grep('^rad\\.',names(radsphis),val=T);
    #logenv$fits$phinames <- phinames <- grep("^phi[0-9]+$",names(logenv$fits$radsphis),val=T);
    # formula for modeling number of simulations needed to converge as a function 
    # of angle. This will get updated to create all the formulas for predicting
    # the radii for the respective ptpnl_ functions
    nsimsfrm<-update(as.formula(paste('nsims~('
                                      ,paste0(sprintf('%1$s+cos(%1$s)+sin(%1$s)',logenv$names$phinames)
                                              ,collapse='+'),')^2')),.~.);
    logenv$fits$frms <- frms <- c(nsims=nsimsfrm
                                ,sapply(logenv$names$radnames
                                        ,function(xx) update(nsimsfrm
                                                             ,as.formula(paste0(xx,'~.')))));
  logenv$fits$models <- sapply(frms,function(xx) {
      oo<-lm(formula=xx,logenv$fits$radsphis); oo$call$formula <- xx; oo;
      },simplify=F);
  });
};

#' To run each time a new set of phis has been completed
env_fitupdt <- function(logenv,...){
  logenv$fits$radsphis <- pt2df(ptenv = logenv,fields=logenv$names$fields,...);
  logenv$fits$models <- sapply(logenv$fits$models,update,data=logenv$fits$radsphis,simplify=F);
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
  if(is.null(phinames<-logenv$names$phinames)||is.null(snames<-logenv$names$snames)||is.null(fnames<-logenv$names$fnames)){
    stop('logenv must contain valid phinames, snames, and fnames vectors inside its names list so we know which columns to use');
  }
  if(example) return(summary(logenv$fits$radsphis[,phinames]));
  if(missing(newdata)) newdata <- logenv$fits$radsphis[,phinames];
  oo<-do.call(data.frame,sapply(logenv$fits$models
                            ,function(xx) with(predict(xx,newdata=newdata,se=T)
                                               ,cbind(fit=fit,se=se.fit)),simplify=F));
  if(!is.null(maxrad)){
    #fnames<-paste0(logenv$fits$radnames,'.fit');
    #snames<-paste0(logenv$fits$radnames,'.se');
    #ex <- apply(oo[,fnames<-paste0(logenv$fits$radnames,'.fit'),drop=F],2,function(xx) xx<0|xx>maxrad);
    ex <- apply(oo[,fnames,drop=F],2,function(xx) xx>maxrad);
    oo[fnames][ex] <- NA;
    oo[snames][ex] <- NA;
    # apparently especially with small samples you can get negative nsims!
    oo$nsims.fit[oo$nsims.fit<0] <- NA;
  }
  oo;
};

env_state <- function(logenv,coords='coords',summ='summ',fits='fits'
                      ,radsphis = 'radsphis'
                      ,fitslist = c('fields','frms','models','phinames','radnames')
                      ,...){
  if(!coords %in% names(logenv) ||
     # below 30 results we get wierd lm behavior given the number of interaction
     # terms we have. The absolute hard limit is 20!
     length(logenv[[coords]]) <= 30 ||
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

#' nphis is only required if logenv has no data yet
#' logenv, npoints, maxs, mins are always required
#' 
#' Why did successful convergence rate go up so dramatically after implementing 
#' this? Possibly because we have pre-empted some bug in the old limit 
#' initialization that was being done inside phi_radius() and now we are
#' properly enforcing valid radii limits. But it may also be that increasing 
#' the number of radii (to 50) and decreasing the bounding box (from 20's to 2's
#' with one 0.5) increases the chances that the first set of radii will give 
#' enough of a starting result for glm to work with.
#' 
#' Definitely looks like setting limits such that they actually get encountered
#' causes a lot of failures
make_phis <- function(logenv,npoints,maxs,mins,phiprefix='phi',bestfrac=0.5,numse=2,fresh=F,...){
  # Note: do not change phi prefix without also changing the fields argument of
  # pt2df() or the generation of subsequent rounds of phis will break
  # if no data obtained yet or fresh manually set to T then phiprefix and nphis
  # are required arguments' logenv and npoints are always required
  # so are maxs and mins until/unless they get saved in logenv as well...
  # bestfrac is the quantile (of the filtering criterion) above which to keep 
  # the candidate phis... in order to target the most informative and least
  # computationally expensive parts of the current parameter space
  if((nphis<-length(phinames<-logenv$names$phinames))==0||
     is.null(snames<-logenv$names$snames)||
     is.null(fnames<-logenv$names$fnames)){
    stop('logenv must contain valid phinames, snames, and fnames vectors inside its names list so we know which columns to use');
  }
  switch(env_state(logenv,...)
         ,needsdata={print('Need to accumulate data pre initialization.');fresh<-T}
         ,needsinit={print('Initializing logenv$fits.');env_fitinit(logenv)}
         ,needsupdate={print('Updating logenv$fits');env_fitupdt(logenv)});
  #nphis<-length(phinames <- logenv$names$phinames);
  #oo<-data.frame(matrix(runif((nphis-1)*npoints,0,pi),nrow=npoints),runif(npoints,0,2*pi));
  # I wonder if having the first dimension be the 0-2pi one will improve coverage? Or maybe blow things up again 
  # and waste valuable time, save for later
  #oo<-data.frame(runif(npoints,0,2*pi),matrix(runif((nphis-1)*npoints,0,pi),nrow=npoints));
  # changing from uniform cartesian below to normal cartesian let's see...
  #oo<-data.frame(crt2pol(gencartlims(maxs,mins,round(npoints/nphis+1))));
  oo <- data.frame(crt2pol(gencartnorm(maxs,mins,npoints)));
  colnames(oo) <- c('maxrad',phinames); #paste0(phiprefix,seq_len(nphis));
  # temporary note: this was hotpatched on phicycle 46 of local instance
  oo[,'maxrad'] <- apply(oo[,-1],1,pollim,maxs=maxs,mins=mins); 
  if(!fresh){
    #browser();
    fp <- env_fitpred(logenv,newdata = oo,maxrad = oo$maxrad);
    #snames <- logenv$names$snames; fnames <- logenv$names$fnames;
    bestfrac<- Inf #disabling this thing, it might be part of the problem
    if(bestfrac<1){
      filter <- rowMeans(apply(fp[,snames],2,rank,na.last = 'keep'),na.rm = T); #/fp$nsims.fit;
      filterkeep <- filter>quantile(filter,bestfrac,na.rm=T);
      # in env_fitpred() we turn illegal nsims.fit values (<0) into NAs so need
      # to not make that break the filtering...
      filterkeep <- filterkeep & !is.na(filterkeep);
      if((addback<-bestfrac*npoints-sum(filterkeep))>0) filterkeep[!filterkeep][seq_len(addback)]<-T;
      oo <- oo[filterkeep,];
      #maxrad <- maxrad[filterkeep];
      fp <- fp[filterkeep,];
    }
    # because of the addback hack above, we now have phis with NA-only radius
    # predictions getting added back in. That's okay, if they are NA they are 
    # meaningless anyway, so replacing NAs with 0s or maxrad values
    #oo$mins <- pmax(apply(fp[fnames]-numse*fp[snames],1,min,na.rm=T),0);
    oo$mins <- pmax(apply(fp[,fnames,drop=F]-numse*fp[,snames,drop=F],1,min,na.rm=T),0);
    oo$mins[is.infinite(oo$mins)]<- 0; #0; TODO: properly handle this case instead of hardcoded value which will break for non-log cases
    oo$maxs <- pmin(apply(fp[,fnames,drop=F]+numse*fp[,snames,drop=F],1,max,na.rm=T),oo$maxrad);
    # TODO: catch the maxs<0 case further upstream
    badmaxs<-oo$maxs<=0|is.infinite(oo$maxs);
    oo$maxs[badmaxs]<-oo$maxrad[badmaxs];
  } else {oo$mins <- 0; oo$maxs <- oo$maxrad};
  oo; #cbind(oo,maxrad);
}

#' Any arguments in ... are treated as expressions to evaluate in the context of
#' ptenv$coords[[XXX]][[summname]]
pt2df <- function(ptenv,summname='summ'
                  ,fields=alist(rad=preds['radest',],phi=unname(phi)
                                ,conv=preds['conv',],maxrad=maxrad,lims=lims
                                ,nsims=nsims,time=time
                                ,cycle=cycle,phicycle=phicycle)
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
#' tfresp is a vector of T/F, radii is a numeric vector
#' power is the target hit-rate (i.e. we are predicting the radius at which 
#' tfresp should be T power fraction of the time) and narate is what fraction
#' of the values in tfresp may be NA before this is treated as an error
resp_preds <- function(tfresp,radii,power=0.8,narate=0.5){
  # think about what to do when this function fails
  #if(missing(glmfit)) {glmfit <- glm(tfresp~radii,family='binomial') } else {
  #  glmfit <- update(glmfit) };
  # system.time shows that creating the fit de-novo is faster than updating
  outnames <- c('radest','radse','respect','respse','conv');
  outerror <- setNames(c(rep(NA,length(outnames)-1),-1),outnames);
  okay <- !is(try(nas<-mean(is.na(tfresp))),'try-error');
  if(okay & nas < narate) {
    okay <- is(glmfit <- try(glm(tfresp~radii,family='binomial')),'glm');
  } else return(outerror);
  if(okay) {
    okay <- is(out <- try(MASS::dose.p(glmfit,p=power)),'glm.dose');
  } else return(outerror);
  if(okay){
    return(setNames(c(out[1],attr(out,'SE')
                      ,with(predict(glmfit,newdata=data.frame(radii=out[1])
                                    ,type='resp',se.fit=T),c(fit[1],se.fit[1]))
                      ,glmfit$converged),outnames));
    } else return(outerr);
  # fits glmfit object to tfresp ~ radius or updates optionally passed glmfit
  # returns vector that contains reverse predictor for radius, se for reverse
  # predictor, and se of prediction of glmfit when evaluated at reverse predictor
}

#' given a vector of radii and data.frame or matrix of responses with the same 
#' number of columns as the length of radii, res_preds() should iterate over each
#' column

# used to have tol at 0.01
preds_lims <- function(preds,tol=1,limit=1e6,numse=2,...){
  # takes the result of iterating res_preds() over each column of tfresps returned
  # by the panel
  # determines which panels failed on latest round 
  # even one the log scale negative values should still be disallowed, so 
  # uncommenting the next line and dcommenting out the one I had for a while
  notfailed <- preds['conv',]==1 & preds['radest',]>0 & preds['radest',] < limit;
  #notfailed <- preds['conv',]==1 & preds['radest',] < limit;
  # determines which panels not yet finished::
  #     sepred*setol > tol  | min < 0 | max > limits
  notdone <- preds['respse',] > tol;
  # if all panels failed, catch failure and signal accordingly
  if(!any(notfailed,na.rm=T)) return(c(min=NA,max=NA,status=-1,notfailed=notfailed,notdone=notdone)); #TODO: log failure
  # if all non-failed panels converged, return final estimates with failed ones flagged
  # if there are any non-failed non-converged panels
  # returns a single max and min for the next round of radii based on those
  #pnstatus <- 
  status<-try(!any(select <- notfailed & notdone));
  if(is(status,'try-error')) browser();
  if(status) select <- notfailed;
  return(c(
     min=max(min(preds['radest',select] - preds['radse',select]*numse),0) # used to be 0
    ,max=min(max(preds['radest',select] + preds['radse',select]*numse),limit)
    ,status=status,notfailed=notfailed,notdone=notdone));
}

#' Not needed, one liner: runif(nn,lims[1],lims[2]);
#lims_radii <- function(lims,phis,opts,...){
  # using limits as provided by preds_lims() and a set of phis, uniformly sample
  # radii within those limits
  # returns a vector of radii
#}

#radii_res <- function(radii,phis,opts,ptsim,pnlst,...){
  # using radii as provided by lims_radii() and phis, a vector of angles
  # simulates data at each value of radii and the entire vector of phis using
  # the ptsim() function and iterates over all the functions in pnlst()
  # returns a data structure which for each value of radii has the verdict returned
  # by each function of pnlst() (T/F)
#}

#testenv <- new.env();
#load('data/testdata.rda',envir = testenv);

#' A function to take a series of pre-generated coordinates, do simulations, and 
#' run the panel on each one. This is different from phi_radius, below, because 
#' that function operates on set of phis at a time, generates its own radii, and
#' keeps cycling until it either converges on optimal radii or experiences a
#' decisive failure. Also unlike phi_radius(), the coordinates are expected to
#' already be on cartesian scale
#' 
#' The purpose of powergrid() at the moment is to experiment with various coordinate
#' sampling algorithms and maxs/mins settings, etc. in order to find a region in
#' the parameter space where the ptpnl_ and ptsim functions at least work at all
#' before embarking on the harder task of finding optima
powergrid <- function(coords,refcoords,ptsim,pnlst
                      ,phicols,cartcols
                      ,pneval=sapply(pnlst,attr,which='eval')
                      ,backtrans=identity,nn=250,logenv=logenv,...){
  # create coords to pass to ptsim from the CARTESIAN part
  iis<-seq_len(nrow(btcoords<-backtrans(as.matrix(coords[,cartcols]))));
  # create philabels from the POLAR part
  philabels<-paste0('phi_',apply(round(coords[,phicols],3),1,paste0,collapse='_'));
  for(ii in iis){
    iidat <- ptsim(btcoords[ii,],nn=nn,refcoords=refcoords,...);
    sapply(pnlst[pneval],function(xx) any(xx(iidat,btcoords[ii,],logenv,philabel=philabels[ii],...)));
  }
};

phi_radius <- function(phi,maxrad,pnlst,pnlph,refcoords
                       ,pneval=sapply(pnlst,attr,which='eval')
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
                       ,phicycle=0,backtrans=identity
                       ,...){
  # The function which will plug the above modules into each other and test them
  # jointly
  
  # for debugging
  logenv$state$phi_radius <- environment();
  #on.exit(.GlobalEnv$phi_radius_env <- phi_radius_env);
  #first_fail <- first_success <- T;
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
  btrefcoords <- backtrans(refcoords);
  list_tfresp <- list_radii <- list();
  tfoffset <- 0;
  t0 <- Sys.time();
  while(lims['status']==0){
    # if the coords have been transformed, backtrans puts them back on the scale
    # that ptsim() expects
    cyclecoords <- backtrans(
      pol2crt(cbind(
        # save untransformed randomly generated radians for later
        # if we sample a root of a runif, less biased toward small distances
        list_radii[[cycle]] <- runif(nrads,lims['min'],lims['max'])^0.75
        # turn the static coordinate vector into matrix with one column for each
        # phi and nrads rows
        ,rbind(phi)[rep_len(1,nrads),])));
    for(ii in 1:nrads){
      # TODO: pre-calculate lcoords, lvars, and refcoords and pass to ptsim
      iidat <- ptsim(cyclecoords[ii,],refcoords=btrefcoords,...);
      list_tfresp[[tfoffset+ii]] <- sapply(pnlst[pneval],function(xx) any(xx(iidat,iicoords)));
      #if(any(is.na(list_tfresp[[tfoffset+ii]]))) list_radii[[cycle]][ii] <- NA;
    }
    # then fit models on the panel verdicts (T/F), tfresp
    # na.omits might be unnecessary
    testtf<-data.frame(do.call(rbind,list_tfresp));
    testrd<-unlist(list_radii);
    #if(length(testrd)!=nrow(testtf)) browser();
    preds <- try(sapply(testtf,resp_preds,radii=testrd));
    if(class(preds)[1]=='try-error') browser();
    #preds<-sapply(na.omit(data.frame(do.call(rbind,list_tfresp))),resp_preds,radii=na.omit(unlist(list_radii)));
    lims <- preds_lims(preds,limit=maxrad,numse = numse);
    #if(any(na.omit(lims[c('min','max')])<0|na.omit(lims[c('min','max')])>maxrad)) browser(text='Invalid limits!');
    # Currently, we give up on this entire set of phis and exit from phi_radius()
    # the first cycle when glm fails for all panel functions regardless of why.
    # TODO: this doesn't go far enough-- find a way to detect > X% NAs in the 
    #       verdicts and also disqualify. Perhaps within the panel function
    # TODO: this may run too long-- set a cycle-limit exceeding which would also 
    #       end the attempt to estimate the radius for these phis
    # DONE: implement logging at the phi level-- don't retain model results or 
    #       stats about the current simulated dataset, just the global part of 
    #       ptpnl_summary especially the phi
    # TODO: implement error logging for ptpnl_ functions
    # DONE: in normal result logging for ptpnl_ functions retain the verdict as well
    # DONE: consider not bothering with na.omit since glm probably does it
    #       internally anyway. Makes it easier to find fraction NA too (above TODO)
    cycle <- cycle+1; tfoffset <- length(list_tfresp);
    if(file.exists(savetrigger)) {
      print('  Saving. ');
      save(logenv,file=savefile);
      file.remove(savetrigger);
    }
    if(file.exists(debugtrigger)){
      browser();
    }
    if(file.exists(sourcepatch)) {
      print('  Patching  ');
      source(sourcepatch,local = T);
      file.rename(sourcepatch,paste0(sourcepatch,'.bak'));
    }
  };
  pnlph(ppdat,preds['radest',],logenv=logenv,time=Sys.time()-t0);
  if(lims['status']==1){
    cat('Success: ');
    #if(first_success){first_success<-F; browser();}
    #for(pp in pnfit[preds['conv',]==1]) {
    # it's not enough to test for convergence-- we have to also make sure the 
    # converged results fall within the limits. Not and-ing the below with the
    # convergence criteria because IIRC for us to even get this far, they both
    # have to have converged
    for(pp in pnfit[lims[paste0('notfailed.',pnfit)]==1]) {
      ppcoords <- backtrans(pol2crt(c(preds['radest',pp],phi)));
      ppdat <-ptsim(ppcoords,refcoords=btrefcoords,...); 
      # TODO: not sure we actually needs to explicitly specify index, they seem
      # to be properly handling it anyway
      pnlst[[pp]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,pp));
      # for each pp (verdict-returning panel function) we generate a separate dataset therefore
      # we need to iterate over all the summary-only non-verdict functions for each of these datasets
      for(qq in pninfo) pnlst[[qq]](ppdat,preds['radest',pp],logenv=logenv,index=c('coords',philabel,qq));
    };
  } else cat('Failure: ');    #if(first_fail){first_fail<-F; browser();}}
  logenv$allpoints[[philabel]] <- data.table(rad=testrd,rbind(phi),testtf);
  if(isTRUE(logenv$state$powertrip$console_log)) cat('radii= ',try(preds['radest',]),'\tphi= ',try(phi),'\tlims= ',c(maxrad,lims[-3]),'\n');
  # DONE: add back in the dynamic script execution and the external exit directive
  # DONE: add a phi -> radius prediction step (multivariate, whole parameter space)
  # DONE: for prioritizing phis, also model the runtime to get most uncertainty
  #       per second of runtime or per simulation
  # TODO: also, model the 'dead-zones' -- places where we had to give up --
  #       and exclude them
  # DONE?: figure out why negative radii are being allowed and fix
  # DONE?: figure out why plots look so wierd... is pol2crt wrong?
  # DONE: finalize outer function
  # DONE?: launch the linear model version
  # TODO: try to resurrect simsurve and survwrapper
}

#' Sobering thought: right now this averages about 0.006 seconds per simulation
#' opsurv() can do about 13 simulations per second
# > system.time(.junk<-replicate(200,opsurv(foo,foo,model='lm')))
# user  system elapsed 
# 15.477   0.008  15.734 
# > 200/15.5
# [1] 12.90323
#' ...which means that it tacks on 0.0775 or ~0.08 seconds per simulation 
#' meaning that is multiplies runtimes 14.3-fold!

powertrip<-function(logenv=logenv,refcoords
                    #,maxs=c(2,4.5,6),mins=c(-3.1,-1.3,-6)
                    ,maxs=c(20,20,20),mins=c(-20,-20,-20)
                    ,npoints=50,nphis=2,nrads=20,numse=2
                    # which fraction of the most impactful phis should we
                    # model each time?
                    ,wd=paste0(getwd(),'/'),savetrigger=paste0(wd,'pt_savedata')
                    ,debugtrigger=paste0(wd,'pt_debug')
                    ,savefile=paste0(wd,'pt_result.rdata')
                    ,sourcepatch=paste0(wd,'pt_sourcepatch.R')
                    ,topsourcepatch=paste0(wd,'pt_top_sourcepatch.R')
                    ,bestfrac = 0.5
                    # these run at the level of each radius
                    ,pnlst=list(#lm=ptpnl_lm
                                   lm2=update(ptpnl_lm,fname='lm2',frm=yy~(.)*group)
                                   ,tt=ptpnl_tt
                                   ,wx=ptpnl_wx
                                   ,simsm=ptpnl_simsumm)
                    # this runs on each set of phis
                    ,pnlph=ptpnl_phisumm
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
  # maybe this is a good place to initialize logenv for use everywhere else in powertrip
  # so we don't have to keep deriving names in multiple places or passing around too many arguments
  # TODO: consider this being a function with access to its calling environment?
  # names of the functions in pnlist that each return a TRUE/FALSE verdict
  logenv$names$pnfit <- pnfit <- names(pneval_)[pneval_];
  # names of the other functions in pnlist that produce summary statistics for each simulated population
  logenv$names$pninfo <- pninfo <- names(pneval_)[!pneval_];
  # names of the columns that hold the predicted lengths of radii for each pnfit function, for pt2df
  logenv$names$radnames <- radnames <- paste0('rad.',logenv$names$pnfit);
  # names that make_phis will use for radius lengths and their SEs, respectively
  logenv$names$fnames<-paste0(logenv$names$radnames,'.fit');
  logenv$names$snames<-paste0(logenv$names$radnames,'.se');
  # names of the phi columns, for pt2df
  logenv$names$phinames <- phinames <- paste0('phi',seq_len(nphis));
  # names of the notfailed elements in the lims vector, for subsetting which rads to report
  logenv$names$notfailed <- paste0('notfailed.',logenv$names$pnfit);
  # alist to be used in the fields argument of pt2df as invoked within env_fitinit()
  logenv$names$fields <- alist(setNames(ifelse(lims[ptenv$names$notfailed]==T,preds['radest',],NA),ptenv$names$radnames),setNames(phi,ptenv$names$phinames),nsims=nsims);

  # The below is a catch-all way to do runtime patches to most arguments
  logenv$state$powertrip <- environment();
  console_log <- T;
  # Increment phicycle from previous run by 1 if it exists, otherwise start with 1
  phicycle <- c(logenv$state$phicycle,0)[1]+1;
  # export final phicycle when exiting
  on.exit(logenv$state$phicycle<-phicycle);

  while(phicycle < maxphicycle){
    if(file.exists(topsourcepatch)) {
      print('  Patching Top Level ');
      source(topsourcepatch,local = T);
      file.rename(topsourcepatch,paste0(topsourcepatch,'.bak'));
    }
    actualpoints <- 0; phis <- c();
    while(actualpoints < 30){
      phis <- rbind(phis,make_phis(logenv=logenv,npoints = npoints,maxs=maxs,mins=mins
                      ,nphis = nphis,numse = numse));
      phis <- subset(phis,maxrad>0);
      if(!is(phis,'data.frame')||nrow(phis)<1||nrow(subset(phis,maxs<0|mins<0|mins>=maxs|maxs>maxrad))>0) {print('Jacked phis created!');browser();}
      #logenv$temp$maxrads <- maxrads <- apply(phis,1,pollim,maxs=maxs,mins=mins);
      actualpoints <- nrow(phis);
    }
    # DEBUG BREAK
    #debugonce(make_phis);
    for(ii in seq_len(actualpoints)){
      if(console_log) cat(phicycle,'\t',ii,'\t');
      phi_radius(phi=unlist(phis[ii,phinames]),maxrad=phis[ii,'maxrad'],pnlst=pnlst
                 ,refcoords = refcoords
                 ,pnlph=ptpnl_phisumm
                 ,pneval=pneval_,pnfit=pnfit,pninfo=pninfo
                 ,logenv=logenv,max=phis[ii,'maxs'],min=phis[ii,'mins'],nrads=nrads
                 ,numse=numse
                 ,ptsim=ptsim
                 ,wd=wd
                 ,savetrigger=savetrigger
                 ,debugtrigger=debugtrigger
                 ,savefile=savefile
                 ,sourcepatch=sourcepatch
                 ,phicycle=phicycle,...);
    }
    #if(phicycle==1) browser();
    phicycle<-phicycle+1
  }
  data.frame(t(sapply(logenv$coords,function(xx) with(xx$summ[[1]],c(
    phi=phi,radii=ifelse(preds['conv',]==1,preds['radest',],NA)
    ,time=time,nsims=nsims,cycle=cycle)))));
  #browser();
};
# Timing stopped at: 6573 5411 6740 ... about 1.87 hours to try 1000 pairs of phis
