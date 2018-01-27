#' Here we will have the various simulation modules and 
#' evaluation panels.
#' 
#' 
#' Simulations
ptsim_binom <- function(coords,...){
  1-rbinom(1,1,mvtnorm::dmvnorm(coords)/mvtnorm::dmvnorm(rep_len(0,length(coords))));
}

#' Title ptsim_2lin
#'
#' Garden variety linear model, with a categorical and a numeric predictor
#' @param coords    numeric vector of coordinates on cartesian scale
#' @param nn        sample size
#' @param refcoords coordinates of control group
#' @param ... 
#'
#' @return A data frame
#' @export
#'
#' @examples simdata <- ptsim_2lin(c(4,-3),40);
#' 
#' This function works with coords vector > 2, silently ignoring the extras
#' and the test_harness() function that calls this works also
ptsim_2lin <- function(coords,nn=100,refcoords=rep_len(0,length(coords)),...){
  coords <- coords + refcoords;
  ooc <- data.frame(group='control',xnum=rnorm(nn));
  ooc$yy <- with(ooc,refcoords[1]+refcoords[2]*xnum+rnorm(nn));
  oot <- data.frame(group='treated',xnum=rnorm(nn));
  oot$yy <- with(oot,coords[1]+coords[2]*xnum+rnorm(nn));
  return(rbind(ooc,oot));
}

ptsim_nlin <- function(coords,nn=100,lcoords=length(coords),lvars=lcoords-1,refcoords=coords*0,...){
  # notice that refcoords automagically are the right length regardless of whether or not the above if()
  # statement is triggered-- I suspect this is due to lazy evaluation-- refcoords remain a call until 
  # they are used in the following expression...
  coords <- coords + refcoords;
  xs <- matrix(rnorm(nn*lvars*2),nrow=nn*2,ncol=lvars);
  data.frame(group=rep(c('control','treated'),each=nn),xs
                    ,yy=ifelse(seq_len(2*nn)<=nn,refcoords[1]+xs %*% refcoords[-1],coords[1]+xs %*%coords[-1])+rnorm(2*nn));
}

#' ### Sketch for future ptsim_lm function, or more generally, location-scale so LS...
#' 
#' The first third of the coordinates are the effect sizes of treatment group membership, 
#' the second third are offsets to the *predictors* (i.e. biased samples, set to 0 if not
#' simulating bias), and the final third are offsets to standard deviation (i.e. 
#' heteroscedasticity, set to 0 if not desired). If one coordinate is left over it's used
#' for sample size, otherwise nn argument is used. If two are left over, it's an error.
#' The dists should be a list of rnorm()-like functions that take sample size as first
#' parameter, location as second, and scale as third
ptsim_ls <- function(coords,nn=100,nc=length(coords)%/%3,refcoords=rep_len(1,nc*3)
                           ,refdists=list(rnorm)[rep_len(1,nc)],dists=refdists){
  cseq <- seq_len(ntot<-3*nc);
  if((diffnc <- length(coords) - ntot)==2) {
    stop('The length of the coords argument should be either divisible by 3 or one extra if simulating sample sizes');
  } else if(diffnc==1) {nn <- coords[ntot+1]; coords <- coords[cseq];}
  # we will split the coordis using the csplt vector with rc the refcoords and tc as the coords
  csplt <- (cseq-1)%/%3; rc <- split(refcoords,csplt); tc <- split(refcoords + coords,csplt);
  ooref <- cbind(mapply(function(fn,lc,sc){fn(nn,lc,sc)},refdists,rc[[2]],rc[[3]],...),group='control');
  ootrt <- cbind(mapply(function(fn,lc,sc){fn(nn,lc,sc)},dists,tc[[2]],tc[[3]],...),group='treated');
  ooref$yy <- ooref[,seq_len(nc)]%*%rc[[1]];
  ootrt$yy <- ootrt[,seq_len(nc)]%*%tc[[1]];
  # oops... but then we need to also have a random error term rather than the errors-in-variables implied
  # by above... also, coords not necessarily monotonically increasing... oh well, like I said, just a sketch
  # so I don't forget, not an actual tested, working prototype.
}

#' *Note*: The next few functions which use simsurv() require coords to have 
#' parameters in the following order: 
#' 
#' 'IMR' (or 'a'), 'RoA' (or 'b'), 'EH' (or 'c'), if applicable 'LL' (or 's')
#' and if applicable sample size last

#' Needs: rlgst.R, rlogmake.R, and the eha package
ptsim_surv <- function(coords,nn=100,refcoords=c(2.433083e-05, 0.005, 3e-11, 0.0015),type='lm',...){
  # 
  coords<- coords*refcoords;
  out <- try(data.frame(group=rep(c('control','treated'),each=nn)
                    ,yy=c(simsurv(nn,type,refcoords)
                          ,simsurv(nn,type,coords))
                    ,cc=1));
  if(is(out,'try-error')) return(expand.grid(group=c('control','treated'),yy=-1,cc=-1)) else out;
}

#' This is like ptsim_surv except it uses the last coord as the
#' sample size (the switch statement in first line of body shows exactly how many
#' params each model takes). The refcoords should be the same length as coords
#' and appropriate to the model and the defaults are just a set of values that 
#' happen to work with Logistic-Makeham and are for demonstration purposes only. 
#' Not going to slow down something that runs so often by putting in more input 
#' checking here.
ptsim_srvn <- function(coords,refcoords=c(2.433083e-05, 0.005, 3e-11, 0.0015,1)
                     ,type=c('e','g','gm','lm'),...){
  lc <- switch(type<-match.arg(type),e=1,g=2,gm=3,lm=4); 
  out <- NULL;
  # the coords start out as offsets from the refcoords, and here they are 
  # turned into actual values. The [lc+1] value is the sample size and is the same
  # between the control/reference group and the treatment group
  coords <- coords*refcoords; refcoords[lc+1] <- coords[lc+1];
  if(round(coords[lc+1])>=round(refcoords[lc+1])){
    # coords[lc+1] is shared as the sample size for both groups and apparently the 
    # proper co... what was I going to write here? Don't know. Darn.
    # TODO: replace cc with 
    # cc = sample(1:0,coords[lc+1],rep=T,prob=c(coords[lc+2],1-coords[lc+2]))
    # TODO: add recruitment age (lc+2)
    out <- try(data.frame(group=rep(c('control','treated'),each=coords[lc+1])
                          ,yy=c(simsurv(coords[lc+1],type,refcoords[1:lc])
                                ,simsurv(coords[lc+1],type,coords[1:lc]))
                          ,cc=1));
  } else cat(' sim-err: N below minimum '); 
  if(is.data.frame(out)) return(out) else {
    if(is(out,'try-error')) cat(out[1]);
    return(expand.grid(group=c('control','treated'),yy=-1,cc=-1));
  }
}

simsurvagenroll <- function(params,nn,agenroll=0,lc=length(params),type){
  if(missing(type)) type <- switch(lc,'e','g','gm','lm');
  shp <- if(agenroll>0) {
    do.call(srvshp,c(x=agenroll,setNames(as.list(params),c('a','b','c','s')[1:lc]),model=type));
  } else 1;
  out <- simsurv(nn/shp,type,params)-agenroll;
  while((ll<-length(out<-out[out>0]))<as.integer(nn)){
    out<-c(out,simsurv((nn-ll)/shp,type,params));
    };
  sample(out,nn);
}
#' ## example... try recruiting at every month of age, 0 to 960
# coords <- c(3.02495622562167e-06, 0.00766970877053115, 1.97042457165941e-05,30, 0, 20)'
# lc<-3;
# foo<-sapply(0:(80*12),function(ii) simsurvagenroll(coords[1:lc],1000,ii));
# plot(survfit(Surv(foo[,1])~1),conf.int=F,las=1,yscale=100,bty='n',xscale=12
#      ,xlab='Survival After Recruitment (years)',ylab='% Alive'
#      ,main='Survival Curves at Various Recruitment Ages\n(red = 0, blue = 80)')
# for(ii in 2:961) lines(survfit(Surv(foo[,ii])~1)
#                        ,col=rainbow(961,start=0,end=4/6,alpha=0.2)[ii],conf.int=F); 
#' or year of age
# for(ii in seq(12,960,by=12)) lines(survfit(Surv(foo[,ii])~1)
#                                    ,col=rainbow(961,start=0,end=4/6,alpha=0.7)[ii],conf.int=F)
#' TODO: Replicate a bunch, run opsurv on each one, and discover relationship between hazard 
#' param residuals and the age at recruitment

#simsurv (n, type = "g", p = c(2.433083e-05, 0.005, 3e-11, 0.0015)) 
#' The following works!
#' 
#ptsim_nlin(pol2crt(c(1.06,5.124,2,-0.5))) -> foo;
  
#' TODO: in addition to a logenv object, return a character vector specifying a
#' path through that logenv
#' TODO: optional formula argument
#' TODO: look in dots for additional arguments (with )
new.ptpnl <- function(fname,fit,result,eval.,...){
  # we use substitute on the fit and result arguments to turn them into 
  # unevaluated calls from whatever fragile state they are before they are 
  # first evaluated
  fit <- substitute(fit); result <- substitute(result); 
  # we also need to do that to the eval argument but that one is optional, if 
  # missing it means the user wants this function to return summary results but
  # no T/F evaluation so we first check if it's missing
  if(!missing(eval.)) eval. <- substitute(eval.);
  dots <- list(...);
  #if(is.null(frm)) warning('If your model fit function uses a formula, please consider setting that argument to "frm" and then setting the "frm" argument of your new.ptpnl() call to the actual formula. This will allow you to update just the formula without rebuilding the whole function.');
  # below is our template for the function-- the commented parts will get
  # inserted or altered
  # TODO: get rid of the pninfo feature, redundant with attributes
  oo <- function(data,coords,logenv=NULL,errenv=NULL,index=c(paste(c('x',coords),collapse ='_'),fname),...){
    if(F){
      return(NA); # metadata about function goes here -- fname and eval T/F
    }
    NA # fit goes here
    if(is(fit,'try-error')){
      # the first element of index, which is a character vector, is the object
      # in the environment errenv; the rest are a path to the place where the
      # output will be stored. It should always work to use indexes of length
      # 2, otherwise may error if specified nodes do not already exist
      if(!is.null(errenv)) errenv[[index[1]]][[index[-1]]] <- fit;
      # return(NA) goes here if !missing(eval)
    } else {
      NA  # summary result goes here if(!is.na(logenv)) logenv[[index]]
    }
  };
  # we make the default value for the index variable be that of fname
  formals(oo)$callingframe <- substitute(parent.frame());
  formals(oo)$philabel <- substitute(c(callingframe$philabel,'unknown')[1]);
  formals(oo)$index <- substitute(c('coords',philabel,fname));
  for(ii in names(dots)) formals(oo)[[ii]] <- dots[[ii]];
  if(!'frm' %in% names(formals(oo))&&!missing(eval.)) warning(
   "  You did not specify the optional argument 'frm' when calling the new.ptpnl()
    function. This is allowed, but if your fit argument relies on a formula it's
    better to use frm instead of an explicit formula ( e.g. lm(frm,data) ) and 
    specify the formula in the frm argument to new.ptpnl(). That way, you will 
    be able to swap out what formula your fit function uses without having to 
    build a separate function just for that."
  );
  # below is the metadata our function will return if its pninfo argument is 
  # set to TRUE
  # body(oo)[[2]][[3]][[2]][[2]] <- list(fname=fname,eval=!missing(eval),call=match.call());
  # Below is where the fit for this function gets created. This is done in two
  # steps so that only the second occurrence of fit is substituted
  body(oo)[[3]] <- substitute(FIT <- try(fit));
  body(oo)[[3]][[2]] <- quote(fit);
  # If the function is supposed to return a TRUE/FALSE value, return NA when 
  # there is an error, otherwise, if it isn't supposed to return anything then
  # don't do anything 
  body(oo)[[4]][[3]][[3]] <- if(missing(eval.)) quote(return(NULL)) else quote(return(NA));
  # Here we make the function return the results of its eval statement to a named
  # entry in the logenv environment... but only if it exists
  #body(oo)[[4]][[4]][[2]] <- substitute(if(!is.null(logenv)) logenv[[index[1]]][[index[-1]]] <- result);
  body(oo)[[4]][[4]][[2]] <- substitute(if(!is.null(logenv)) {val.<-eval(result);deepassign(obj = logenv,path = index,val=val.)});
  # we return our TRUE/FALSE result, if this is a function that returns a value
  if(missing(eval.)) body(oo)[[5]] <- quote(return(NULL)) else {
    body(oo)[[5]] <- substitute(return(eval.));
  }
  attr(oo,'call') <- match.call();
  attr(oo,'eval') <- !missing(eval.);
  attr(oo,'fname') <- fname;
  class(oo) <- c('ptpnl_fn',class(oo));
  oo;
}
#' Example use:
#' 
#' panel function that returns a decision
#ptpnl_passthru <- new.ptpnl("passthru"
#                            , fit = data[[1]][1]
#                            , result = list(summary = fit, coords = coords, call = match.call(expand.dots = T))
#                            , eval = fit[1] == 1);
#' panel function to be used just for generating a summary result

ptpnl_qntile <- new.ptpnl("qntile"
                          , fit = quantile(data[, 1], c(0.05, 0.1, 0.5, 0.9, 0.95))
                          , result = list(summary = fit, coords = coords, call = match.call(expand.dots = T)));

# DONE: ptpnl_popsummary, ptpnl_phisummary

ptpnl_simsumm <- new.ptpnl('simsm'
                              ,fit = split(data.frame(data),data.frame(data)[,1])
                              ,result = c(summaries=sapply(fit,function(xx) sapply(xx,summary,simplify=F),simplify=F)));
# TODO: write a version of ptpnl_simsumm that also does integrals under survival curve at certain timepoints (e.g. 480 and 1080 months)
                              
ptpnl_phisumm <- new.ptpnl('summ'
                           ,fit=c()
                           #,fit = split(data.frame(data),data.frame(data)[,1])
                           ,time=quote(Sys.time())
                           # only put something in literals if its small enough that you want to keep it in its 
                           # entirety
                           ,literals=c('cycle','phi','preds','lims','maxrad','phicycle')
                           ,result = c(sapply(intersect(literals,names(callingframe)),function(xx) callingframe[[xx]],simplify=F)
                                       ,list(nsims=length(callingframe$list_tfresp)
                                       # cbinding to NA so there would consistently always be an NA row,
                                       # so remember to subtract 1 from this if using
                                       ,cyclecoords=apply(rbind(callingframe$cyclecoords,NA),2,summary)
                                       ,cycleradii=summary(callingframe$testrd)
                                       ,lastsim=summary(rbind(callingframe$iidat,NA))
                                       # subtract 1 from the first and last fields
                                       ,testtf=table(interaction(rbind(callingframe$testtf,T,F)))
                                       ,ts=Sys.time()
                                       ,time=time))
);
#' NOTE: avoid creating variables that match the regexp "^phi[0-9]$" because 
#' env_fitinit() will mistake them for phinames
#' 
# ptpnl_summary <- new.ptpnl('summ.old'
#                            ,fit = split(data.frame(data),data.frame(data)[,1])
#                            #,cycles.=quote(cycle)
#                            #,nsims.=quote(length(list_tfresp))
#                            #,phi.=quote(phi)
#                            ,time=quote(Sys.time())
#                            ,literals=c('coords','cycle','phi','preds','lims','maxrad','phicycle')
#                            #,philabel_= quote(callingframe$philabel)
#                            #,preds.=quote(preds)
#                            #,philabel_=quote(philabel)
#                            ,result = c(summaries=sapply(fit,function(xx) sapply(xx,summary,simplify=F),simplify=F)
#                                        ,sapply(intersect(literals,names(callingframe)),function(xx) callingframe[[xx]],simplify=F)
#                                        ,nsims=length(callingframe$list_tfresp)
#                                        ,time=time
#                                        ,tstamp=Sys.time())
#                                           # ,preds=preds,phi=phi,time=time,cycles=cycles,nsims=nsims)
#                            #,index=substitute(c('coords',`philabel_`,'summ')));
# );

# TODO: change all the below functions to automatically tack on a detect=eval(eval.)
# onto their result argument
ptpnl_wx <- new.ptpnl("wx"
                      , fit = wilcox.test(frm, data)
                      , result = list(model=broom::glance(fit),detect=eval(eval.))
                      , eval. = fit$p.value < psig
                      , frm = yy ~ group, psig = 0.05);

ptpnl_tt <- new.ptpnl("tt"
                      , fit = t.test(frm, data)
                      , result = list(model=broom::glance(fit),detect=eval(eval.))
                      , eval. = fit$p.value < psig
                      , frm = yy ~ group, psig = 0.05);

ptpnl_lm <- new.ptpnl("lm"
                      , fit = lm(formula = frm, data)
                      , result = list(model=broom::glance(fit),detect=eval(eval.))
                      , eval. = p.adjust(with(summary(fit), coefficients[, "Pr(>|t|)"][rownames(coefficients) %in% 
                                                                                                matchterm])) < psig
                      , frm = yy ~ ., psig = 0.05, matchterm = substitute(paste0('grouptreated',c('',paste0(':',names(data))))));

#' survreg, i.e. Weibull accelerated failure time
ptpnl_sr <- new.ptpnl('sr'
                      ,fit=survreg(formula=frm,data)
                      ,result = list(model=broom::glance(fit),detect=eval(eval.))
                      ,eval.= summary(fit)$table[matchterm,'p']<psig
                      # should also be able to handle Surv(yy,cc)~.
                      ,frm=Surv(yy)~.,psig=0.05
                      ,matchterm='grouptreated');
#' Detect differences between quantiles. A measure of *practical*, as opposed to
#' statistical, significance. Defaults to 12 months difference in either direction
#' for the third quantile, but can be set to different values with optional arguments
#' `which` (which part of summary) or `cutoff`
ptpnl_diff <- new.ptpnl(fname = "diff"
                        ,fit = sapply(split(data.frame(data)$yy,data.frame(data)[, 1])
                                      , summary, simplify = F)
                        ,result = c(detect = eval(eval.), unlist(fit))
                        ,eval. = abs(fit[[2]] - fit[[1]])[which] > cutoff
                        ,which = "3rd Qu.", cutoff = 12);
#' coxph
ptpnl_cx <- new.ptpnl('cx'
                      ,fit=coxph(formula=frm,data)
                      ,result = list(model=broom::glance(fit),detect=eval(eval.))
                      ,eval.= summary(fit)$coef[matchterm,'Pr(>|z|)']<psig
                      # should also be able to handle Surv(yy,cc)~.
                      ,frm=Surv(yy)~group,psig=0.05
                      ,matchterm='grouptreated');

ptpnl_gm <- new.ptpnl('gm'
                      ,fit = list(h0=opsurv(x=data$yy[data$group==matchterm],y=data$yy[data$group!=matchterm],model='gm',tlog=T,par=c(3e-6,8e-3,2e-5,0),cons=c(0,0,0,0))
                                ,h1=opsurv(x=data$yy[data$group==matchterm],y=data$yy[data$group!=matchterm],model='gm',tlog=T,par=c(3e-6,8e-3,2e-5,0)))
                      ,result = list(model=fit, detect=eval(eval.))
                      ,eval.= pchisq(2 * (fit$h1$maximum - fit$h0$maximum), df = 3,lower.tail = F) < psig
                      ,frm=Surv(yy)~group,psig=0.05
                      ,matchterm='control');

#' getCall method for a ptpnl_fn function.
getCall.ptpnl_fn <- function(xx,...) attr(xx,'call');

ptpnl_2lin <- function(data,coords=NULL,termrxp='group',sig=0.05,...){
  fits <- list(add=try(lm(yy~.,data))
              ,intx=try(lm(yy~(.)^5,data))
              ,addaov=try(aov(yy~.,data))
              ,intxaov=try(aov(yy~(.)^5,data)));
  outcomes<-sapply(fits
                   ,function(xx) tryCatch(any(p.adjust(
                     subset(broom::tidy(xx),grepl(termrxp,term))$p.value
                     )<sig)
                     ,error=function(ee) NA)
                   );
  results<-lapply(fits,function(xx) if(class(xx)[1]=='try-error') xx else {
    try(broom::glance(xx));
  });
  list(outcome=outcomes
       ,coords=coords
       ,results=results
       ,call=sys.call(sys.parent())
       );
}

deepassign <- function(obj,path,val,append=T){
  # assuming that obj is list-like and path is a path through nested named branches and val gets assigned at the terminal branch
  # in other words, given a path into a hierarchical list-like object, and a value
  # when this function is done running that node will exist and the value will be
  # assigned to it. This has been tested on environments.
  # The optional append argument means that by default, the val will be appended 
  # as an anonymous sub-object to the terminal node of path preserving any previous
  # value/s. If append is set to F, then val will create or replace the terminal
  # node of the path. Unexpected things can happen if you create nodes with append=F
  # and later populate them with append=T because you'll be appending to the existing
  # value of the node rather than appending a new sub-node. Can't think of a clean
  # way to detect and warn of this right now, so this comment is your warning!
  if(length(path)<=1) parsepath='obj' else {
    if(!path[1] %in% names(obj)) obj[[path[1]]] <- list();
    parsepath <- paste0("obj[['",path[1],"']]");
    for(ii in seq_along(path[-length(path)])[-1]) {
      newparsepath <- paste0(parsepath,"[['",path[ii],"']]");
      if(!path[ii] %in% names(eval(parse(text=parsepath)))) {
        eval(parse(text=paste0(newparsepath,'<- list()')));
      }
      parsepath <- newparsepath;
    }
  }
  # at this point we have an obj that should be capable of taking an 
  # assignment to the final name of path regardless
  finalpath <- paste0(parsepath,"[['",tail(path,1),"']]");
  if(append){
    lorig <- eval(parse(text=paste0('length(',finalpath,')')));
    eval(parse(text=paste0(finalpath,"[[",lorig+1,"]]<-val")));
  } else eval(parse(text=paste0(finalpath,"<-val")));
}
