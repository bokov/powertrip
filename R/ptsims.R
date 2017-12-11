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
  out <- data.frame(group=rep(c('control','treated'),each=nn)
                    ,xs
                    ,yy=ifelse(seq_len(2*nn)<=nn,refcoords[1]+xs %*% refcoords[-1],coords[1]+xs %*%coords[-1])+rnorm(2*nn));
}

#' Needs: rlgst.R, rlogmake.R, and the eha package
ptsim_surv <- function(coords,nn=100,refcoords=c(2.433083e-05, 0.005, 3e-11, 0.0015),type='gm',...){
  coords<- coords+refcoords;
  out <- data.frame(group=rep(c('control','treated'),each=nn)
                    ,yy=c(simsurv(nn,type,refcoords)
                          ,simsurv(nn,type,coords))
                    ,cc=1);
}
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
ptpnl_passthru <- new.ptpnl("passthru"
                            , fit = data[[1]][1]
                            , result = list(summary = fit, coords = coords, call = match.call(expand.dots = T))
                            , eval = fit[1] == 1);
#' panel function to be used just for generating a summary result

ptpnl_qntile <- new.ptpnl("qntile"
                          , fit = quantile(data[, 1], c(0.05, 0.1, 0.5, 0.9, 0.95))
                          , result = list(summary = fit, coords = coords, call = match.call(expand.dots = T)));

# TODO: ptpnl_popsummary, ptpnl_phisummary

ptpnl_simsumm <- new.ptpnl('simsm'
                              ,fit = split(data.frame(data),data.frame(data)[,1])
                              ,result = c(summaries=sapply(fit,function(xx) sapply(xx,summary,simplify=F),simplify=F)));
                              
ptpnl_phisumm <- new.ptpnl('summ'
                           ,fit=c()
                           #,fit = split(data.frame(data),data.frame(data)[,1])
                           ,time=quote(Sys.time())
                           ,literals=c('coords','cycle','phi','preds','lims','maxrad','phicycle')
                           ,result = c(sapply(intersect(literals,names(callingframe)),function(xx) callingframe[[xx]],simplify=F)
                                       ,nsims=length(callingframe$list_tfresp)
                                       ,time=time)
);
#' NOTE: avoid creating variables that match the regexp "^phi[0-9]$" because 
#' env_fitinit() will mistake them for phinames
#' 
ptpnl_summary <- new.ptpnl('summ.old'
                           ,fit = split(data.frame(data),data.frame(data)[,1])
                           #,cycles.=quote(cycle)
                           #,nsims.=quote(length(list_tfresp))
                           #,phi.=quote(phi)
                           ,time=quote(Sys.time())
                           ,literals=c('coords','cycle','phi','preds','lims','maxrad','phicycle')
                           #,philabel_= quote(callingframe$philabel)
                           #,preds.=quote(preds)
                           #,philabel_=quote(philabel)
                           ,result = c(summaries=sapply(fit,function(xx) sapply(xx,summary,simplify=F),simplify=F)
                                       ,sapply(intersect(literals,names(callingframe)),function(xx) callingframe[[xx]],simplify=F)
                                       ,nsims=length(callingframe$list_tfresp)
                                       ,time=time)
                                          # ,preds=preds,phi=phi,time=time,cycles=cycles,nsims=nsims)
                           #,index=substitute(c('coords',`philabel_`,'summ')));
);

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
                      ,eval= summary(fit)$table[matchterm,'p']<psig
                      # should also be able to handle Surv(yy,cc)~.
                      ,frm=Surv(yy)~.,psig=0.05
                      ,matchterm='grouptreated');
#' coxph
ptpnl_cx <- new.ptpnl('cx'
                      ,fit=coxph(formula=frm,data)
                      ,result = list(model=broom::glance(fit),detect=eval(eval.))
                      ,eval= summary(fit)$coef[matchterm,'Pr(>|z|)']<psig
                      # should also be able to handle Surv(yy,cc)~.
                      ,frm=Surv(yy)~group,psig=0.05
                      ,matchterm='grouptreated');

#' coxph
#' The following work (after running the ptsim_nlin example near top of script):
#' 
#' Group alone:
#ptpnl_lm(foo,c(1.06,5.124,2,-0.5),matchterm = c('grouptreated','grouptreated:X1','grouptreated:X2','grouptreated:X3'));
#' With interactions:
#ptpnl_lm(foo,c(1.06,5.124,2,-0.5),frm=yy~(.)*group, matchterm = c('grouptreated','grouptreated:X1','grouptreated:X2','grouptreated:X3'));

#' there can be a list of these, and they can repeat
#baz<-c(ptpnl_passthru,ptpnl_qntile,ptpnl_passthru);
#' you can generate names for them automatically
#names(baz) <- make.unique(sapply(baz,attr,'fname'));
#' you can generate a logical vector indicating which of them should return a TRUE/FALSE decision
#sapply(baz,attr,'eval');
#' you can see how each of them were created
#sapply(baz,attr,'call');
#' Of course you can evaluate them all on the same dataset, the main purpose of them
#logenv <- new.env();
#for(ii in names(baz)) baz[[ii]](mtcars,c(-1,1),logenv=logenv,index=ii);

#' getCall method for a ptpnl_fn function.
getCall.ptpnl_fn <- function(xx,...) attr(xx,'call');
#' This enables the following to work:
#ptpnl_quartl <- eval(update(ptpnl_qntile,fname='quartl',fit=quantile(data[,1],c(0,.25,.5,.75,1))));
#' Panels
# ptpnl_passthru <- function (data, coords, logenv = NULL, errenv = NULL, index, 
#           pninfo = F, ...) 
# {
#   if (pninfo) {
#     return(list(fname = "passthru", eval = TRUE))
#   }
#   fit <- try(data[1])
#   if (is(fit, "try-error")) {
#     if (!is.null(errenv)) 
#       errenv[[index]] <- fit
#     return(NA)
#   }
#   else {
#     if (!is.null(logenv)) 
#       logenv[[index]] <- list(summary = fit[1], coords = coords, 
#                               call = sys.call(sys.parent()))
#   }
#   return(fit[1])
# }
# 
#' 
ptpnl_passthru.bak <- function(data,coords=NULL,...){
  list(outcome=data
       ,coords=coords
       ,call=sys.call(sys.parent()));
}

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
