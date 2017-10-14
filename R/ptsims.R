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
ptsim_2lin <- function(coords,nn=100,refcoords=rep_len(0,length(coords)),...){
  coords <- coords + refcoords;
  ooc <- data.frame(group='control',xnum=rnorm(nn));
  ooc$yy <- with(ooc,refcoords[1]+refcoords[2]*xnum+rnorm(nn));
  oot <- data.frame(group='treated',xnum=rnorm(nn));
  oot$yy <- with(oot,coords[1]+coords[2]*xnum+rnorm(nn));
  return(rbind(ooc,oot));
}

new.ptpnl <- function(fname,fit,result,eval,...){
  # we use substitute on the fit and result arguments to turn them into 
  # unevaluated calls from whatever fragile state they are before they are 
  # first evaluated
  fit <- substitute(fit); result <- substitute(result); 
  # we also need to do that to the eval argument but that one is optional, if 
  # missing it means the user wants this function to return summary results but
  # no T/F evaluation so we first check if it's missing
  if(!missing(eval)) eval <- substitute(eval);
  # below is our template for the function-- the commented parts will get
  # inserted or altered
  # TODO: get rid of the pninfo feature, redundant with attributes
  oo <- function(data,coords,logenv=NULL,errenv=NULL,index,pninfo=F,...){
    if(pninfo){
      return(NA); # metadata about function goes here -- fname and eval T/F
    }
    NA # fit goes here
    if(is(fit,'try-error')){
      if(!is.null(errenv)) errenv[[index]] <- fit;
      # return(NA) goes here if !missing(eval)
    } else {
      NA  # summary result goes here if(!is.na(logenv)) logenv[[index]]
    }
  };
  # we make the default value for the index variable be that of fname
  formals(oo)$index <- fname;
  # below is the metadata our function will return if its pninfo argument is 
  # set to TRUE
  body(oo)[[2]][[3]][[2]][[2]] <- list(fname=fname,eval=!missing(eval),call=match.call());
  # Below is where the fit for this function gets created. This is done in two
  # steps so that only the second occurrence of fit is substituted
  body(oo)[[3]] <- substitute(FIT <- try(fit));
  body(oo)[[3]][[2]] <- quote(fit);
  # If the function is supposed to return a TRUE/FALSE value, return NA when 
  # there is an error, otherwise, if it isn't supposed to return anything then
  # don't do anything 
  body(oo)[[4]][[3]][[3]] <- if(missing(eval)) quote(return(NULL)) else quote(return(NA));
  # Here we make the function return the results of its eval statement to a named
  # entry in the logenv environment... but only if it exists
  body(oo)[[4]][[4]][[2]] <- substitute(if(!is.null(logenv)) logenv[[index]] <- result);
  # we return our TRUE/FALSE result, if this is a function that returns a value
  if(missing(eval)) body(oo)[[5]] <- quote(return(NULL)) else {
    body(oo)[[5]] <- substitute(return(eval));
  }
  attr(oo,'call') <- match.call();
  attr(oo,'eval') <- !missing(eval);
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
#' there can be a list of these, and they can repeat
#baz<-c(ptpnl_passthru,ptpnl_qntile,ptpnl_passthru);
#' you can generate names for them automatically
#names(baz) <- make.unique(sapply(baz,attr,'fname'));
#' you can generate a logical vector indicating which of them should return a TRUE/FALSE decision
#sapply(baz,attr,'eval');
#' you can see how each of them were created
#sapply(baz,attr,'call');
#' Of course you can evaluate them all on the same dataset, the main purpose of them
logenv <- new.env();
for(ii in names(baz)) baz[[ii]](mtcars,c(-1,1),logenv=logenv,index=ii);

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
