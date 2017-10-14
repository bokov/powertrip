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
  fit <- substitute(fit); result <- substitute(result); 
  if(!missing(eval)) eval <- substitute(eval);
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
  body(oo)[[2]][[3]][[2]][[2]] <- list(fname=fname,eval=!missing(eval));
  body(oo)[[3]] <- substitute(FIT <- try(fit));
  body(oo)[[3]][[2]] <- quote(fit);
  if(!missing(eval)) body(oo)[[4]][[3]][[3]] <- quote(return(NA));
  body(oo)[[4]][[4]][[2]] <- substitute(if(!is.null(logenv)) logenv[[index]] <- result);
  if(!missing(eval)) body(oo)[[5]] <- substitute(return(eval));
  oo;
}

#' Panels
ptpnl_passthru <- function(data,coords=NULL,...){
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
