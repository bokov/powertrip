qsrv <- function(y, x, q=.5) {
  tolerance <- .Machine$double.eps^0.5;
  q <- 1-q;
  keep <- (!is.na(y) & y < (q + tolerance));
  if (!any(keep)) NA else {
    x <- x[keep];
    y <- y[keep];
    if (abs(y[1] - q) < tolerance && any(y < y[1])){
      (x[1] + x[min(which(y < y[1]))])/2 } else x[1]
  }
}
##         med <- minmin(surv, time)
##         if (!is.null(upper)) {
##             upper <- minmin(upper, time)
##             lower <- minmin(lower, time)


smean <- function(x,cx=NULL,start.time=NULL,end.time=NULL,var=F,sd=F,q=F,ci=.95,cenlab=c('0','1')) {
  # x$time is just the unique time-points for x, ordered
  time <- sort(unique(x));
  # x$n.event is the number of non-censored events at each time point
  if(length(unique(cx))==1) cx <- NULL;
  if(!is.null(cx)){
    tc <- table(x,cx);
    n.event <- tc[,cenlab[2]];
    ntotal <- n.event+tc[,cenlab[1]];
  } else {
    ntotal <- n.event <- table(x);
  }
  # x$n.risk is the is n.event + the number of censored events at each time point
  n.risk <- rev(cumsum(rev(ntotal)));
  trisk <- pmax(n.risk,1);
  surv <- cumprod((trisk - n.event)/trisk);
  if(is.null(start.time)) start.time <- min(0,time);
  if(is.null(end.time)) end.time <- max(time);
  keep <- which(time <= end.time);
  hh <- c(ifelse((k <- n.event/n.risk)==1, 0, k/(n.risk - n.event))[keep],0);
  temptime <- c(time[keep], end.time);  n <- length(temptime);
  tempsurv <- c(surv[keep], surv[max(keep)]);
  # difference between each time point and the next
  delta <- diff(c(start.time, temptime))
  # area under the step-function, leaving out the final point
  rectangles <- delta * c(1, tempsurv[-n]);
  mean <- sum(rectangles) + start.time;
  if(!var&!sd&!q) return(mean) else{
    out <- c(mean=mean);
    varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1]);
    if(var) out <- c(out,var=varmean);
    if(sd) out <- c(out,sd=sqrt(varmean));
    if(any(q)){
      # usual stderrs, log conf intervals
      stderr <- sqrt(unlist(cumsum(n.event/(trisk * (trisk - n.event)))));
      zval <- qnorm(1 - (1 - ci)/2);
      xx <- pmin(surv,1);
      upper <- ifelse(surv==0,NA,pmin(exp(log(xx) + zval * stderr),1));
      lower <- ifelse(surv==0,NA,exp(log(xx) - zval * stderr));
      for(i in q){
        out <- c(out,lci=qsrv(lower,time,i),qsrv(surv,time,i),uci=qsrv(upper,time,i));
        names(out)[(length(out)-2):length(out)] <- paste(i*100,c('%, LCI','%','%, UCI'),sep='');
      }
    }
    return(out);
  }
}

              
  #if(is.null(start.time)) {
  #  if(!is.null(x$start.time)) start.time <- x$start.time else start.time <- min(0,x$time);
  #}
  #if(is.null(end.time)) end.time <- max(x$time);
  # find: n.risk, n.event
  # (all above are passed as arguments)
  # hh is some kind of hazard-- (# deaths / # at-risk) per surviving subject
  #hh <- with(x,ifelse((n.risk - n.event) == 0, 0, n.event/(n.risk * (n.risk - n.event))));
  #keep <- which(x$time <= end.time);
  # timepoints before a designated end-point
  #temptime <- c(x$time[keep], end.time);
  # corresponding survivorships for those timepoints
  #tempsurv <- c(x$surv[keep], x$surv[max(keep)]);
  # hazards for those timepoints
  #hh <- c(hh[keep], 0);
  # length of those vectors
  #n <- length(temptime);
  # difference between each time point and the next
  #delta <- diff(c(start.time, temptime))
  # area under the step-function, leaving out the final point
  #rectangles <- delta * c(1, tempsurv[-n]);
 # sum of cumulative sum of the squares of the rectangles multiplied by the corresponding hazards, in reverse order
##   varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1]);
##   mean <- sum(rectangles) + start.time;
##   if(!var&!sd) return(mean) else{
##     out <- c(mean=mean);
##     if(var) out <- c(out,var=varmean);
##     if(sd) out <- c(out,sd=sqrt(varmean));
##     return(out);
##   }
## }

## srvm <- function (x, scale = 1, rmean) 

##     if (!is.null(x$start.time)) start.time <- x$start.time else start.time <- min(0, x$time)
##     pfun <- function(nused, time, surv, n.risk, n.event, lower, 
##         upper, start.time, end.time) {
#         minmin <- function(y, x) {
#             tolerance <- .Machine$double.eps^0.5
#             keep <- (!is.na(y) & y < (0.5 + tolerance))
#             if (!any(keep)) 
#                 NA
#             else {
#                 x <- x[keep]
#                 y <- y[keep]
#                 if (abs(y[1] - 0.5) < tolerance && any(y < y[1])) 
#                   (x[1] + x[min(which(y < y[1]))])/2
#                 else x[1]
#             }
#         }
##         if (!is.na(end.time)) {
##             hh <- ifelse((n.risk - n.event) == 0, 0, n.event/(n.risk * 
##                 (n.risk - n.event)))
##             keep <- which(time <= end.time)
##             temptime <- c(time[keep], end.time)
##             tempsurv <- c(surv[keep], surv[max(keep)])
##             hh <- c(hh[keep], 0)
##             n <- length(temptime)
##             delta <- diff(c(start.time, temptime))
##             rectangles <- delta * c(1, tempsurv[-n])
##             #browser();
##             varmean <- sum(cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
##             mean <- sum(rectangles) + start.time
##         }
##         else {
##             mean <- 0
##             varmean <- 0
##         }
##         med <- minmin(surv, time)
##         if (!is.null(upper)) {
##             upper <- minmin(upper, time)
##             lower <- minmin(lower, time)
##             c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean), 
##                 sqrt(varmean), med, lower, upper)
##         }
##         else c(nused, max(n.risk), n.risk[1], sum(n.event), sum(mean), 
##             sqrt(varmean), med, 0, 0)
##     }
##     stime <- x$time/scale
##     if (is.numeric(rmean)) 
##         rmean <- rmean/scale
##     surv <- x$surv
##     plab <- c("records", "n.max", "n.start", "events", "*rmean", 
##         "*se(rmean)", "median", paste(x$conf.int, c("LCL", "UCL"), 
##             sep = ""))
##     ncols <- 9
##     if (is.null(x$strata)) {
##         if (rmean == "none") 
##             end.time <- NA
##         else if (is.numeric(rmean)) 
##             end.time <- rmean
##         else end.time <- max(x$time)
##         if (is.matrix(surv)) {
##             out <- matrix(0, ncol(surv), ncols)
##             for (i in 1:ncol(surv)) {
##                 if (is.null(x$conf.int)) 
##                   out[i, ] <- pfun(x$n, stime, surv[, i], x$n.risk, 
##                     x$n.event, NULL, NULL, start.time, end.time)
##                 else out[i, ] <- pfun(x$n, stime, surv[, i], 
##                   x$n.risk, x$n.event, x$lower[, i], x$upper[, 
##                     i], start.time, end.time)
##             }
##             dimnames(out) <- list(dimnames(surv)[[2]], plab)
##         }
##         else {
##             out <- matrix(pfun(x$n, stime, surv, x$n.risk, x$n.event, 
##                 x$lower, x$upper, start.time, end.time), nrow = 1)
##             dimnames(out) <- list(NULL, plab)
##         }
##     } else {
##         nstrat <- length(x$strata)
##         stemp <- rep(1:nstrat, x$strata)
##         last.time <- (rev(x$time))[match(1:nstrat, rev(stemp))]
##         if (rmean == "none") 
##             end.time <- rep(NA, nstrat)
##         else if (is.numeric(rmean)) 
##             end.time <- rep(rmean, nstrat)
##         else if (rmean == "common") 
##             end.time <- rep(median(last.time), nstrat)
##         else end.time <- last.time
##         if (is.matrix(surv)) {
##             ns <- ncol(surv)
##             out <- matrix(0, nstrat * ns, ncols)
##             if (is.null(dimnames(surv)[[2]])) 
##                 dimnames(out) <- list(rep(names(x$strata), rep(ns, 
##                   nstrat)), plab)
##             else {
##                 cname <- outer(dimnames(surv)[[2]], names(x$strata), 
##                   paste, sep = ", ")
##                 dimnames(out) <- list(c(cname), plab)
##             }
##             k <- 0
##             for (i in 1:nstrat) {
##                 who <- (stemp == i)
##                 for (j in 1:ns) {
##                   k <- k + 1
##                   if (is.null(x$lower)) 
##                     out[k, ] <- pfun(x$n[i], stime[who], surv[who, 
##                       j], x$n.risk[who], x$n.event[who], NULL, 
##                       NULL, start.time, end.time[i])
##                   else out[k, ] <- pfun(x$n[i], stime[who], surv[who, 
##                     j], x$n.risk[who], x$n.event[who], x$lower[who, 
##                     j], x$upper[who, j], start.time, end.time[i])
##                 }
##             }
##         }
##         else {
##             out <- matrix(0, nstrat, ncols)
##             dimnames(out) <- list(names(x$strata), plab)
##             for (i in 1:nstrat) {
##                 who <- (stemp == i)
##                 if (is.null(x$lower)) 
##                   out[i, ] <- pfun(x$n[i], stime[who], surv[who], 
##                     x$n.risk[who], x$n.event[who], NULL, NULL, 
##                     start.time, end.time[i])
##                 else out[i, ] <- pfun(x$n[i], stime[who], surv[who], 
##                   x$n.risk[who], x$n.event[who], x$lower[who], 
##                   x$upper[who], start.time, end.time[i])
##             }
##         }
##     }
##     if (is.null(x$lower)) 
##         out <- out[, 1:7, drop = F]
##     if (rmean == "none") 
##         out <- out[, -(5:6), drop = F]
##     list(matrix = out[, , drop = T], end.time = end.time)
## }




survsumm <- function(x,y=NULL,cx=NULL,cy=NULL,group=NULL,names=NULL,qs=c(.5,.9),nboot=1000,qmethod="PengHuang",tmax=NULL,quiet=F,adjust='holm'){
  if(is.null(y)&&is.null(group)) stop('The `y` and `group` arguments cannot both be left blank.');
  if(is.null(y)){
    xy <- x; cxy <- cx;
    x <- split(xy,group)[[1]]; y <- split(xy,group)[[2]];
    cx <- split(cxy,group)[[1]]; cy <- split(cxy,group)[[2]];
    group <- as.factor(group);
    if(is.null(names)) names <- levels(group);
  } else if(is.null(group)){
    if(is.null(names)) names <- c(deparse(substitute(x)),deparse(substitute(y)));
    group <- as.factor(rep(names,c(length(x),length(y)))); xy <- c(x,y);
  }
  if(is.null(tmax)) tmax <- max(xy);
  if(is.null(cx)) cx <- sign(x);
  if(is.null(cy)) cy <- sign(y);
  cxy <- c(cx,cy);
  xout <- c(N=length(x),censored=sum(1-cx),smean(x,cx,end.time=max(xy),sd=T));
  yout <- c(N=length(y),censored=sum(1-cy),smean(y,cy,end.time=max(xy),sd=T));
  pout <- c(N=NA,censored=NA,Mean=empnull(cbind(x,cx),cbind(y,cy),function(q,r,endt=tmax) {
    sq<-smean(q[,1],q[,2],end.time=endt,var=T);sr<-smean(r[,1],r[,2],end.time=endt,var=T);
    return((sq[1]-sr[1])/sqrt(sq[2]/nrow(q)+sr[2]/nrow(r)));
  },keepperm=T)$p,SD=NA);
  if(!is.null(qs)){
     crqout <- crq(Surv(xy,cxy)~group,taus=qs,method=qmethod);
    # slap on an extra quantile to work around a really silly limitation
    # of crq that causes it to output in a different format depending on
    # whether there is one tau or more than one
     if(quiet){
       junk <- capture.output(crqoutsumm <- try(summary(crqout,taus=c(.1,qs))));
     } else crqoutsumm <- try(summary(crqout,taus=c(.1,qs)));
     crqs <- {if(class(crqoutsumm)!='try-error') sapply(crqoutsumm,"[[",'tau') else {print(crqoutsumm);NA}};
     
    browser();
    # fix this to adapt to the actual quantiles crq/summary.crq can come up with
     qlabels <- apply(expand.grid(c('%','% LCI','% UCI'),qs*100)[2:1],1,paste,collapse='');
    for(i in qs){
      xout <- c(xout,smean(x,cx,q=i)[-1]);
      yout <- c(yout,smean(y,cy,q=i)[-1]);
      ip <- {if(j <- match(i,crqs,nomatch=0)) crqoutsumm[[j]]$coefficients[2,6] else NA} 
      #yout <- c(yout,quantile(y,iqs),qci(y,iqs));
      pout <- c(pout,NA,ip,NA);
    }
  }
  out <- cbind(xout,yout,pout);
  #
  #if(!is.null(qs)) rownames(out)[-(1:3)] <- qlabels;
  colnames(out) <- c(names,'P');
  #browser();
  if(!is.na(adjust)) out<-cbind(out,AdjP=p.adjust(out[,'P'],adjust)); 
  out <- list(Summary=out,cph=coxph(Surv(xy,cxy)~group)); out[['Log Rank']] <- summary(out$cph)$sctest; out;
  #list(Summary=out,`Log Rank`=surv2.logrank(Surv(xy,cxy),group));
}
