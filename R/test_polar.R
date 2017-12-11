#' ---
#' title: "Test Functions for PowerTrip"
#' author: "Alex F. Bokov"
#' date: "10/03/2017"
#' ---

#' Rescale cartesian parameters to span the space 
#' between maxs and mins (both of them vectors equal to
#' number of dimensions over which we are simulating)
#' xx is a matrix of cartesian coordinates with as many
#' columns as there are dimensions
#' Apparently we should just leave the ctr parameter alone
rescale <- function(xx,maxs=0,mins=0,ctr=(maxs+mins)/2,smaxs=1,smins=-1,sctr=0){
  sranges <- smaxs - smins;
  ranges <- maxs - mins;
  # first we center and scale to unit circle (cartesian)
  # then by multiplying by ranges and subtracting ctr we get
  # it to the new scale and center
  oo <- t(apply(xx,1,function(zz) ranges*(zz - sctr)/sranges+ctr));
}

#' Generates a random binomial value whose probability of being
#' 1 is a bivariate normal distribution scaled such that it is unity
#' at 0,0 and diminishes asymptotically in all directions from there
#' The purpose of this is to generate an underlying smooth, 
#' "radially monotonic" distribution for which only binomial outcomes
#' be observed, for testing powertrip functionality in isolation from
#' any specific statistical tests, etc.
#' @param xx numeric, the closer to 0, the more likely to return 1
#' @param yy numeric, the closer to 0, the more likely to return 1
gen_binom<-function(xx,yy){
  rbinom(1,1,mvtnorm::dmvnorm(c(xx,yy))/mvtnorm::dmvnorm(c(0,0)));
}

## convert cartesian to polar coordinates
## (why isn't this part of base R?!)
## returns a matrix of same dimension
crt2pol <- function(xx,center=kmeans(xx,1)$centers,sort=c('angle','radius','none')){
  ## xx = cartesian coordinates matrix [x,y]
  ## center = the center relative to which the polar coordinates will be returned
  ## sort = how/whether to sort the coordinates
  ##        (so the output matrix forms a closed shape when used by lines(...) )
  cxx<-scale(xx,center=center,scale=F);
  ## note the multiplication by sign--
  ## if not done, results will only be correct for x > center[1]
  oo<-cbind(r=sqrt(rowSums(cxx^2))*ifelse(sign(cxx[,1]),sign(cxx[,1]),1),theta=atan(cxx[,2]/cxx[,1]));
  ## retain the center for later use
  attr(oo,'center')<-center;
  ## preliminary testing suggests the below will work even in the degenerate case where
  ## all oo[,1] has the same sign, so don't worry about that for now
  oo[]<-switch(match.arg(sort),
               angle={
                 ooa<-oo[oo[,1]>=0,];oob<-oo[oo[,1]<0,];
                 ooa<-ooa[order(ooa[,2],decreasing=T),];
                 oob<-oob[order(oob[,2],decreasing=T),];
                 rbind(oob,ooa);},
               radius=oo[order(oo[,1]),],
               none=oo);
  oo;
}

## convert polar coordinates back to cartesian
## returns a matrix of same dimension
pol2crt.old <- function(pp
                    ,center=if('center'%in%names(attributes(pp))) as.numeric(attr(pp,'center')) else rep(0,len=ncol(pp))
                    #,maxes=Inf,mins=-Inf,denom=2
                    ){
  ## pp = polar coordinate matrix [r,theta]
  ## center = point on which the cartesian result should be centered
  oo<-scale(cbind(x=pp[,1]*cos(pp[,2]),y=pp[,1]*sin(pp[,2])),center= -center,scale=F);attr(oo,'center')<-center;
  # if(!isTRUE(abs(maxes)==Inf)&!isTRUE(abs(mins)==Inf)){
  #   # find the out-of-range values
  #   oor<-apply(oo,1,function(xx) any(xx>maxes)||any(xx<mins));
  #   oo[oor,]<-oo[sample(which(!oor),sum(oor),rep=T),];
  #   cat('\nReplacing ',sum(oor),' out-of-range values via resample')
  #   rngs <- (maxes - mins);
  #   cat('\nScaling...\n');
  #   oo <- t(apply(oo,1,function(xx) (xx-mins)/rngs));
  # }
  oo;
}

pol2crt <- function(coords,center=0,names){
  # Try to save some effort by accepting both matrices and vectors
  if(isTRUE(nrow(coords)>0)) {
    thisfun <- match.fun(as.character(sys.call()[1]));
    return(t(apply(coords,1,thisfun
              ,center=center,names=names)));
    } else {
      # we assume that the first coordinate is the radial one
      if(missing(names)) names<-make.names(seq_along(coords));
      setNames(coords[1] * 
                 # the rest are assumed to be thetas
                 c(cos(coords[-1]),1) *
                 cumprod(c(1,sin(coords[-1]))) + 
                 center,names);
    }
}

#' Title Find sampling limits on r for a given set of thetas
#' Returns the maximum value of 'r' such that upon conversion
#' to cartesian all coordinates will fall within the specified
#' limits. Note: 0 probably won't work as a limit
#'
#' @param coords  thetas from a set of polar coordinates
#' @param maxs    cartesian maximum values; indicate absence of limits with Inf
#' @param mins    cartesian minimum values; indicate absence of limits with -Inf
#' @param ... 
#'
#' @return        the maximum permitted value for 'r' at those thetas
#' @export
#'
#' @examples
pollim <- function(coords,maxs=Inf,mins=-Inf,...){
  if(length(maxs)!=length(coords)+1 && length(maxs)>1) stop('maxs should be 1 more than length of coords');
  if(length(mins)!=length(coords)+1 && length(mins)>1) stop('mins should be 1 more than length of coords');
  maxs <- rep_len(maxs,length(coords)+1); mins <- rep_len(mins,length(coords)+1);
  oo <- c(maxs,mins)/(c(cos(coords),1)*cumprod(c(1,sin(coords))));
  min(pmax(oo,0));
  # aha! hello and goodbye, bug! If oo is all negative, the below becomes Inf
  #min(oo[oo>0]);
}
#' Okay, the above works! The following generates a plot of radially sampled
#' points strictly bounded by the values of `maxs` and `mins`
#' which as of this moment are `c(4,3)` and `c(-2,-4)` respectively
# plot(
#   pol2crt(
#     cbind(
#       r=apply(baz
#               ,1,function(xx) runif(1,0,pollim(xx,maxs=maxs,mins=mins)))
#       ,theta=baz)));
#' Let's try iterating over the r's for one theta
sample_polar <- function(nn=1000,rmax=1,thetawrap=0.1
                         ,ntheta=1,center=c(0,0),maxs,mins){
  oo <- cbind(r=runif(nn,0,rmax),theta=runif(nn,0,(2+thetawrap)*pi));
  oocr <- rescale(pol2crt(oo)
                  ,maxs=maxs,mins=mins,smaxs=rmax,smins=-rmax);
  return(list(spol=oo,scrt=oocr));
}

#' This is like `sample_polar2()` but only generates valid samples, doesn't run
#' sims or tests
sample_polar3 <- function(thetas,nn=20
                          ,maxs=maxs,mins=mins
                          ,maxlims=pollim(thetas,maxs=maxs,mins=mins)
                          ,minlims=0
                          ,rpred=NA # here we will put the result of dose.p()
                          ,est=0,se=Inf,n.se=2
                          ,simfun=ptsim_binom
                          ,panfun=ptpnl_passthru
                          ,...){
  if(!missing(rpred)) {est<-rpred[1];se=attr(rpred,'SE')};
  lbound <- max(est-n.se*se,minlims,0);
  ubound <- min(est+n.se*se,maxlims);
  if(lbound>ubound) return(matrix(NA,nrow=0,ncol=length(thetas)+1));
  oo<-cbind(r=runif(nn,lbound,ubound)
            ,rbind(thetas)[rep_len(1,nn),]);
}
#' This is a row-wise operation
#' 
sample_polar2 <- function(thetas,nn=20
                          ,maxs=maxs,mins=mins
                          ,maxlims=pollim(thetas,maxs=maxs,mins=mins)
                          ,minlims=0
                          ,rpred=NA # here we will put the result of dose.p()
                          ,est=0,se=Inf,n.se=2
                          ,simfun=ptsim_binom
                          ,panfun=ptpnl_passthru
                          ,...){
  if(!missing(rpred)) {est<-rpred[1];se=attr(rpred,'SE')};
  lbound <- max(est-n.se*se,minlims,0);
  ubound <- min(est+n.se*se,maxlims);
  if(lbound>ubound) return(matrix(NA,nrow=0,ncol=length(thetas)+1));
  oo<-cbind(r=runif(nn,lbound,ubound)
            ,rbind(thetas)[rep_len(1,nn),]);
  if(!is.null(simfun)){
    data <- result <- list();
    oocr <- pol2crt(oo);
    for(ii in 1:nn){
      result[[ii]] <- panfun(simfun(oocr[ii,]),oocr[ii,],full=F);
    }
    # this might not always work, but it at least seems to work 
    oo<- cbind(oo,res=sapply(result,function(xx) xx[[1]]));
  }
  oo;
}

#' xx is a list with a polar and a cartesian version of the 
#' same coordinates
#' logenv is an environment that could contain an iter variable
#' if it doesn't, the variable is created and set to 1
simpoints <- function(xx,logenv
                      ,simfun=ptsim_binom,panelfun=ptpnl_passthru,...){
  if('iter' %in% ls(logenv)) myiter <- logenv$iter + 1 else myiter <- 1;
  rn <- rownames(xx$spol) <- rownames(xx$scrt) <- apply(
    xx$scrt,1,function(ii) paste0(c('_',ii),collapse='_'));
  names(rn) <- NULL;
  data <- sapply(rn,function(ii) simfun(xx$scrt[ii,],...)
                 ,simplify=F);
  logenv[[itername<-paste0('i',myiter)]] <- list();
  outcome<-c();
  for(ii in rn) {
    logenv[[itername]][[ii]]<-panelfun(data[[ii]],coords=xx$scrt[ii,]);
  }
  logenv$iter <- myiter;
  
  #spol <- rbind(spol,cbind(newsmp$spol,res=apply(newsmp$scrt,1,function(zz) gen_binom(zz[1],zz[2]))));
  cbind(xx$spol
        ,res=sapply(logenv[[itername]],function(ii) ii$outcome[1])
        ,iter=myiter
        );
}

#' select coordinates from a sample of the coordinate
#' space (xx) (in polar coordinates) using various criteria.
selcoords <- function(xx,withpreds=F
                      ,type=c('range','serange','quantile','topn')
                      ,target=0.8,se=1,quantile=0.1,range=target+c(-.1,.1),topn=100
                      ,model
                      ,...
){
  cols <- colnames(xx)[!colnames(xx) %in% c('res','iter')];
  # make predictions
  if(missing(model)) model <- glm(res~r*cos(theta)*sin(theta)
                                  ,data.frame(xx)
                                  ,family='binomial');
  preds <- predict(model,data.frame(xx),type='response',se.fit=T);
  nrmpreds <- (preds$fit - target)^2/preds$se.fit  
  type<-match.arg(type,several.ok = T);
  oo <- list();
  if('range' %in% type)
    oo$range <- preds$fit>range[1]&preds$fit<range[2];
  if('serange' %in% type)
    oo$serange <- abs(preds$fit-target)<se*preds$se.fit;
  if('quantile' %in% type)
    oo$quantile <- nrmpreds<quantile(nrmpreds,quantile);
  if('topn' %in% type)
    oo$topn <- rank(nrmpreds)<=topn;
  out <- apply(do.call(cbind,oo),1,all);
  if(withpreds) return(cbind(keep=out,preds=preds$fit,se.fit=preds$se.fit,normpreds=nrmpreds)) else return(out);
}

estimate_rs.new <- function(res,rr,glmfit,power=0.8,saveglm=F){
  if(missing(glmfit)) {init<-T; glmfit <- glm(res~rr,family='binomial') } else {
    init <- F; glmfit <- update(glmfit) };
  out <- MASS::dose.p(glmfit,p=power);
  out <- c(radest=as.vector(out[1])
           ,radse=attr(out,'SE')
           ,with(predict(glmfit,newdata=data.frame(rr=out[1]),type='resp',se.fit=T)
                 ,c(respest=fit,respse=se.fit)));
  if(init&&saveglm) attr(out,'glmfit') <- glmfit;
  out;
}

estimate_rs <- function(data,rlist=list(),power=0.8,...){
  if(!'iiglm' %in% rlist) {
    rlist$iiglm <- glm(res~r,data.frame(data),family='binomial');
  } else rlist$iiglm <- update(rlist$iiglm);
  rlist$iiprinv <- MASS::dose.p(rlist$iiglm,p=power);
  # TODO: check for rlist$iiprinv[1] being infinite
  rlist$iiprinv <- c(rlist$iiprinv[1],attr(rlist$iiprinv,'SE'));
  rlist$iipr <- with(predict(rlist$iiglm,newdata = data.frame(r=rlist$iiprinv[1])
                       ,type='resp',se.fit=T)
               ,c(fit,se.fit));
  rlist;
}

pt_plot <- function(dataenv,coords=c('cartesian','polar'),iter1col='red',iter1pch='+',...){
  trfun <- switch(match.arg(coords),cartesian=pol2crt,polar=identity); 
  dat<- cbind(trfun(cbind(dataenv$rs,dataenv$phis)),iters=dataenv$iters);
  plot(dat[,1:2],...);
  points(dat[dat[,'iters']==1,1:2],col=iter1col,pch=iter1pch)
}

samplephis <- function(dataenv,logenv=new.env(),errenv=new.env()
                       ,samples=200,dims=1,nworst=100
                       # at each phi, keep gathering simulations
                       # on sampled r's until tolse*se.fit < tol
                       ,power=0.8,tol=0.01,tolse=1
                       ,maxs=c(2,4.5),mins=c(-3.1,-1.3)
                       # named list of functions to run on each sim, the pt-panel
                       # ptpnl_qntile doesn't yet work
                       ,pnlst=list(lm=ptpnl_lm,lm2=update(ptpnl_lm,fname="lm2",frm=yy~(.)^2))
                       ,wd=paste0(getwd(),'/')
                       # when a file having the name specified by this variable is found, 
                       # the dataenv,logenv, and errenv objects are saved
                       ,savetrigger=paste0(wd,'pt_savedata')
                       # name of file to which to save when savetrigger encountered
                       ,savefile=paste0(wd,'pt_result.rdata')
                       # name of file that will be sourced (once and then moved) if found
                       ,sourcepatch=paste0(wd,'pt_sourcepatch.R')
                       ,...){
  tt <- proc.time();
  if(missing(dataenv)) dataenv <- new.env();
  # identify which functions on the panel list return T/F results vs. summary only
  pneval <- sapply(pnlst,attr,'eval');
  on.exit(save(dataenv,logenv,errenv,file=savefile));
  # create a matrix of randomly sampled angles 
  phis0 <- matrix(runif(samples*dims,0,2*pi),nrow=samples,ncol=dims);
  colnames(phis0) <- paste0('phi',seq_len(dims));
  # rownames(phis) <- NULL;
  # TODO: catch length mismatches and missing variables from dataenv
  if(length(setdiff(c('rs','phis','r_ses','phi_ses','iters','iilm','calls','times'),ls(dataenv)))==0){
    phipr <- predict(dataenv$iilm,data.frame(phis0),se.fit=T);
    phis0 <- phis0[phikeep <- which(rank(1-phipr$se.fit)<=nworst),,drop=F];
    phipr <- lapply(phipr,function(xx) xx[phikeep]);
  } else {
    dataenv$rs <- dataenv$phis <- dataenv$r_ses <- dataenv$phi_ses <- dataenv$cycle <- c();
    dataenv$iters <- 0;
    dataenv$calls <- dataenv$times <- list();
    phipr <- list(fit=rep_len(0,nrow(phis0)),se.fit=rep_len(Inf,nrow(phis0)));
    iter <- 1;
  }
  iter <- tail(dataenv$iters,1) + 1;
  #itername <- paste0('iter',iter);
  pb <- txtProgressBar(0,nrow(phis0),style=3);
  for(ii in 1:nrow(phis0)){
    iiphis <- phis0[ii,];
    iiname <- paste0('x_',paste0(iiphis,collapse='_'));
    maxlims0 <- pollim(iiphis,maxs=maxs,mins=mins);
    # here is where we may diverge from the current ptsim and ptpnl
    # This is a hardcoded test panel... the panel should be an argument
    iisims <- sample_polar2(iiphis,maxlims=maxlims0
                            ,est = phipr$fit[ii]
                            ,se = phipr$se.fit[ii]
                            # to make it just return coordinates with no results...
                            ,simfun = NULL 
                            ,n.se = 3);
    # here we need to iterate over iisims, for each row creating a sim and running
    # a panel on it
    failed <- F;
    cycle0 <- 0;
    jjres <- list();
    if((jjnn<-nrow(iisims))==0) {failed <- T;iilist<-NA} else {
      browser();
      for(jj in 1:jjnn){
        jjcoords <- iisims[jj,];
        jjname <- paste0('r_',paste(jjcoords,collapse='_'));
        jjdat <- ptsim_2lin(pol2crt(jjcoords));
        jjres[[jj]] <- sapply(pnlst[pneval],function(xx) any(
          xx(jjdat,jjcoords
             #,logenv=logenv
             ,errenv=errenv
             ,index=c(jjname,attr(xx,'fname')))));
      };
      browser();
      # here the test code stops, and the original code may break
      # do.call(rbind,jjres) should return a matrix where each column is a result
      # (yy) for a different one of the evaluating panel functions (i.e. have their
      # eval attribute set to TRUE) and the iisims[,'r'] column is the predictor
      # for a binomial model whose inverse prediction should approximate the values
      # needed for the desired detection rate (the inverse prediction is currently
      # handled by the estimate_rs() function)
      # initialize the glm models...
      sapply(data.frame(do.call(rbind,jjres)),estimate_rs.new,rr=foo$r,simplify=F) -> ban;
      # update them
      mapply(function(aa,bb) estimate_rs.new(aa,rr=foo$r,glmfit=attr(bb,'glmfit')),data.frame(do.call(rbind,jjres)),ban);
      iilist <- estimate_rs(iisims,power=power);
      # put next two lines inside sample_polar2?
      #colnames(iisims) <- c('r',colnames(phis),'res');
      #rownames(iisims) <- NULL;
      while(iilist$iipr[2]*tolse>tol&&!failed){
        # TODO: if the estimate is out of bounds, give up and move on
        iilist<-estimate_rs(
          iisims<-rbind(iisims,sample_polar2(iiphis,maxlims=maxlims0,rpred=iilist$iiprinv))
          ,power=power);
        cycle0 <- cycle0 + 1;
        if(cycle0 > 1 && (iilist$iiprinv[1]>maxlims0||iilist$iiprinv[1]<0)) failed <- T;
      }
    }

    # Now, actually retain the last result
    if(failed) {
      errenv[[iiname]] <- list(cycle=cycle0,phis=iiphis
                               ,preds=iilist,iter=iter);
      } else {
      dataenv$rs[iiname] <- iilist$iiprinv[1];
      dataenv$r_ses[iiname] <- iilist$iiprinv[2];
      dataenv$iters[iiname] <- iter;
      dataenv$cycle[iiname] <- cycle0;
      dataenv$phis <- rbind(dataenv$phis,iiphis);
      rownames(dataenv$phis)[ncol(dataenv$phis)] <- iiname;
    }
    setTxtProgressBar(pb,ii);
    if(file.exists(savetrigger)) {
      save(dataenv,logenv,errenv,file=savefile);
      file.remove(savetrigger);
    }
    if(file.exists(sourcepatch)) {
      source(sourcepatch,local = T);
      file.rename(sourcepatch,paste0(sourcepatch,'.bak'));
    }
  }
  close(pb);
  if('iilm' %in% ls(dataenv)) dataenv$iilm <- with(dataenv,update(iilm)) else {
    #environment(rformula) <- dataenv;
    #dataenv$iilm <- 
    with(dataenv,{
      phinames <- colnames(phis);
      rformula <- formula(paste0('rs~('
                                 ,paste0('sin('
                                         ,phinames
                                         ,')+cos('
                                         ,phinames,')'
                                         ,collapse='+')
                                 ,')^2'));
      iilm <- lm(rformula,data=data.frame(phis))});
  }
  dataenv$calls[paste0('iter',iter)] <- match.call();
  dataenv$times[paste0('iter',iter)] <- proc.time() - tt;
  return(dataenv);
}



#' ### Here we try it
#' 
#maxes <- c(1.5,1.8); mins<-c(-1.65,-2.3);
# maxs <- c(4,3); mins <- c(-2,-4);
#ctr <- c(0.3,-0.2); # okay, so the center parameter works...
#ctr <- c(0,0);
#' First, create the dataenv
# foo<-samplephis(maxs=maxs,mins=mins);
#' Sample from a polar space immediately converting to cartesian
#' TODO: create a sample_polar function
#thetas<-cbind(theta=runif(1000,0,2*pi));
#ii <- 3;
#ban <- (cbind(rlim=apply(thetas,1,pollim,maxs=maxs,mins=mins),theta=thetas));
#baz <- cbind(r=runif(100,0,ban[ii,'rlim']),theta=ban[ii,'theta']);
#baz <- cbind(baz,res=apply(baz,1,function(xx) ptsim_binom(pol2crt(xx))));
# baz <- sample_polar2(runif(1,0,2*pi));
# colnames(baz) <- c('r','theta','res');
# rownames(baz) <- NULL;
# bax <- glm(res~r,data.frame(baz),family='binomial');
# prdinv <- dose.p(bax,p=.8);
# prdr<-with(predict(bax,newdata = data.frame(r=prdinv[1])
#                    ,type='resp',se.fit=T)
#            ,c(fit,se.fit));
# while(prdr[2]*1>0.01){
#   # sim0 <- cbind(r=runif(20,max(0,sum(print(prdinv)*c(1,-2))),sum(print(prdinv)*c(1,2))),theta=ban[ii,'theta']);
#   # sim0 <- cbind(sim0,res=apply(sim0,1,function(xx) ptsim_binom(pol2crt(xx))));
#   baz <- rbind(baz,sample_polar2(baz[ii,'theta'],rpred=prdinv));
#   prdinv <- MASS::dose.p(bax<-update(bax),p=.8);
#   prdr<-with(predict(bax,newdata = data.frame(r=prdinv[1])
#                      ,type='resp',se.fit=T)
#              ,c(fit,se.fit));
# }
#' View the sampled points and the outcome
# dim(baz);
# plot(predict(bax,type='response')~baz[,'r'],col='red',ylim=c(0,1),pch='.');
# points(baz[,'r'],baz[,'res'],pch='+');
# abline(h=0.8,v=prdinv[1]+c(-2,2)*print(prdinv)[2],lty=2);
# spol <- simpoints(sample_polar(maxs=maxs,mins=mins,thetawrap=0.1),new.env());
#' Fit a model
#prmod <- glm(res~r*sin(theta)*cos(theta),data.frame(spol),family='binomial');
#twopi <- cbind(0,rep_len(2*pi,ncol(spol)-3),0,0);
# prmod <- glm(res~r*cos(theta)*sin(theta)
#              ,data.frame(rbind(spol
#                                # ,spol + twopi[rep_len(1,nrow(spol)),]
#                                # ,spol - twopi[rep_len(1,nrow(spol)),]
#                                ))
#              ,family='binomial');
# logenv <- new.env();
#' The following parts get repeated many times
# for(jj in 1:100){
#   newsmp <- sample_polar(nn=1000,rmax=1
#                          ,maxs=maxs,mins=mins,thetawrap=0.1);
#   bestfit <- selcoords(newsmp$spol
#                        ,type=c('quantile','range')
#                        ,model=update(prmod)
#                        ,quantile=.05
#                        );
#   newsmp <- sapply(newsmp,function(xx) xx[bestfit,],simplify=F);
#   spol <- rbind(spol,simpoints(newsmp,logenv));
# }
#' End repeat/update part
#' 
#' Below are the visualizations that can be done on any iteration
# plot(spol[,1:2],pch='.',col='#00000050'); #,xlim=c(1.5,2.5));
# points(spol[spol[,'iter']==max(spol[,'iter']),1:2],pch='+',col='red');
#' How good is the fitted contour?
#bar <- abs(predict(update(prmod),type='response')-0.2);
# foo<-rescale(pol2crt(spol[selcoords(spol,type='serange',se=2),1:2]),maxs=maxs,mins=mins);plot(foo,pch='.',col='#00000099');dim(foo);
#' Run the following once only after the first few iteration, for reference
# bar.bak <- bar; 
# foo.bak <- foo;
# points(foo.bak,col='red',pch='+');
#' How does this look on polar coordinates?
# plot(spol[spol[,'iter']<100&spol[,'iter']>0,1:2],pch='+',col='red',xlim=range(spol[,1],na.rm=T),ylim=range(spol[,2],na.rm=T));
# points(spol[selcoords(spol,type='quantile'),1:2],pch='*');
#' Note: if the distribution is not logistic, so far it's in a way
#' that does not cause it to be zero-inflated, overdispersed, or 
#' heteroscedastic
#' 
#' Let's see if this holds up when we are using real p-values.
#' But first we need to get bounds and centering to work.
#' 


#' Next steps...
#' 
#' Scaling just doesn't look like it will work. Or maybe that's
#' an artifact of the process not really being a statistical test?
#' 
#' The other alternative is to ...
#' 
#' 1. DONE Randomly sample from the thetas
#' 2. DONE At each point on thetas 
#'   1. DONE Exclude those going beyond limits, -RETAIN: smallest excluded r- (no longer necessary)
#'   2. DONE If too few left, resample within above limit
#'   3. DONE Simulate from 0 to 1. RETAIN: (best) result
#'   4. DONE Fit a univariate logistic model with r as predictor
#'   5. DONE Select all predicted values within 3*se.fit of 0.8
#'   6. DONE RETAIN: max and min r's corresponding to those values
#' 3. DONE At subsequent iterations...
#'   1. DONE At each r, predict the min-r and the smaller of max-r or exclusion zone
#'   2. DONE Now only simulate between these two points
#'   3. DONE RETAIN: result, smallest excluded, max-r, min-r
