require(survival);

# doquant and the qsf hack below are temporary workarounds to bugs 15110 and 15111

doquant <- function (p, time, surv, upper, lower) 
{
    findq <- function(x, y, p) {
        if (all(is.na(y))) 
            return(rep(NA, length(p)))
        toss <- duplicated(y)
        if (any(toss)) {
            newy <- y[!toss]
            newx <- x[!toss]
        }
        else {
            newy <- y
            newx <- x
        }
        indx <- approx(c(0, 1 - newy), c(0:length(newy)), p)$y
        indx2 <- ceiling(indx)
        result <- newx[indx2]
        if (any(!is.na(indx) & indx == indx2)) {
            special <- which(indx == indx2)
            upper <- c(newx, max(x))[indx2[special] + 1]
            result[special] <- (result[special] + upper)/2
        }
        result
    }
    qq <- findq(time, surv, p)
    if (missing(upper)) 
        qq
    else rbind(qq, findq(time, lower, p), findq(time, upper, 
        p))
}

qsf <- survival:::quantile.survfit;
environment(qsf) <- .GlobalEnv;

# This is to be a replacment of the survsumm() function. Completely
# refactored, no longer using any contributed packages (but might
# implement the option to use crq later), and taking input in a more
# convenient form: xx is a vector of survival times, cc is a vector
# of right-censoring indicators, formula is a right-hand formula,
# data is a data.frame or similar, and qs is a vector of quantiles
# in the range of 0 to 1. In the future, might write a
# survsummng.survfit. Maybe the immediate future...

survsummng.survfit <- function(sf,data,qs,doquantbug=T,rmean='individual',...){
  # doquantbug needs to be set until the survival packages is patched
  # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15110
  lsf <- length(sf$strata);
  qs <- sort(qs);
  sflist <- sapply(names(sf$strata),
                   function(ii) {
                     oo <- sf[ii];
                     # Need to manually assign strata because a bug
                     # in quantile.survfit prevents individual strata
                     # from being analyzed and that in turn needlessly
                     # limits how many of the quantiles are done when
                     # some strata achieve the quantiles and others do
                     # not. If we don't assign strata, a single survfit
                     # object has no strata property.
                     # https://bugs.r-project.org/bugzilla3/show_bug.cgi?id=15111
                     oo$strata <- sf$strata[ii];
                     oo;
                   },
                   simplify = F);
  out <- do.call(rbind,
                 c(list(t(rbind(rbind(summary(sf, rmean = rmean, print.rmean = T)$table)[, -c(2:3, 7:9)]))),
                          ## rbind(summary(sf,
                          ##         rmean='common',
                          ##         print.rmean=T))$table[,-c(2:3,7:9)])),
                   lapply(qs,function(jj) sapply(sflist,function(ii) {
                     if(doquantbug){
                       ii$strata <- c(X=ii$n);
                       oo<-tryCatch(qsf(ii,jj),error=function(e) return(list(NA,NA,NA)));
                     } else {
                       oo<-tryCatch(quantile(ii,jj),
                                    error=function(e) return(list(NA,NA,NA)));
                     }
                     oo<-c(oo[[2]],oo[[1]],oo[[3]]);
                     dbg<-try(names(oo)<-paste(
                                  c('Lower Bound','Estimated','Upper Bound'),
                                  sprintf('%02.0f%% Mortality',jj*100)));
                     if(class(dbg)[1]=='try-error') browser();
                     oo;
                   }))));
  out <- rbind(out[1,,drop=F],censored=out[1,]-out[2,],out[-1,]);
  return(out);
}

survsummng.Surv <- function(sv,formula,data,qs,...){
  sv;
  newf<-update(formula,sv~.);
  sf <- survfit(newf,data=cbind(data,sv=sv));
  survsummng.survfit(sf,data,qs,...);
}

survsummng <- function(xx,formula,data,qs=c(.5,.9),cc=rep(1,length(xx)),...){
  switch(class(xx),
         survfit = survsummng.survfit(xx,data,qs,...),
         Surv = survsummng.Surv(xx,formula,data,qs,...),
         {
           xxerror <- 'The xx argument must be numeric and if specified the cc argument must be a vector of 0s and 1s equal in length to xx';
           if(is(xx,'character')){
             if(xx[1]%in%colnames(data)){
               warning('The xx argument given as a character. Obtaining event ages from data, but possibly ignoring censoring.');
               xx <- data[,xx[1]];
               if(length(cc)!=length(xx)) cc <- rep(1,length(xx));
             } else stop(xxerror);
           }
           if(length(cc)!=length(xx)||!is.numeric(xx)||any(!cc%in%0:1)) stop(xxerror);
           sv <- Surv(xx,cc);
           survsummng.Surv(sv,formula,data,qs,...);
         });         
}


survsummboot <- function(formula,data,qs=c(.5,.9),nperm=200,contnames=c('main.1','main.2','interaction'),
                         cnt=list(a=c(-1,1,0,0),b=c(0,0,-1,1),c=c(-1,1,1,-1)),
                         bootci=F,rowselect='^\\*rmean$|^Estimated',...){
  # formula = Surv(event,censor) ~ covariate + covariate + ... etc.
  #           Treats everything as a categorical covariate
  #          YOU DO NOT NEED TO SCREW AROUND WITH LHS AND RHS FORMULAS
  #          FOR THIS FUNCTION, JUST A STANDARD Surv(y,c)~x1+x2+... WILL DO
  # data = data frame
  # qs = quantiles of interest
  # nperm = number of permutations to use
  # contnames = character vector same length as cnt below
  # cnt = list of vectors in (-1,0,1), each the same length as the number of distinct combinations of covariates
  # bootci = use bootstrap confidence intervals
  # rowselect = search pattern to use for finding the values of interest
  # usage, two-group case:
  #  print(result01<-survsummboot(Surv(time, status) ~ sex,lung,contnames='Sex',cnt=list(c(-1,1))))
  # usage, two-way design case:
  #  lung$cage<-cut(lung$age,2); print(result02<-survsummboot(Surv(time, status) ~ cage+sex,lung,contnames=c('Age','Sex','Sex by Age')))
  dm <- model.frame(formula,data=data);
  nr <- nrow(dm);
  names(cnt) <- contnames;
  #perms <- replicate(nperm,sample(1:nr,nr,rep=T));
  cat('\nInitializing permutations.\n');
  perms <- replicate(nperm,do.call(ave,c(list(x=1:nr),dm[,-1],list(FUN=function(ii) {
    cat('.');sample(ii,length(ii),rep=T)}))));
  cat('\nPermutations ready.\n');
  orig<-survsummng.survfit(survfit(formula,data),data=data,qs,...);
  #orig <- survsummng.Surv(dm[,1],formula,dm,qs);
  #orig <- survsummng.Surv(dm[,1],formula,data,qs);
  rows <- grep(rowselect,rownames(orig));
  outnrows <- length(rows);
  outncols <- ncol(orig);
  origdiffs <- sapply(cnt,function(ii){
    rowSums(orig[rows,ii!=0]*matrix(ii,nrow=outnrows,ncol=outncols,byrow=T)[,ii!=0])});
  cat('\nStarting permutations:\n');
  sm <- sapply(1:nperm,function(ii) {
    cat('.');
    survsummng.survfit(survfit(formula,data[perms[,ii],]),data=data[perms[,ii]],qs,...);
    #survsummng.Surv(dm[perms[,ii],1],formula,dm[perms[,ii],],qs);
  },simplify='array');
  cat('\nPermutations complete.\n');
  diffs <- sapply(1:dim(sm)[3],function(jj) sapply(cnt,function(ii) {
    rowSums(sm[rows,ii!=0,jj]*matrix(ii,nrow=outnrows,ncol=outncols,byrow=T)[,ii!=0]);
  }),simplify='array');
  # Below, we test to make sure the differences aren't all NA or all the same
  # May need a more stringent criterion than that
  ps <- apply(diffs,1:2,function(ii) if(length(unique(na.omit(ii)))<2) NA else {
    2*min((1:0)+c(-1,1)*mean(ii>0,na.rm=T))
  });
  colnames(ps) <- paste('p ',names(cnt),sep='');
  # preparing output
  comparisons <- cbind(origdiffs,ps);
  if(bootci){
    # TODO: validation, like that done for p-values above;
    #       maybe make length(unique(na.omit(ii)))<2 a local function?
    # TODO: Also, probably more efficient to do prob=c(0.025,0.5,0.975) and then subset the object into orig
    orig[grep('Lower Bound ',rownames(orig)),]<-apply(sm[grep('Estimated ',dimnames(sm)[[1]]),,],
                                                      1:2,quantile,prob=0.025,na.rm=T);
    orig[grep('Upper Bound ',rownames(orig)),]<-apply(sm[grep('Estimated ',dimnames(sm)[[1]]),,],
                                                      1:2,quantile,prob=0.975,na.rm=T);
    orig[grep('Estimated ',rownames(orig)),]<-apply(sm[grep('Estimated ',dimnames(sm)[[1]]),,],
                                                    1:2,quantile,prob=0.5,na.rm=T);
  }
  out <- cbind(orig,matrix(NA,nrow=nrow(orig),ncol=ncol(comparisons),
                            dimnames=list(rownames(orig),colnames(comparisons))));
  out[rows,colnames(comparisons)]<-comparisons;
  format(data.frame(out,check.names=F),drop0trailing=T,quote=F,digits=5,na.encode=T);
  out;
}
