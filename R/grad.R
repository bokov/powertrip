## automatically generate functions



sfninit<-function(){
  fns<-list(
         g=expression(log(aa*exp(bb*xx-aa*(exp(bb*xx)-1)/bb))),
         gm=expression(log((cc+aa*exp(bb*xx))*exp(-aa*(exp(bb*xx)-1)/bb-cc*xx))),
         l=expression(log(aa*exp(bb*xx)*(1+ss*aa*(exp(bb*xx)-1)/bb)^(-(ss+1)/ss))),
         lm=expression(log((cc+aa*exp(bb*xx)/(1+ss*aa*(exp(bb*xx)-1)/bb))*(exp(-cc*xx)*(1+ss*aa*(exp(bb*xx)-1)/bb)^(-1/ss))))
         );
  pars<-list(g=c(a='aa',b='bb'),gm=c(a='aa',b='bb',c='cc'),l=c(a='aa',b='bb',s='ss'),lm=c(a='aa',b='bb',c='cc',s='ss'));
  d1<-mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),fns,pars);
  d2<-mapply(function(pp,qq) outer(pp,qq,function(xx,yy) mapply(function(ii,jj) D(ii,jj),xx,yy)),d1,pars);
  grs<-sapply(d1,sapply,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii));oo;});
  hss<-sapply(d2,function(hh) apply(hh,1:2,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii[[1]]));oo;}));
  options(Survomatic_fns=fns); options(Survomatic_grs=grs); options(Survomatic_hss=hss);
}

hazsll<-list(g=expression(log(exp(kk)*exp(bb*xx))),
           gm=expression(log(cc+exp(kk)*exp(bb*xx))),
           l=expression(log((exp(kk)*exp(bb*xx))/(1+ss*exp(kk)*(exp(bb*xx)-1)/bb))),
           lm=expression(log(cc+(exp(kk)*exp(bb*xx))/(1+ss*exp(kk)*(exp(bb*xx)-1)/bb))));
srvsll<-list(g=expression(-exp(kk)*(exp(bb*xx)-1)/bb),
           gm=expression(-exp(kk)*(exp(bb*xx)-1)/bb-cc*xx),
           l=expression(log((1+ss*exp(kk)*(exp(bb*xx)-1)/bb)^(-1/ss))),
           lm=expression(log(exp(-cc*xx)*(1+ss*exp(kk)*(exp(bb*xx)-1)/bb)^(-1/ss))));


hazs<-list(g=expression(log(aa*exp(bb*xx))),
           gm=expression(log(cc+aa*exp(bb*xx))),
           l=expression(log((aa*exp(bb*xx))/(1+ss*aa*(exp(bb*xx)-1)/bb))),
           lm=expression(log(cc+(aa*exp(bb*xx))/(1+ss*aa*(exp(bb*xx)-1)/bb))));
srvs<-list(g=expression(-aa*(exp(bb*xx)-1)/bb),
           gm=expression(-aa*(exp(bb*xx)-1)/bb-cc*xx),
           l=expression(log((1+ss*aa*(exp(bb*xx)-1)/bb)^(-1/ss))),
           lm=expression(log(exp(-cc*xx)*(1+ss*aa*(exp(bb*xx)-1)/bb)^(-1/ss))));
pars<-list(g=c(a='aa',b='bb'),gm=c(a='aa',b='bb',c='cc'),l=c(a='aa',b='bb',s='ss'),lm=c(a='aa',b='bb',c='cc',s='ss'));
parsll<-list(g=c(a='kk',b='bb'),gm=c(a='kk',b='bb',c='cc'),l=c(a='kk',b='bb',s='ss'),lm=c(a='kk',b='bb',c='cc',s='ss'));
## d1h<-mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),hazs,pars);
## Holy cow! The below is a simpler version because D() doesn't always simplify properly
## Not only that, but the extra complexity may be causing rounding errors
## The xx^0 is necessary so that the sum of reciprocals is returned instead of just one 
d1h<-c(list(g=alist(a=xx^0/aa,b=xx)),mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),hazs[-1],pars[-1]));
d1s<-mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),srvs,pars);
d2h<-mapply(function(pp,qq) outer(pp,qq,function(xx,yy) mapply(function(ii,jj) D(ii,jj),xx,yy)),d1h,pars);
d2s<-mapply(function(pp,qq) outer(pp,qq,function(xx,yy) mapply(function(ii,jj) D(ii,jj),xx,yy)),d1s,pars);
grh<-sapply(d1h,sapply,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii));oo;});
grs<-sapply(d1s,sapply,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii));oo;});
hsh<-sapply(d2h,function(hh) apply(hh,1:2,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii[[1]]));oo;}));
hss<-sapply(d2s,function(hh) apply(hh,1:2,function(ii) {oo<-function(xx,aa,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii[[1]]));oo;}));

d1hll<-mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),hazsll,parsll);
d1sll<-mapply(function(ii,jj) sapply(jj,function(kk) {D(ii,kk)}),srvsll,parsll);
d2hll<-mapply(function(pp,qq) outer(pp,qq,function(xx,yy) mapply(function(ii,jj) D(ii,jj),xx,yy)),d1hll,parsll);
d2sll<-mapply(function(pp,qq) outer(pp,qq,function(xx,yy) mapply(function(ii,jj) D(ii,jj),xx,yy)),d1sll,parsll);
grhll<-sapply(d1hll,sapply,function(ii) {oo<-function(xx,kk,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii));oo;});
grsll<-sapply(d1sll,sapply,function(ii) {oo<-function(xx,kk,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii));oo;});
hshll0<-sapply(d2hll,function(hh) apply(hh,1:2,function(ii) {oo<-function(xx,kk,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii[[1]]));oo;}));
hssll0<-sapply(d2sll,function(hh) apply(hh,1:2,function(ii) {oo<-function(xx,kk,bb,cc,ss){};body(oo)<-substitute(sum(jj),env=list(jj=ii[[1]]));oo;}));


hssll<-list(
         g=function(xx,kk,bb,cc,ss) {
           ebx<-exp(bb*xx); ebx1kb<-(ebx-1)*exp(kk)/bb; ebxekx <- ebx*exp(kk)*xx;
           matrix(c(
                    sum(-ebx1kb),rep(sum(-(ebxekx-ebx1kb)/bb),2),sum(-(bb*ebxekx*xx-2*ebxekx+2*ebx1kb)/bb^2)
                    ),nrow=2)
           ## matrix(c(sum(-(exp(bb*xx+kk)-exp(kk))/bb),
           ##          rep(sum(-((bb*exp(kk)*xx-exp(kk))*exp(bb*xx)+exp(kk))/bb^2),2),
           ##          sum(-((bb^2*exp(kk)*xx^2-2*bb*exp(kk)*xx+2*exp(kk))*exp(bb*xx)-2*exp(kk))/bb^3)),nrow=2)
         },
         gm=function(xx,kk,bb,cc,ss) {
           ekebx1b <- (exp(kk)*(exp(bb*xx)-1))/bb;
           xebxkb <- (xx*exp(bb*xx+kk))/bb;
           dkdb <- ekebx1b/bb - xebxkb; 
           matrix(c(sum(-ekebx1b),sum(dkdb),0,
                    sum(dkdb),sum(-xx*xebxkb+2*xebxkb/bb-2*ekebx1b/bb^2),0,
                    0,0,0),nrow=3);
         },
         l=function(xx,kk,bb,cc,ss) {matrix(nrow=3,ncol=3)},
         lm=function(xx,kk,bb,cc,ss) {matrix(nrow=4,ncol=4)}       
         );

hshll<-list(g=function(xx,kk,bb,cc,ss) matrix(0,nrow=2,ncol=2),
            gm=function(xx,kk,bb,cc,ss) {
              ebxk <- exp(bb*xx+kk); ebxkc2 <- (ebxk+cc)^-2;
              ebxkebxkc2 <- ebxk*ebxkc2;
              dkdk <- cc*ebxkebxkc2; dkdb <- dkdk*xx; dkdc <- -ebxkebxkc2; dbdb <- dkdb*xx; dbdc <- -ebxkebxkc2*xx; dcdc <- -ebxkc2;
              matrix(c(sum(dkdk),sum(dkdb),sum(dkdc),sum(dkdb),sum(dbdb),sum(dbdc),sum(dkdc),sum(dbdc),sum(dcdc)),nrow=3);
            },
            l=function(xx,kk,bb,cc,ss) {matrix(nrow=3,ncol=3)},
            lm=function(xx,kk,bb,cc,ss) {matrix(nrow=4,ncol=4)}
            );

survgrad<-function(pars,xx,ce=T,model='g',aa=pars[1],bb=pars[2],cc=ifelse(length(pars)>2,pars[3],0),ss=tail(pars,1)){
  mapply(function(ii,jj) ii(xx,aa=aa,bb=bb,cc=cc,ss=ss)+jj(xx[as.logical(ce)],aa=aa,bb=bb,cc=cc,ss=ss),grs[[model]],grh[[model]]);
}

survhess<-function(pars,xx,ce=rep(1,len(xx)),model='g',aa=pars[1],bb=pars[2],cc=ifelse(length(pars)>2,pars[3],0),ss=tail(pars,1),sr=hss,hz=hsh){
  apply(hss[[model]],1:2,function(ii) {ii[[1]](xx,aa,bb=bb,cc=cc,ss=ss)})+apply(hsh[[model]],1:2,function(jj) {jj[[1]](xx[as.logical(ce)],aa,bb=bb,cc=cc,ss=ss)});
}

confint.survomatic <- function(object,parm,data,level=0.95,xx,ce,models='available',...){
  ## object = an object fitted by survomatic
  ## parm = for which parameters should confidence intervals be reported? optional, all by default, currently ignored
  ## level = confidence level, optional
  ## xx,ce = event times and censoring indicator. optional, if missing looks in object, if fails, error
  ## models = models for which to return ci's; if 'available' (default), returns cis for all models available in object
  ##          ...but currently ignored, and CIs calculated for 'g'
  if(missing(data)) data<-object$data;
  zz <- qnorm((1-level)/2,lower=F);
  mf<-model.frame(object$formula,data=data);
  ## to prevent error for univariate cases
  if(ncol(mf)==1) mf$`_X` <- 1;
  ## when gm, l, and lm CIs ready, fix the hardcoded 'g' below
  prs<-lapply(object$unconstrained.fits[['g']],function(ii) c(log(ii$par[1]),ii$par[2:4]));
  mapply(function(vv,ww) hessllconfint(xx=vv[,1],ce=vv[,2],pars=ww,zz=zz,model='g',alreadyLoggedLambda=T),split(mf[,1],mf[,-1]),prs,SIMPLIFY=F);
  ## extract data from object
  ## identify which models needed
  ## extract params from object lapplying over models
  ## mapply hessllconfint
}

## low-level function for calculating confidence intervals on a single dataset, using log transform of lambda
hessllconfint <- function(xx,ce,zz=qnorm(.975),pars,model='g',alreadyLoggedLambda=F){
  ## xx, ce = event times and censoring indicator, respectively
  ## zz = the number of standard errors corresponding to the desired conf.int (.95 two-tailed by default)
  ## pars = parameter vector of length 4, as extracted from a survomatic object
  ## model = string scalar, 'g','gm','l', or 'lm'
  ## alreadyLoggedLambda = if true, assume the lambda parameter is already log transformed, otherwise (default) log it
  if(!alreadyLoggedLambda) pars[1]<-log(pars[1]);
  ses<-sqrt(diag(-qr.solve(
                    hssll[[model]](xx,kk=pars[1],bb=pars[2],cc=pars[3],ss=pars[4])+
                    hshll[[model]](xx[ce!=0],kk=pars[1],bb=pars[2],cc=pars[3],ss=pars[4]))));
  lower<-upper<-pars; lower[!is.na(pars)]<-lower[!is.na(pars)]-zz*ses; upper[!is.na(pars)]<-upper[!is.na(pars)]+zz*ses;
  pars[1]<-exp(pars[1]); lower[1]<-exp(lower[1]); upper[1]<-exp(upper[1]);
  oo <- cbind(est=pars,lower=lower,upper=upper);
  attr(oo,'ses')<-ses; attr(oo,'zz')<-zz; attr(oo,'model')<-model;
  oo;
}

## WARNING: Censoring ignored at this time!
## BUT can be fixed by having separate loghazard and logsurv fns, doing log surv on all
## and loghaz on just the non-censored and adding the sums for the final grad or hess

## Returns a gradient for ONE population and set of parameters
survgrad.bak<-function(pars,xx,model='g',aa=pars[1],bb=pars[2],cc=ifelse(length(pars)>2,pars[3],0),ss=tail(pars,1)){
  ## xx = event times
  ## ...because if the above are not in pars, and the pars are correct, the model won't need them anyway
  ## cc<-ss<-0;
  ## switch(model,gm={cc<-pars[3]},l={ss<-tail(pars,1)},lm={cc<-pars[3];ss<-pars[4]});
  sapply(getOption('Survomatic_grs')[[model]],function(ii) ii(xx,aa=aa,bb=bb,cc=cc,ss=ss));
}

## Returns a Hessian for ONE population and set of parameters
survhess.bak<-function(pars,xx,model=c('g','gm','l','lm'),aa=pars[1],bb=pars[2],cc=ifelse(length(pars)>2,pars[3],0),ss=tail(pars,1)){
  ## ...because if the above are not in pars, and the pars are correct, the model won't need them anyway
  ## cc<-ss<-0;
  ## switch(model,gm={cc<-pars[3]},l={ss<-tail(pars,1)},lm={cc<-pars[3];ss<-pars[4]});
  apply(getOption('Survomatic_hss')[[model]],1:2,function(ii) {ii[[1]](xx,aa=aa,bb=bb,cc=cc,ss=ss)});
}

## Uses survhess to calculate standard errors for ONE population and set of parameters
survhse<-function(pars,xx,ce=rep(1,length(xx)),model=c('g','gm','l','lm'),cens=T,slvfn=qr.solve,sr=hss,hz=hsh) {
  hss<-substitute(hss);hsh<-substitute(hsh);
  if(cens) sqrt(diag(slvfn(-survhess(pars,xx,ce,model,sr=sr,hz=hz)))) else sqrt(diag(slvfn(-survhess.bak(pars,xx,model))));
}

demogng <- function(frm,data){
  with(survfit(frm,data),{
    group<-rep(names(strata),strata);
    px<-do.call(c,by(surv,group,function(ii) ii/c(NA,ii)[seq_along(ii)]));
    ux<-ifelse(px==0,NA,-log(px));
    data.frame(time=time,deaths=n.event,censored=n.censor,atrisk=n.risk,
               SE=std.err,lx=surv,px,ux,
               lnux=ifelse(ux%in%c(-Inf,Inf),NA,log(ux)),
               group=gsub('^.*=','',group))
  });
}
## if foo is the output from survwrapperng, then the following will work whether or not keepdata=T
## demogng(foo$formula,eval(foo$data))

print.survomaticng<-function(xx,what=c('fits','hypotheses'),models='best',...){
  print(summary.survomaticng(xx,what=what,models=models,...));
  invisible(xx);
}

summary.survomaticng<-function(xx,what=c('modelsearch','hypotheses','fits','params','aic','ll'),models=c('all','best','chosen'),...){
  what<-match.arg(what,several.ok=T); models<-match.arg(models);
  oo<-list();
  if(!xx$force&&'modelsearch'%in%what) oo$Model.Search<-setNames(data.frame(do.call(rbind,xx$uncLRT),xx$modSigs),c('Chisq','df','p'));
  showmodels<-switch(models,all=names(xx$unconstrained.fits),best=xx$models.selected[1],chosen=xx$models.selected);
  if('fits'%in%what) {
    oo$Fits<-list();
    for(ii in showmodels) {
      iifits<-sapply(xx$unconstrained.fits[[ii]],`[[`,'par');iinc<-ncol(iifits);
      iises<-if(is.null(xx$SEs[[ii]])) matrix(NA,nrow=4,ncol=iinc) else apply(xx$SEs[[ii]],2,modparsng,min=ii,mout='lm',pad=NA);
      rownames(iifits)<-c('a','b','c','s'); dimnames(iises)<-dimnames(iifits);
      colidx<-rep(seq.int(iinc),each=2)+rep(c(0,iinc),len=2*iinc);
      oo$Fits[[ii]]<-data.frame(iifits,SE=iises)[,colidx];
      if(!xx$hypotheses) names(oo$Fits[[ii]])<-c('X','X.SE');
    }
  }
  if(xx$hypotheses&&'hypotheses'%in%what){
    oo$Hypothesis.Tests<-list();
    for(ii in if('best'%in%models) showmodels else names(xx$partially.constrained.fits)) {
      iiht<-t(sapply(xx$partially.constrained.fits[[ii]],function(jj) c(jj$LL,length(na.omit(jj$par))-ncol(jj$par)+1)));
      if(ii %in% names(xx$constrained.fits)) iiht<-rbind(iiht,any=c(xx$constrained.fits[[ii]]$LL,length(na.omit(xx$constrained.fits[[ii]]$par))));
      ## iiunc<-do.call(`+`,lapply(xx$unconstrained.fits[[ii]],function(kk) c(kk$LL,length(na.omit(kk$par)))));
      ## above: another >2 groups bug
      iiunc <- rowSums(sapply(xx$unconstrained.fits[[ii]],function(kk) c(kk$LL,length(na.omit(kk$par)))));
      iiht<-t(apply(iiht,1,`-`,iiunc))%*%-diag(2:1);
      iiht<-cbind(iiht,pchisq(iiht[,1],iiht[,2],lower=F)); colnames(iiht)<-c('Chisq','df','P');
      oo$Hypothesis.Tests[[ii]]<-iiht;
    }
  }
  if('params'%in%what){
    oo$Parameters<-lapply(xx$unconstrained.fits[showmodels],sapply,`[[`,'par');
  }
  if('ll'%in%what||'aic'%in%what){
    iill<-cbind(t(sapply(xx$unconstrained.fits,sapply,function(ii) c(LL=ii$LL, DF=iidf<-sum(!is.na(ii$par)),AIC=2*(-ii$LL+iidf)))));
    ## browser();
    colnames(iill)<-apply(expand.grid(c('LL','DF','AIC'),xx$groups)[,2:1],1,paste,collapse='.');
    oo$LL.and.AIC<-cbind(iill,
                         jnt.LL=rowSums(iill[,seq(1,ncol(iill),by=3),drop=F]),
                         jnt.DF=rowSums(iill[,seq(2,ncol(iill),by=3),drop=F]),
                         jnt.AIC=rowSums(iill[,seq(3,ncol(iill),by=3),drop=F]));
  }
  oo;
}
## TODO: Make generic, so it can use matrices, vectors, etc
## DONE: maybe implement forcemodel?
## DONE: retain the names and order of predictor groups
## TODO: retain sample size
## TODO: try nlminb for optimization (someday)
## TODO: implement weibull?
## TODO: implement fake exponential
## TODO: more customizable constraints
## TODO: summary.survomatic, print.survomatic, print.summary.survomatic, plot.survomatic
##       ...with minimalist/batch summary options
## Here it is!!!
survwrapperng<-function(frm,data,testmodels=NULL,forcemodels=NULL,doses=c('selectedonly','all','none'),keepdata=F,loglambda=F,pscale,...){
  ## frm = formula (with Surv(...) on the left side
  ## data = corresponding data.frame
  ## testmodels = optional character vector of models on which to test hypotheses even if they are not selected by the model search
  ## forcemodels = optional character vector of the ONLY models on which to test hypotheses; this overrides testmodels and bypasses model selection
  ## doses = which SEs to calculate: selectedonly is only for the models where hypothesis tests were done
  ##                                 all is for all (unconstrained) models
  ##                                 none is to skip ses
  ## keepdata = whether to return the actual data in the output or just the call to it
  ## Note: used to have a subset option, then renamed it to sbs, but both broke model.frame in envsetup, so gave up
  ## Just manually subset the data if necessary
  if(class(frm)[1]=='call') frm <- eval(frm);
  uenv<-envsetup(frm,data,cons=c(T,T),model='g',loglambda=loglambda); uncLRT<-list(); force <- length(forcemodels)>0;
  deps<-c('g','gm','l');
  if(force){
    testmodels<-forcemodels; deps<-unique(c(testmodels,if('lm'%in%testmodels) deps else 'g'));}
  sftemp<-survfit(frm,data);
  pgstart<- lapply(seq.int(max(1,length(sftemp$strata))),function(ii) with(sftemp[ii],{
    pars<-lm(lnux~time,data.frame(lnux=log(-log(c(surv[-1],NA)/surv)),time),subset=!is.infinite(lnux))$coef;
    if(!loglambda) pars[1]<-exp(pars[1]); setNames(pars,c('a','b'))
  }));
  if(missing(pscale)) pscale<-c(if(!loglambda) 1e-6 else .1,1e-3,1e-6,1e-2);
  ## So, here we do the log-linear regression for initial 'g' params the data comes ready-made as formulae!
  ##g_startx<-with(survfit(Surv(x,cx)~1),lm(lnux~time,subset(data.frame(lnux=log(-log(c(surv[-1],NA)/surv)),time),!is.infinite(lnux)))$coef)
  ## Then, fit g models on that series of params and uenv
  umods<-list(g=mapply(function(ii,jj) opsurvng(ii,mpars=jj,pscale=pscale[1:2]),uenv,pgstart,SIMPLIFY=F));
  ## Update uenv to gm, modpars the g params to gm, and fit gm models
  if('gm'%in%deps){
    pgmstart <- lapply(umods$g,function(ii) modparsng(na.omit(ii$par[,1]),min='g',mout='gm'));
    updateuncenv(uenv,'gm',mvars=c('a','b','c'),loglambda=loglambda);
    umods$gm<-mapply(function(ii,jj) opsurvng(ii,mpars=jj,pscale=pscale[1:3]),uenv,pgmstart,SIMPLIFY=F);
    if(!force){
      ## uncLRT$gm.vs.g<-(do.call(`+`,lapply(umods$gm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
      ##                  do.call(`+`,lapply(umods$g,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      uncLRT$gm.vs.g<-(rowSums(sapply(umods$gm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
                       rowSums(sapply(umods$g,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
    }
  }
  ## Update uenv to l, modpars the g params to l, and fit on l models
  if('l'%in%deps){
    plstart <- lapply(umods$g,function(ii) modparsng(na.omit(ii$par[,1]),min='g',mout='l'));
    updateuncenv(uenv,'l',mvars=c('a','b','s'),loglambda=loglambda);
    umods$l<-mapply(function(ii,jj) opsurvng(ii,mpars=jj,pscale=pscale[c(1:2,4)]),uenv,plstart,SIMPLIFY=F);
    if(!force){
      ## uncLRT$l.vs.g<-(do.call(`+`,lapply(umods$l,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
      ##                 do.call(`+`,lapply(umods$g,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      uncLRT$l.vs.g<-(rowSums(sapply(umods$l,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
                       rowSums(sapply(umods$g,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
    }
  }
  ## Update uenv to lm, modpars (the g? the gmean of gm and l?) to lm, and fit on lm models
  ## TODO: conditional statement for whether to evaluate the lm model
  if(!force) modSigs<-sapply(uncLRT,function(ii) pchisq(ii[1],ii[2],lower=F));
  if('lm'%in%c(testmodels,deps)||(!force&&any(modSigs<.05))){
    plmstart <- mapply(function(jj,kk) exp((log(jj)+log(kk))/2),
                       lapply(umods$gm,function(hh) pmax(1e-16,hh$par[,1],na.rm=T)),
                       lapply(umods$l,function(ii) pmax(1e-16,ii$par[,1],na.rm=T)),SIMPLIFY=F);
    updateuncenv(uenv,'lm',mvars=c('a','b','c','s'),loglambda=loglambda);
    umods$lm<-mapply(function(ii,jj) opsurvng(ii,mpars=jj,pscale=pscale),uenv,plmstart,SIMPLIFY=F);
    if(!force){
      ## uncLRT$lm.vs.gm<-(do.call(`+`,lapply(umods$lm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
      ##                   do.call(`+`,lapply(umods$gm,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      ## uncLRT$lm.vs.l<-(do.call(`+`,lapply(umods$lm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
      ##                  do.call(`+`,lapply(umods$l,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      uncLRT$gm.vs.g<-(rowSums(sapply(umods$lm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
                       rowSums(sapply(umods$gm,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      uncLRT$gm.vs.g<-(rowSums(sapply(umods$lm,function(ii) c(ii$LL,sum(!is.na(ii$par))))) -
                       rowSums(sapply(umods$l,function(ii) c(ii$LL,sum(!is.na(ii$par))))))*2:1;
      modSigs<-c(modSigs,sapply(tail(uncLRT,2),function(ii) pchisq(ii[1],ii[2],lower=F)));
    }
  }
  ## Model selection done here; at the moment hardcoded
  ## modSigs<-sapply(uncLRT,function(ii) pchisq(ii[1],ii[2],lower=F))<.05;
  if(!force){
    names(modSigs)<-gsub("\\..+","",names(modSigs));
    modelsSelected<-unique(names(modSigs)[modSigs<.05]);
    ## Note: unique seems to be faster than union
    if(length(intersect(modelsSelected,c('gm','l')))==0) modelsSelected<-'g' else {
      if(length(modelsSelected)>1) modelsSelected<-unique(c('lm',modelsSelected))};
    modelsSelected<-unique(c(modelsSelected,testmodels));
  } else modelsSelected <- forcemodels;
  ## ses disabled-- the right way to get dispersion information from now on is to do confint on this object
  ## the stub below is temporarily here to prevent breaking other stuff
  ses <- list();
  ## ses<-switch(match.arg(doses),
  ##             selectedonly=sapply(modelsSelected,function(ii) mapply(function(jj,kk)
  ##               tryCatch(with(kk$data[[1]],survhse(na.omit(jj$par[,1]),xx=srv[,1],ce=srv[,2],model=ii)),
  ##                        error=function(ee) rep(NA,length(na.omit(jj$par[,1])))),umods[[ii]],uenv),simplify=F),
  ##             all=sapply(names(umods),function(ii) mapply(function(jj,kk)
  ##               tryCatch(with(kk$data[[1]],survhse(na.omit(jj$par[,1]),xx=srv[,1],ce=srv[,2],model=ii)),
  ##                        error=function(ee) rep(NA,length(na.omit(jj$par[,1])))),umods[[ii]],uenv),simplify=F),
  ##             none=list());
  ## TODO: check to make sure there are more than 1 groups before doing hypothesis tests
  ## make constraint matrix for param differences
  oo<-list(formula=frm,data=if(!keepdata) substitute(data) else data, survfit=if(keepdata) sftemp else NULL, force=force,
           ## the formatversion number is the date when the format being used here was implemented
           ## this is to dispatch obsolete versions of the survomatic object to the obsolete methods
           SEs=ses,models.selected=modelsSelected,formatversion=13061700,call=sys.call(),
           ## below used to say sapply(umods,names)[,1] ... in univariate case that is an error... will this now be an error in multi?
           unconstrained.fits=umods,groups=sapply(umods,names,simplify=F)[[1]]);
  ## browser();
  if(!force) {oo$uncLRT<-uncLRT; oo$modSigs<-modSigs;}
  if(length(uenv)>1){
    pcons<-1-diag(4); dimnames(pcons)<-list(c('a','b','c','s'),c('a','b','c','s')); cumods <- list();
    ## Do modelinfo on each of modelsSelected, capture pnames, do pcons[pnames,pnames] when building env
    ## TODO: see if this will perform better vectorized
    selmodvars<-sapply(modelsSelected,modelinfo,'vars',simplify=F);
    for(ii in modelsSelected) {
      iivars<-selmodvars[[ii]]; iicons<-pcons[iivars,iivars]; iistart<-na.omit(unlist(lapply(umods[[ii]],`[[`,'par')));
      for(jj in seq.int(nrow(iicons))) {
        if(exists('cuenv')) envsetup(cuenv,cons=iicons[jj,],model=ii,loglambda=loglambda) else {
          cuenv<-envsetup(frm,data,cons=iicons[jj,],model=ii,loglambda=loglambda);
        }
        ## browser();
        cumods[[ii]][[iivars[jj]]]<-opsurvng(cuenv[[1]],pars=iistart);}};
    ## Now, fully constrained models
    ## cstarts<-lapply(umods[modelsSelected],function(jj) c(na.omit(exp(do.call(`+`,lapply(jj,function(ii) log(ii$par)))/length(jj)))));
    ## above was breaking under > 2 groups, below is an attempt to fix
    if(loglambda){
      cstarts <- lapply(umods[modelsSelected], function(jj) c(na.omit(rowSums(sapply(jj,function(ii) ii$par/length(jj))))));
    } else {
      cstarts <- lapply(umods[modelsSelected], function(jj) c(exp(na.omit(rowSums(sapply(jj,function(ii) log(ii$par)/length(jj)))))));
    }
    cmods<-list();
    for(ii in modelsSelected){
      iivars<-selmodvars[[ii]]; iistart<-cstarts[[ii]];
      if(exists('cenv')) envsetup(cenv,cons=rep(0,len=length(iivars)),model=ii,loglambda=loglambda) else {
        cenv<-envsetup(frm,data,cons=rep(0,len=length(iivars)),model=ii,loglambda=loglambda);
      }
      cmods[[ii]]<-opsurvng(cenv[[1]],pars=iistart);
    }
    oo$hypotheses<-T;
    oo$partially.constrained.fits<-cumods; oo$constrained.fits<-cmods;
  } else oo$hypotheses<-F;

  attr(oo,'class')<-'survomatic';
  oo;
  ## the below was for testing
  ## gvars<-modelinfo('g','vars'); gcons<-pcons[gvars,gvars];gstart<-na.omit(unlist(lapply(umods$g,`[[`,'par')));
  ## gmvars<-modelinfo('gm','vars'); gmcons<-pcons[gmvars,gmvars];gmstart<-na.omit(unlist(lapply(umods$gm,`[[`,'par')));
  ## lvars<-modelinfo('l','vars'); lcons<-pcons[lvars,lvars];lstart<-na.omit(unlist(lapply(umods$l,`[[`,'par')));
  ## lmvars<-modelinfo('lm','vars'); lmcons<-pcons[lmvars,lmvars];lmstart<-na.omit(unlist(lapply(umods$lm,`[[`,'par')));
  ## cuenv<-envsetup(frm,data,cons=gcons[1,],model='g');
  ## if(length(uenv)>1) ...
  ## testmodels <- union(best joint model, testmodels)
  ## do hypothesis tests on all of testmodels (note the best model is always first)
  ## probably by lapplying over a list of constraints to generate a list of envs,
  ## generating a list of gmean-ed mpars, and lapplying opsurvng over those
  ## remember to split the globbed fully constrained estimates out
  ## ...and then extracting LL's as before, no need to sum this time
  ## Maybe have the various opsurvng outputs have additional model info for convenience? Such as next/prev?
  ## return unconstrained lists, and list of lists of constrained hypothesis test outputs
  ## ...after adding a version vector and setting class
  ## Example of how to get SEs, but might be a bug in the survhse code
  ## mapply(function(jj,kk) {jj$SE<-survhse(na.omit(jj$par[,1]),xx=kk$srv[,1],ce=kk$srv[,2],model='lm');jj},umods$lm,uenv[[1]]$data,SIMPLIFY=F)
}

## To save a little time, updating specifically environment lists for unconstrained models
updateuncenv<-function(env,model,lb=c(1e-14,0,0,0),ub=c(.5,.5,.5,.5),mvars=modelinfo(model,'vars'),loglambda=F){
  ## env = environment previously created by envsetup
  ## model = the usual g, gm, l, lm character string
  ## should really inherid loglambda from env if not specified, but below will hopefully work for now
  if(loglambda) {lb[1]<-log(lb[1]); ub[1]<-log(ub[1]);}
  mlen<-length(mvars); mseq<-setNames(seq_along(mvars),mvars); ofun<-paste0('of',model); ub<-ub[mseq]; lb<-lb[mseq];
  idx<-modparsng(mseq,min=model,mout='lm',pad=NA);
  lapply(env,function(ii) {ii$model<-model;ii$ofun<-ofun; ii$ub<-ub; ii$lb<-lb; if(loglambda) ii$ll <- T;
                           ii$cons<-matrix(mseq,ncol=1,dimnames=list(mvars,names(ii$data[1])));
                           ii$data[[1]]$idx<-idx; ii});
}

## takes a Surv() formula, data, a constraint model, and a hazard model name ('g','gm','l','lm')
## and returns an environment for use with future objf and multisurvgrad
## NOT KNOWN what the hell happens if you do something obnoxious like use strata/cluster/frailty
## Suggest that you DON'T, and don't use numeric predictors... pretend this is survdiff or survfit
envsetup <- function(frm,data,cons,model,lb=c(1e-14,0,0,0),ub=c(.5,.5,.5,.5),loglambda=F){
  ## frm = either a formula, an environment previously returned by envsetup, or a list of such environments
  ##       if a formula, must have Surv(...) on the left side; either no censoring or right-censored only
  ## data = data from which the formula gets its variables
  ## cons = 1/0 , T/F or matrix with nrow = number of params, ncol = number of groups, passed to valcons
  ## model = 'g','gm','l','lm'
  ## DONE: look for redundant columns in cons and glob corresponding models
  ## DONE: look for non-overlapping groups of columns and split into separate envs
  ## DONE: have the idx/srv object be a list within an environment that also stores cons and some kind of tracking indexes for globbed groups
  ## DONE: update objfng and multisurvgrad to the new env format
  ## below: if loglambda is specified, reparametrize 
  if(loglambda) {lb[1]<-log(lb[1]);ub[1]<-log(ub[1]);}
  lb<-modparsng(lb,'lm',model); ub<-modparsng(ub,'lm',model);
  if(class(frm)[1]=='call') frm <- eval(frm);
  if(class(frm)[1]=='formula'){
    oo <- model.frame(frm,data); if(ncol(oo)==1) oo$X<-1;
    oo <- split(oo[,1,drop=F],oo[,-1]);
    if(all(cons)){
      ## If it's a fully unconstrained model, don't bother with the general stuff, just split it
      mvars<-modelinfo(model,'vars'); mlen<-length(mvars); mseq<-setNames(seq_along(mvars),mvars);
      ofun<-paste0('of',model);ub<-ub[mseq]; lb<-lb[mseq];  idx=modparsng(mseq,min=model,mout='lm',pad=NA);
      ## For some lbs and ubs where they are not all equal, the following may break for the 'l' model
      return(sapply(names(oo),function(ii) as.environment(list(model=model,ofun=ofun,ub=ub,lb=lb,
                                                               cons=matrix(mseq,ncol=1,dimnames=list(mvars,ii)),ll=if(loglambda) T else NULL,
                                                               data=setNames(list(list(srv=oo[[ii]],idx=idx,size=nrow(oo[[ii]]))),ii)))));
    }
    cons <- valcons(cons,names(oo),model);
    ## split the groups that have nothing to do with each other
    splits<-unique(apply(cons,2,function(ii) apply(cons==(replicate(ncol(cons),ii)),2,any)));
    splits<-lapply(seq(nrow(splits)),function(ii) seq(ncol(splits))[splits[ii,,drop=F]]);
    oo<-lapply(splits,function(ii) list(data=oo[ii],cons=cons[,ii,drop=F]));
    ## glob together groups that are identical, and save a toglob object for returning the correct number and order of estimates
    ## unless there is only group to begin with, in which case return just that object
    oo<-lapply(oo,function(ii) with(ii,if(ncol(cons)>1) {
      globm<-unique(apply(cons,2,function(jj) apply(cons==replicate(ncol(cons),jj,simplify='matrix'),2,all)));
      toglob<-lapply(seq(nrow(globm)),function(jj) setNames(seq(ncol(globm))[globm[jj,]],colnames(globm)[globm[jj,]]) );
      names(toglob)<-lapply(toglob,function(jj) paste0(colnames(globm)[jj],collapse='_'));
      gout<-list(data=lapply(toglob,function(jj) do.call(rbind,data[jj])),cons=cons[,sapply(toglob,`[`,1),drop=F],toglob=toglob);
      gout;
    } else ii));
    ## refresh each split-group's cons object, insert idx objects into the data, and make each split-group into an environment
    ## browser();
    oo<-lapply(oo,function(ii) {
      ii$cons<-valcons(ii$cons,groups=names(ii$data),model=model);
      if(loglambda) ii$ll <- T;
      ii$data<-mapply(function(jj,kk) list(
                                        srv=jj,idx=modparsng(ii$cons[,kk],min=model,mout='lm',pad=NA),
                                        size=nrow(jj)),
                      ii$data,seq(ncol(ii$cons)),SIMPLIFY=F);
      ii$lb=expandpars(rep(lb,ncol(ii$cons)),ii$cons,loglambda=!is.null(ii$ll));ii$ub=expandpars(rep(ub,ncol(ii$cons)),ii$cons,loglambda=!is.null(ii$ll));
      ii$ofun<-paste0('of',model); ii$model<-model; if(loglambda) ii$ll <- T;
      as.environment(ii);});
    return(oo);
  } else if(is.environment(frm)) oo <- list(frm) else if(is.list(frm)) oo<-frm;
  ## If we have been passed an environment or list of environments, we do the below instead
  oo<-lapply(oo,function(ii){
    ii$cons<-oocons<-valcons(cons,groups=names(ii$data),model=model);
    for(jj in seq.int(ncol(oocons))) ii$data[[jj]]$idx<-modparsng(oocons[,jj],min=model,mout='lm',pad=NA);
    ii$lb<-expandpars(rep(lb,ncol(oocons)),oocons,loglambda=!is.null(ii$ll)); ii$ub<-expandpars(rep(ub,ncol(oocons)),oocons,loglambda=!is.null(ii$ll));
    if(ii$model!=model) {ii$ofun<-paste0('of',model); ii$model<-model;}
  });
}

## function to convert a single-group starting param vector to the number of groups and constraints specified by a valid cons matrix (from valcons)
## here too is a hack to accomodate loglambda-- when loglambda is true, ALL parameters get arithmetic mean instead of geometric mean
expandpars<-function(pars,cons,loglambda=F){
  if(loglambda) unlist(lapply(split(pars,cons),function(ii) mean(ii))) else {
    unlist(lapply(split(pars,cons),function(ii) exp(mean(log(ii)))));
  }
  ## mc<-max(cons);
  ## oo<-rep(NA,len=mc); 
  ## for(ii in seq.int(ncol(cons))) oo[cons[,ii]]<-pars;
  ## oo;
}

## simple function for combining together Surv objects
rbind.Surv<-function(...){
  targets<-list(...);
  oo<-do.call(rbind,lapply(targets,unclass));
  type<-unique(unlist(sapply(targets,attr,'type')));
  if(is.null(type)||length(type)>1) warning('Censor type mismatch when binding Surv objects');
  attr(oo,'type')<-type[1];
  class(oo)<-'Surv';
  oo;
}

## New objf!
## Now sometimes works. Need to figure out why for test data says "non-finite value supplied by optim"
## Okay, DOES work, once valid starting pars are found
## Needs to be aware of upper and lower bounds
## TODO: implement cutoffs, and compare to opsurv with two groups
## TODO: test
## Outputs a VECTOR of LL's one for each group, should be wrapped in sum() if used with optim()
objfng <- function(pars,env){
  ## pars, env = passed by optim and generated by envsetup as usual
  ## ofun = 'ofg','ofgm','ofl','oflm'
  if(any(pars<env$lb|pars>env$ub)) return(-Inf);
  if(!is.null(env$ll)) ltrans <- exp else {
    ## to permit log(lambda), a.k.a. kk
    ltrans <- identity;
  }
  sapply(env$data,function(ii) with(ii,{try(.C(env$ofun,
                                                 a=as.double(ltrans(pars[idx[1]])),
                                                 b=as.double(pars[idx[2]]),
                                                 c=as.double(if(is.na(idx[3])) 0 else pars[idx[3]]),
                                                 d=as.double(if(is.na(idx[4])) 0 else pars[idx[4]]),
                                                 x = as.integer(srv[,1][,1]), # am I breaking something else by doing this?
                                                 size = as.integer(size),
                                                 cens = as.integer(srv[,1][,2]), # am I breaking something else by doing this?
                                                 ans = double(1),
                                                 PACKAGE = "Survomatic")$ans)}));
}

## A wrapper for calling survgrad from optim
## TODO: test returning the right sized value given new cons format
multisurvgrad<-function(pars,env){
  ## pars = numeric vector passed by optim
  ## env = environment as created by envsetup
  ## model = character scalar: 'g','gm','l', or 'lm'
  ## below: if the ll flag in the env is set, instead of using the lambda parameter, use exp(log(kk))
  if(!is.null(env$ll)) grfun <- survgradll else grfun <- survgrad;
  grtemp <- sapply(env$data,function(ii) with(ii,try(grfun(pmax(pars[idx],0,na.rm=T),srv[,1],srv[,2],model=env$model))));
  constemp <- na.omit(env$cons);
  oo <- pars;
  for(ii in seq_along(oo)) oo[ii]<-sum(grtemp[constemp==ii]);
  oo;
}

## The long-awaited replacement for opsurv.
## Hopefully a minimalist wrapper for optim
opsurvng<-function(env,pars,mpars=expandpars(pars,env$cons,loglambda=!is.null(env$ll)),
                   method1='Nelder-Mead',method2=c('BFGS','Nelder-Mead'),pscale,debug=F){
  ## pars = param vector appropriate to the model specified within env
  ## env = environment created by envsetup
  ## pscale = added to ctrl which is in turn passed to optim's control argument
  ## method1 = optim method to use on first iteration
  ## method2 = a character vector of TWO methods to alternate between for subsequent iterations
  ##           to not alternate, make them the same
  ## mpars = auto-calculated in the arg, but pre-generated version can be passed for batch mode, etc.
  ##
  ## initialize ctrl and functions
  if(missing(pscale)) {
    pscale<-c(if(is.null(env$ll)) 1e-6 else .1,1e-3,1e-6,1e-2);
    pscale<-expandpars(rep(pscale,ncol(env$cons)),env$cons,loglambda=!is.null(env$ll));
  }
  fn<-function(pp) sum(objfng(pp,env)); gr<-function(pp) multisurvgrad(pp,env);
  ctrl<-getOption('ctrl'); ctrl$parscale<-pscale;
  ## initialize loop
  oo<-optim(mpars,fn,gr,method=method1,control=ctrl);
  oldval<- -Inf; whichm2 <- c(T,F);
  ## if keeplog, run the slow version else, just run the loop
  ## Deliberately repetitive, so that the debug logic doesn't slow down the non-debug loop
  if(debug) {
    oplog<-list(oo);
    while(oldval<oo$value){
      oldval<-oo$value;
      oo<-optim(oo$par,fn,gr,method=method2[whichm2],control=ctrl);
      whichm2<-!whichm2;
      oplog<-c(oo,oplog);
    }
    attr(oo,oplog)<-oplog; return(oo);
  } else while(oldval<oo$value) {
    oldval<-oo$value;
    oo<-optim(oo$par,fn,gr,method=method2[whichm2],control=ctrl);
    whichm2<-!whichm2
  };
  ## return the parts that matter
  return(list(LL=oo$value,par=sapply(env$data,function(ii) oo$par[ii$idx])));
}
