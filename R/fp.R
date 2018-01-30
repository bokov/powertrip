# NOTE: censoring (cx and cy) is BROKEN!! Fix!

fp<-function(x,y=NULL,d=NULL,cx=NULL,cy=NULL,pars=NULL,cons=NULL,perms=NULL,pf=NULL,tlog=NULL,groupnames=NULL,checkinput=T,recalc=F,ynull=F,addpers=F,comparemodels=T,lx=length(x),ly=length(y), cullcons=T, modsig=NULL, modelchoose=c('manual','pChi','BIC','AIC','all'),verbose=0){
modelchoose<-match.arg(modelchoose);
maxdf<-4;

if(checkinput){
  if(modelchoose=='BIC'| modelchoose=='AIC'){
    warning('Sorry, BIC and AIC model selection not yet implemented. Defaulting to pChi model selection for now.'); 
    modelchoose<-'pChi';
  }
  
  switch(class(x)[1],
         list=,
         fp={d<-x; x<-NULL;},
         integer=,
         numeric={},
         data.frame=,
         array=,
         matrix=,
         stop('Support for multidimensional x variables not yet implemented')
         );
  if(is.null(d)) d<-list(); class(d)<-c('fp','list');
  if(!is.null(tlog)) d$tlog<-tlog;
  if(!is.null(groupnames)) d$groupnames<-groupnames;
  if(!is.null(x)) {d$x<-x;recalc<-T;}; if(!is.null(y)) {d$y<-y;recalc<-T}; 
  if(!is.null(cx)) {d$cx<-cx;recalc<-T;}; if(!is.null(cy)) {d$cy<-cy;recalc<-T;};
  if(!is.null(pars)) d$pars<-pars; if(!is.null(cons)) d$cons<-cons; 
  if(!is.null(perms)&(is.null(d$perms)||d$perms<perms)) {d$perms<-perms; addpers<-T;}
  if(!is.null(modsig)) d$modsig<-matrix(modsig,nrow=1,dimnames=list('Significance cutoff for models:',""));
  if(!is.null(pf)) d$pf<-pf; if(!is.null(tlog)) d$tlog<-tlog;
  if(is.null(d$x))
    stop('x is a required parameter that must be either a numeric vector or an fp object');
  if(is.null(d$y))
    stop('y is a required parameter that muse either be specified explicitly or inside an fp object');
  lx<-length(d$x); ly<-length(d$y); lxy<-lx+ly;
  if(is.null(d$pf)) d$pf<-mean;
  if(is.null(d$tlog)) d$tlog<-F;
  if(is.null(d$groupnames)) d$groupnames<-c('Control','Experimental');
  if(is.null(d$cx)) d$cx<-rep(1,lx);
  if(is.null(d$cy)) d$cy<-rep(1,ly);
  if(recalc){d$xy<-c(d$x,d$y); d$cxy<-c(d$cx,d$cy);}
  if(is.null(d$pars)) d$pars<-'iterative';
  if(is.null(d$modsig)) d$modsig<-matrix(0.05,nrow=1,dimnames=list('Significance cutoff for models:',""));
    # ought to convert a vector dpars to a list,
    # ought to make sure names(pars) == names(cons),
    # for now, though, just let it error?
  if(is.null(d$cons)|modelchoose=='all'){ 
   if(modelchoose=='manual'){modelchoose<-'pChi';}
    d$cons<-list(w=matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,0,0,1,1),
                   ncol=maxdf,byrow=T,
                   dimnames=list(NULL,par=c('a','b','c','s'))),
                 g=matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,0,0,1,1),
                   ncol=maxdf,byrow=T,
                   dimnames=list(NULL,par=c('a','b','c','s'))),
                 gm=matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,0,0,0,1),
                   ncol=maxdf,byrow=T,
                   dimnames=list(NULL,par=c('a','b','c','s'))),
                 l=matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,0,1,0),
                   ncol=maxdf,byrow=T,
                   dimnames=list(NULL,par=c('a','b','c','s'))),
                 lm=matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,0,0),
                   ncol=maxdf,byrow=T,
                   dimnames=list(NULL,par=c('a','b','c','s')))
		);
  } else if(length(d$cons)<2) comparemodels<-F;

    # if adaptive, get prereqs else if compare get prereqs
  omodels<-names(d$cons);
  dmodels<-modelinfo(omodels,'deps');
    # NOTE TO SELF: Watch cmodels... modelinfo's output has changed format!
  cmodels<-modelinfo(omodels);
  if(d$pars=='iterative'){
    addmodels<-dmodels;
  } else if(comparemodels) addmodels<-cmodels;

    # make sure all prereqs are met
  for(i in addmodels) if(is.null(d$cons[[i]])){
    d$cons[[i]]<-matrix(c(1,1,1,1),
                        ncol=maxdf,
                        dimnames=list(NULL,par=c('a','b','c','s')));
    }

    # if there is anything unsupported in d$cons, blow it away so it doesn't confuse the script later
  for(i in setdiff(names(d$cons),c('e','w','g','gm','l','lm')))
    d$cons[[i]]<-NULL;

  d$cons<-mapply(function(x,name){
	# make sure it's the right number of columns
    if(ncol(x)<maxdf){
      x<-cbind(x,matrix(1,nrow=nrow(x),ncol=maxdf-ncol(x)));
    } else if(ncol(x)>maxdf) {
      x<-x[,1:maxdf];warning('Too many columns in constraint, truncating.');
    }
	# make sure there is an unconstrained model and that it goes first
    if(!all(as.logical(x[1,]))) x<-rbind(1,x);
	# make sure that the unused params in each matrix = 1 
    x[,-match(modelinfo(name,'v'),c('a','b','c','s'))]<-1;
	# unique; named a,b,c,s; 
    x<-unique(x); colnames(x)<-c('a','b','c','s'); x;
    },d$cons,names(d$cons),SIMPLIFY=F); 
    if(is.null(d$perms)) d$perms<-500;
    if(is.null(d$pf)) d$pf<-mean;
    if(is.null(d$tlog)) d$tlog<-F;
}

if(is.null(d$seed)) d$seed<-round(unclass(Sys.time())*runif(3));
if(is.null(d$unc)|recalc){ 
  d$unc<-array(dim=c(6,6,2), dimnames=list(model=c('e','w','g','gm','l','lm'), data=c('MLE','cycles','a','b','c','s'), group=c('x','y')));
}

if(d$pars=='iterative'){
  cat('\nSelecting hazard models.\n');
  for(j in dimnames(d$unc)[[1]]){
    if(j == 'e'){
      d$unc['e',c('MLE','a'),'x']<-c(sum(log(1/mean(d$x))-d$x/mean(d$x)),1/mean(d$x));
      d$unc['e',c('MLE','a'),'y']<-c(sum(log(1/mean(d$y))-d$y/mean(d$y)),1/mean(d$y));
      next;
    }
	# find a starting parameter; if there are more than one candidate, take average
    compto<-modelinfo(j,ext=F);
    jxpars<-modpars(na.omit(d$unc[compto[1],c('a','b','c','s'),'x']), compto[1],j,pf=pf,trim=T,onegrp=T,nil=0);
    jypars<-modpars(na.omit(d$unc[compto[1],c('a','b','c','s'),'y']), compto[1],j,pf=pf,trim=T,onegrp=T,nil=0);
	# consider doing this inside modpars in the future
	# NOTE: this is a departure from how findpars() chooses starting parameters for lm;
	# in findpars, lm parameters are fitted using first the final parameters of gm as
	# the starting values (with s=0) and then again using the final parameters of l as
	# the starting values (this time with c=0). If the MLE for the gm vs lm comparison
	# that fitted the first set of candidate parameters is less than the MLE for the
	# l vs lm comparison that fitted the second set of parameters, the gm starting 
	# parameters are used for constrained fits; otherwise the l starting parameters are
	# used. Here the lm model is fitted only once, using an *average* of the final 
	# parameters from gm and l. I believe this is a very reasonable decision and one 
	# that may shorten the number of cycles needed, but it's still a source of
	# different results between Survomatic 1.3.0 and earlier and the versions
	# subsequent to this one, and may need to be documented.
    if(length(compto)==2) {
      jxpars<-(jxpars+modpars(na.omit(d$unc[compto[2],c('a','b','c','s'),'x']), compto[2],j,pf=pf,trim=T,onegrp=T,nil=0))/2;
	    jypars<-(jypars+modpars(na.omit(d$unc[compto[2],c('a','b','c','s'),'y']), compto[2],j,pf=pf,trim=T,onegrp=T,nil=0))/2;
    }
	# the magic number assigned to lb fixes an error, but it should really be fixed 
	# in the compiled library
    jxout<-opsurv(d$x,model=j,par=jxpars,tlog=d$tlog,cx=d$cx,lb=c(1e-15,0,0,0),verbose=verbose);
    jyout<-opsurv(d$y,model=j,par=jypars,tlog=d$tlog,cx=d$cy,lb=c(1e-15,0,0,0),verbose=verbose);
    d$unc[j,c('MLE','cycles','a','b','c','s'),'x']<- c(jxout$maximum,jxout$titer,modpars(jxout$estimate,j,trim=F)[1:4]);
    d$unc[j,c('MLE','cycles','a','b','c','s'),'y']<- c(jyout$maximum,jyout$titer,modpars(jyout$estimate,j,trim=F)[1:4]);
  } 
} else for(j in dimnames(d$unc)[[1]]){
# this is probably broken because later on the LRs are hard-coded, and not all the prerequisite functions may exist
  if(j == 'e'){
    d$unc['e',c('MLE','a'),'x']<-c(sum(log(1/mean(d$x))-d$x/mean(d$x)),1/mean(d$x));
    d$unc['e',c('MLE','a'),'y']<-c(sum(log(1/mean(d$y))-d$y/mean(d$y)),1/mean(d$y));
    next;
  }
  jout<-opsurv(d$x,model=j,par=d$pars[j,],tlog=d$tlog,cx=d$cx,lb=c(1e-15,0,0,0),verbose=verbose);
  d$unc[1,j,c('MLE','cycles','a','b','c','s'),1,1]<-
    c(jout$maximum,jout$titer,modpars(jout$estimate,j,trim=F)[1:4]);
  if(!ynull){
    jout<-opsurv(d$y,model=j,par=d$pars[j,],tlog=d$tlog,cx=d$cy,lb=c(1e-15,0,0,0),verbose=verbose);
    d$unc[1,j,c('MLE','cycles','a','b','c','s'),1,2]<-
      c(jout$maximum,jout$titer,modpars(jout$estimate,j,trim=F)[1:4]);
  }
}

# hardcoding, ugh; generalize someday
modtemp<-options('dmodels')[[1]][options('dmodels')[[1]]$compare==T,];
d$lrsy<-d$lrsx<-data.frame(Model=modtemp$complex, modlo=modtemp$simple, df=modtemp$df, LR=NA, Chi=NA, p=NA, sig=NA);

# calculate the log-ratios
for(i in unique(d$lrsx$Model)){
  d$lrsx[d$lrsx$Model==i,'LR']<-d$unc[modelinfo(i),'MLE','x'] - d$unc[i,'MLE','x'];
  d$lrsy[d$lrsy$Model==i,'LR']<-d$unc[modelinfo(i),'MLE','y'] - d$unc[i,'MLE','y'];
}

# also need to calculate the BIC and AIC when get around to it

d$lrsx$Chi<- -2*d$lrsx$LR; d$lrsy$Chi<- -2*d$lrsy$LR;
d$lrsx$p<-pchisq(d$lrsx$Chi,d$lrsx$df,low=F); d$lrsy$p<-pchisq(d$lrsy$Chi,d$lrsy$df,low=F);
d$lrsx$sig<-d$lrsx$p<as.numeric(d$modsig); d$lrsy$sig<-d$lrsy$p<as.numeric(d$modsig);

if(modelchoose=='pChi'){
    # based on LR's, choose best significant models
  lrsxtemp<-subset(d$lrsx,sig); lrsytemp<-subset(d$lrsy,sig);
  lrsxtemp<-choosemodel(paste(lrsxtemp$Model,lrsxtemp$modlo,sep='')); 
  lrsytemp<-choosemodel(paste(lrsytemp$Model,lrsytemp$modlo,sep=''));
    # select the simplest model that encompasses both groups
  d$jointmodels<-matrix(maxmodel(lrsxtemp,lrsytemp),nrow=1,dimnames=list('Models selected:',NULL));
}

if(modelchoose=='manual'|modelchoose=='all'){
  d$jointmodels<-matrix(maxmodel(names(d$cons),""),nrow=1,dimnames=list('Models selected:',NULL));
}


# have similar code-blocks for BIC and AIC

colnames(d$jointmodels)<-rep("",length(d$jointmodels));

if(modelchoose!='manual'){for(i in names(d$cons)){if(! i %in% d$jointmodels){d$cons[[i]]<-NULL;}}}

# d$perc is the array of unconstrained MLEs and parameter estimates to use for comparing with constrained MLES
# and as starting parameter values, respectively
d$perc<-array(dim=c(length(d$jointmodels),maxdf,2,d$perms), dimnames=list(model=as.character(d$jointmodels), params=c('a','b','c','s'), groups=c('x','y'), perms=NULL));
clbl=c('c','u');
# d$pers is the array of constrained log-ratios
d$pers<-array(dim=c(length(d$jointmodels),2,2,2,2,3,1+d$perms), dimnames=list(model=as.character(d$jointmodels), a=clbl, b=clbl, c=clbl, s=clbl, data=c('MLE','LR','df'), perms=NULL));

cat('\nFitting models to original data\n');
for(i in d$jointmodels){
    # sum and copy MLEs from d$unc into d$pers
  d$pers[i,2,2,2,2,'MLE',1]<-sum(d$unc[i,'MLE',]);
  ipars<-modpars(d$unc[i,c('a','b','c','s'),],modeli='lm',modelo=i,trim=T);
  for(j in 2:nrow(d$cons[[i]])){
    # don't know why I created dbg_jcons and dbg_dcons but getting rid of them to stop build warning
    # and because assignment to .GlobalEnv considered harmful
    jcons<-1+d$cons[[i]][j,]; #dbg_jcons<<-jcons; dbg_dcons<<-d$cons[[i]][j,]
    d$pers[i,jcons[1],jcons[2],jcons[3],jcons[4],'MLE',1]<-opsurv(d$x,d$y,i,ipars,jcons-1,cx=d$cx,cy=d$cy, tlog=d$tlog,verbose=verbose)$maximum;
	# populate the degrees of freedom for original AND all permutations, since that information repeats
    d$pers[i,jcons[1],jcons[2],jcons[3],jcons[4],'df',]<- maxdf-sum(jcons-1);
  }
}


if(d$perms>1&length(d$jointmodels)>0){
    set.seed(d$seed[3],kind="Mersenne-Twister",normal.kind="Inversion");
    xyperms<-apply(matrix(nrow=lxy,ncol=d$perms+1),2,function(z) sample(rep(c(T,F),c(lx,ly)),lxy));
    xyperms[,1]<-rep(c(T,F),c(lx,ly));

    # populate the unconstrained fits, using starting params from d$unc;
    cat('\nFitting unconstrained models to permuted data\n');
    for(i in d$jointmodels){
      ipars<-modpars(d$unc[i,c('a','b','c','s'),'x'],modeli='lm',modelo=i,trim=T);
      ipars<-(ipars+modpars(d$unc[i,c('a','b','c','s'),'y'],modeli='lm',modelo=i,trim=T))/2
      for(j in 2:(d$perms+1)){
        ijopsrv<-opsurv(d$xy[xyperms[,j]], d$xy[!xyperms[,j]], i, par=ipars, cx=d$cxy[xyperms[,j]], cy=d$cxy[!xyperms[,j]], verbose=verbose);
        d$perc[i,,,j-1]<-modpars(ijopsrv$estimate, modeli=i,pf=d$pf);
        d$pers[i,2,2,2,2,'MLE',j]<-ijopsrv$maximum;
      }
    }
        
    # populate the constrained fits, using starting params from d$perc[unconstrained]
    cat('\nFitting constrained models to permuted data\n');
    for(i in 2:(d$perms+1)){
      for(j in d$jointmodels){
        jpars<-modpars(as.numeric(d$perc[j,,,i-1]),modeli='lm',modelo=j,trim=T);
	    # the first row of each matrix in d$cons is assumed to be the unconstrained model
        for(k in 2:nrow(d$cons[[j]])){
          kcons<-1+d$cons[[j]][k,];
          d$pers[j,kcons[1],kcons[2],kcons[3],kcons[4],'MLE',i]<-opsurv(d$xy[xyperms[,i]], d$xy[!xyperms[,i]], j, jpars, kcons-1, cx=d$cxy[xyperms[,i]], cy=d$cxy[!xyperms[,i]], tlog=d$tlog,verbose=verbose)$maximum;
        }
      }
	# optionally write out data on each cycle; need to add a 'label' argument to function
    }
}

for(j in d$jointmodels){d$pers[j,,,,,'LR',]<-rep(d$pers[j,'u','u','u','u','MLE',],each=16)-d$pers[j,,,,,'MLE',];}

# OK all that remains now is to make the summary statistics table for the constraints (the statistics being pChi and pEmp)
#* OK remove unused fields from d$pers
# think about configurable d$pers
#* OK write summary function for fp
#* OK write plot function for fp
#* OK clean up opsurv (screen-spam off by default)
# clean up the argument parsing of FP
#* OK use fix() instead of tcl tables in survgui
#* integrate fp into survgui
#	OK update references to fp
#	move plotting to a menu
#	add new fp plotting functions

return(d);
}

