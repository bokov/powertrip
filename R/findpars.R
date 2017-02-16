# former nil and bnil defaults: nil=6e-10,bnil=0.005,
`findpars`<-
function(x, y=NULL, cx=NULL, cy=NULL, nil=0, bnil=0, wbnil=1, pf=mean, label=NULL, summary=F, id=0, tlog=F, digits=22, sig=.05, models=NULL, writeunc=F, writegml=F, verbose=1){
	if(is.null(models)){models<-c('e','w','wa','wb','wab','g','ga','gb','gm','gma','gmb','gmc',
				      'l','la','lb','ls','gmlm','llm','lm','lma','lmb','lmc','lms');}
	xemax<-sum(log(1/mean(x))-x/mean(x));
	yemax<-if(is.null(y)){xemax;}else{sum(log(1/mean(y))-y/mean(y));}
	if(models[1]=='nocons'| is.null(y)){
		modcons<- grep('[w,g,gm,l,lm][a,b,c,s]',models);
		if(length(modcons)>0) {models<-models[-modcons]; }
	}
	# I've been returning double the maximum in the one-group case; the if statement
	# below should now fix it, but might also introduce new bugs
	es.e<-list(maximum=if(is.null(y)){xemax}else{sum(xemax,yemax)},
		   estimate=c(1/mean(x),1/mean(if(is.null(y)) x else y)),group0=list(maximum=xemax),
		   group1=list(maximum=yemax),titer=0,runtime=0);
	if(is.null(y)){onegrp=T;xrep=1;} else {onegrp=F;xrep=2;}
	partemp<-modpars(es.e$estimate,'e','g',nil=nil,trim=T);
	parw<-partemp[c(2,4)]<-bnil;
	es.w<-es.wa<-es.wb<-es.wab<-es.g<-es.ga<-es.gb<-es.gab<-es.gm<-es.gma<-es.gmb<-es.gmc<-es.gmabc<-
	es.l<-es.la<-es.lb<-es.ls<-es.labs<-es.gmlm<-es.llm<-es.lm<-es.lma<-es.lmb<-es.lmc<-es.lms<-es.lmabcs<-
	list(maximum= -Inf,titer=NA,group0=list(maximum=-Inf),group1=list(maximum=-Inf));
	es.wab$estimate<-es.gab$estimate<-rep.int(-Inf,2);
	es.wa$estimate<-es.wb$estimate<-es.ga$estimate<-es.gb$estimate<-es.gmabc$estimate<-es.labs<-rep.int(-Inf,3);
	es.gma$estimate<-es.gmb$estimate<-es.gmc$estimate<-es.la$estimate<-es.lb$estimate<-
	es.ls$estimate<-rep.int(-Inf,5);
	es.gmlm$estimate<-es.llm$estimate<-es.lmabcs$estimate<-rep.int(-Inf,4*xrep)
	es.lma$estimate<-es.lmb$estimate<-es.lmc$estimate<-es.lms$estimate<-rep.int(-Inf,7);
	if(any(c('w','wa','wb','wab') %in% models) ){
		# apparently it's normal to have a huge b value for weibull
		# in fact, below a certain 0 value, weibul errors out; therefore
		# for weibull only we substitute in an empirically discovered 
		# value that does not cause errors; should eventually be moved
		# to the ofw C function, though wbnil=1 generates errors where
		# it should be identical to the exponential model
		parw<-partemp; parw[c(2,4)]<-wbnil;
		es.w<-opsurv(x,y,'w',par=parw,tlog=tlog,ub=c(.5,10,.5,10),cx=cx,cy=cy, verbose=verbose);
		#browser();
		if(match('wa',models,nomatch=0)){
			es.wa<-opsurv(x,y,'w',cons=c(F,T,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.w$estimate,'w',cno=c(F,T,T,T),
						pf=pf,trim=T),tlog=tlog,ub=c(.5,10,.5,10));
			if(es.wa$maximum>es.w$maximum){
				es.w<-opsurv(x,y,'w',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.wa$estimate,'w',cni=c(F,T,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog,
							ub=c(.5,10,.5,10));}
		} else {es.wa$estimate<-rep.int(-Inf,3);}
		if(match('wb',models,nomatch=0)){
			es.wb<-opsurv(x,y,'w',cons=c(T,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.w$estimate,'w',cno=c(T,F,T,T),
						pf=pf,trim=T),tlog=tlog,ub=c(.5,10,.5,10));
			if(es.wb$maximum>es.w$maximum){
				es.w<-opsurv(x,y,'w',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.wb$estimate,'w',cni=c(T,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog,
							ub=c(.5,10,.5,10));}
		} else {es.wb$estimate<-rep.int(-Inf,3);}
		if(match('wab',models,nomatch=0)){
			es.wab<-opsurv(x,y,'w',cons=c(F,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.w$estimate,'w',cno=c(F,F,T,T),
						pf=pf,trim=T),tlog=tlog,ub=c(.5,10,.5,10));
			if(es.wab$maximum>es.w$maximum){
				es.w<-opsurv(x,y,'w',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.wab$estimate,'w',cni=c(F,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog,
							ub=c(.5,10,.5,10));}
		} else {es.wab$estimate<-rep.int(-Inf,2);}
	} else {es.w$estimate<-rep.int(-Inf,2*xrep);}
	if(any(c('g','ga','gb','gm','gma','gmb','gmc','l','la','lb','ls','lm',
		 'lma','lmb','lmc','lms','lmabcs')%in% models)){
          browser();
		es.g<-opsurv(x,y,'g',par=partemp,lb=c(1e-14,0,0,0),tlog=tlog,cx=cx,cy=cy, verbose=verbose);
		if(match('ga',models,nomatch=0)){
			es.ga<-opsurv(x,y,'g',cons=c(F,T,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.g$estimate,'g',cno=c(F,T,T,T),
						pf=pf,trim=T),tlog=tlog);
			if(es.ga$maximum>es.g$maximum){
				es.g<-opsurv(x,y,'g',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.ga$estimate,'g',cni=c(F,T,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.ga$estimate<-rep.int(-Inf,3);}; 
		if(match('gb',models,nomatch=0)){		
			es.gb<-opsurv(x,y,'g',cons=c(T,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.g$estimate,'g',cno=c(T,F,T,T),
						pf=pf,trim=T),tlog=tlog);
			if(es.gb$maximum>es.g$maximum){
				es.g<-opsurv(x,y,'g',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.gb$estimate,'g',cni=c(T,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.gb$estimate<-rep.int(-Inf,3);}; 
		partemp<-modpars(es.g$estimate,'g','gm',nil=nil,trim=T,onegrp=onegrp);
	} else {es.g$estimate<-rep.int(-Inf,2*xrep);}; 
	if(any(c('gm','gma','gmb','gmc','lm','lma','lmb','lmc','lms','lmabcs')%in% models)){
		es.gm<-opsurv(x,y,'gm',par=partemp,lb=c(1e-14,0,0,0),tlog=tlog,cx=cx,cy=cy, verbose=verbose);
		if(match('gma',models,nomatch=0)){		
			es.gma<-opsurv(x,y,'gm',cons=c(F,T,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.gm$estimate,'gm',cno=c(F,T,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.gma$maximum>es.gm$maximum){
				es.gm<-opsurv(x,y,'gm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.gma$estimate,'gm',cni=c(F,T,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.gma$estimate<-rep.int(-Inf,5);}
		if(match('gmb',models,nomatch=0)){		
			es.gmb<-opsurv(x,y,'gm',cons=c(T,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.gm$estimate,'gm',cno=c(T,F,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.gmb$maximum>es.gm$maximum){
				es.gm<-opsurv(x,y,'gm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.gmb$estimate,'gm',cni=c(T,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.gmb$estimate<-rep.int(-Inf,5);}
		if(match('gmc',models,nomatch=0)){		
			es.gmc<-opsurv(x,y,'gm',cons=c(T,T,F,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.gm$estimate,'gm',cno=c(T,T,F,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.gmc$maximum>es.gm$maximum){
				es.gm<-opsurv(x,y,'gm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.gmc$estimate,'gm',cni=c(T,T,F,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.gmc$estimate<-rep.int(-Inf,5);}
		if(any(c('lm','lma','lmb','lmc','lms','lmabcs') %in% models)){
			#print('112');browser();
			es.gmlm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.gm$estimate,'gm','lm',nil=nil,trim=T,
						    onegrp=onegrp),lb=c(1e-14,0,0,0),tlog=tlog);
		} else {es.gmlm$estimate<-rep.int(-Inf,4*xrep);}
	} else {es.gm$estimate<-rep.int(-Inf,3*xrep);};
	if(any(c('l','la','lb','ls','lm','lma','lmb','lmc','lms','lmabcs')%in% models)){
		es.l<-opsurv(x,y,'l',par=partemp,lb=c(1e-14,0,0,0),tlog=tlog,cx=cx,cy=cy, verbose=verbose);
		if(match('la',models,nomatch=0)){		
			es.la<-opsurv(x,y,'l',cons=c(F,T,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.l$estimate,'l',cno=c(F,T,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.la$maximum>es.l$maximum){
				es.l<-opsurv(x,y,'l',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.la$estimate,'l',cni=c(F,T,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.la$estimate<-rep.int(-Inf,5);}
		if(match('lb',models,nomatch=0)){
			es.lb<-opsurv(x,y,'l',cons=c(T,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.l$estimate,'l',cno=c(T,F,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.lb$maximum>es.l$maximum){
				es.l<-opsurv(x,y,'l',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lb$estimate,'l',cni=c(T,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lb$estimate<-rep.int(-Inf,5);}
		if(match('ls',models,nomatch=0)){
			es.ls<-opsurv(x,y,'l',cons=c(T,T,T,F),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.l$estimate,'l',cno=c(T,T,T,F),pf=pf,trim=T),
				tlog=tlog);
			if(es.ls$maximum>es.l$maximum){
				es.l<-opsurv(x,y,'l',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.ls$estimate,'l',cni=c(T,T,T,F),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.ls$estimate<-rep.int(-Inf,5);}
		if(any(c('lm','lma','lmb','lmc','lms','lmabcs') %in% models) ){
			es.llm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
				       par=modpars(es.l$estimate,'l','lm',nil=nil,trim=T,onegrp=onegrp),
				       lb=c(1e-14,0,0,0),tlog=tlog);
		} else {es.llm$estimate<-rep.int(-Inf,4*xrep);}
	} else {es.l$estimate<-rep.int(-Inf,3*xrep);}
	# Might wish to consider getting rid of es.gmlm and have the gm output go straight to es.lm
	if(es.gmlm$max<es.llm$max){es.lm<-es.llm;} else {es.lm<-es.gmlm;}
	if(is.finite(es.lm$maximum)){
		if(match('lma',models,nomatch=0)){
			es.lma<-opsurv(x,y,'lm',cons=c(F,T,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.lm$estimate,'lm',cno=c(F,T,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.lma$maximum>es.lm$maximum){
				es.lm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lma$estimate,'lm',cni=c(F,T,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lma$estimate<-rep.int(-Inf,7);}
		if(match('lmb',models,nomatch=0)){
			es.lmb<-opsurv(x,y,'lm',cons=c(T,F,T,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.lm$estimate,'lm',cno=c(T,F,T,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.lmb$maximum>es.lm$maximum){
				es.lm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lmb$estimate,'lm',cni=c(T,F,T,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lmb$estimate<-rep.int(-Inf,7);}
		if(match('lmc',models,nomatch=0)){
			es.lmc<-opsurv(x,y,'lm',cons=c(T,T,F,T),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.lm$estimate,'lm',cno=c(T,T,F,T),pf=pf,trim=T),
				tlog=tlog);
			if(es.lmc$maximum>es.lm$maximum){
				es.lm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lmc$estimate,'lm',cni=c(T,T,F,T),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lmc$estimate<-rep.int(-Inf,7);}
		if(match('lms',models,nomatch=0)){
			es.lms<-opsurv(x,y,'lm',cons=c(T,T,T,F),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.lm$estimate,'lm',cno=c(T,T,T,F),pf=pf,trim=T),
				tlog=tlog);
			if(es.lms$maximum>es.lm$maximum){
				es.lm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lms$estimate,'lm',cni=c(T,T,T,F),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lms$estimate<-rep.int(-Inf,7);}
		if(match('lmabcs',models,nomatch=0)){
			es.lmabcs<-opsurv(x,y,'lm',cons=c(F,F,F,F),cx=cx,cy=cy, verbose=verbose,
				par=modpars(es.lm$estimate,'lm',cno=c(F,F,F,F),pf=pf,trim=T),
				tlog=tlog);
			if(es.lmabcs$maximum>es.lm$maximum){
				es.lm<-opsurv(x,y,'lm',cx=cx,cy=cy, verbose=verbose,
					par=modpars(es.lms$estimate,'lm',cni=c(F,F,F,F),
							cno=c(T,T,T,T),pf=pf,trim=T),tlog=tlog);}
		} else {es.lmabcs$estimate<-rep.int(-Inf,4);}
	} else {es.lm$estimate<-rep.int(-Inf,4*xrep);};
	maxs<-c(es.e$max,es.w$max,es.wa$max,es.wb$max,es.wab$max,es.g$max,es.ga$max,es.gb$max,es.gm$max,
		es.gma$max,es.gmb$max,es.gmc$max,es.l$max,es.la$max,es.lb$max,es.ls$max,es.gmlm$max,
		es.llm$max,es.lm$max,es.lma$max,es.lmb$max,es.lmc$max,es.lms$max,es.lmabcs$max);
	iters<-c(0,es.w$titer,es.wa$titer,es.wb$titer,es.wab$titer,es.g$titer,es.ga$titer,es.gb$titer,
		es.gm$titer,es.gma$titer,es.gmb$titer,es.gmc$titer,es.l$titer,es.la$titer,
		es.lb$titer,es.ls$titer,es.gmlm$titer,es.llm$titer,es.lm$titer,es.lma$titer,es.lmb$titer,
		es.lmc$titer,es.lms$titer,es.lmabcs$titer);
	lrs<-c(0,es.w$max-es.e$max,es.w$max-es.wa$max,es.w$max-es.wb$max,es.w$max-es.wab$max,es.g$max-es.e$max,
		es.g$max-es.ga$max,es.g$max-es.gb$max,es.gm$max-es.g$max,es.gm$max-es.gma$max,
		es.gm$max-es.gmb$max,es.gm$max-es.gmc$max,es.l$max-es.g$max,es.l$max-es.la$max,
		es.l$max-es.lb$max,es.l$max-es.ls$max,es.gmlm$max-es.gm$max,es.llm$max-es.l$max,
		max(es.gmlm$max-es.gm$max,es.llm$max-es.l$max),
		es.lm$max-es.lma$max,es.lm$max-es.lmb$max,es.lm$max-es.lmc$max,es.lm$max-es.lms$max,
		es.lm$max-es.lmabcs$max);
	difs<-sign(lrs);
	pchi<-pchisq(2*lrs,c(1,1,1,1,3,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,3),low=F);
	issig<-pchi<sig;
	pars<-rbind(modpars(es.e$estimate,'e',nil=nil,trim=F),
			modpars(es.w$estimate,'w',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.wa$estimate,'w',nil=nil,cni=c(F,T,T,T),trim=F),
			modpars(es.wb$estimate,'w',nil=nil,cni=c(T,F,T,T),trim=F),
			modpars(es.wab$estimate,'w',nil=nil,cni=c(F,F,T,T),trim=F),
			modpars(es.g$estimate,'g',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.ga$estimate,'g',nil=nil,cni=c(F,T,T,T),trim=F),
			modpars(es.gb$estimate,'g',nil=nil,cni=c(T,F,T,T),trim=F),
			modpars(es.gm$estimate,'gm',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.gma$estimate,'gm',nil=nil,cni=c(F,T,T,T),trim=F),
			modpars(es.gmb$estimate,'gm',nil=nil,cni=c(T,F,T,T),trim=F),
			modpars(es.gmc$estimate,'gm',nil=nil,cni=c(T,T,F,T),trim=F),
			modpars(es.l$estimate,'l',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.la$estimate,'l',nil=nil,cni=c(F,T,T,T),trim=F),
			modpars(es.lb$estimate,'l',nil=nil,cni=c(T,F,T,T),trim=F),
			modpars(es.ls$estimate,'l',nil=nil,cni=c(T,T,T,F),trim=F),
			modpars(es.gmlm$estimate,'lm',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.llm$estimate,'lm',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.lm$estimate,'lm',nil=nil,trim=F,onegrp=onegrp),
			modpars(es.lma$estimate,'lm',nil=nil,cni=c(F,T,T,T),trim=F),
			modpars(es.lmb$estimate,'lm',nil=nil,cni=c(T,F,T,T),trim=F),
			modpars(es.lmc$estimate,'lm',nil=nil,cni=c(T,T,F,T),trim=F),
			modpars(es.lms$estimate,'lm',nil=nil,cni=c(T,T,T,F),trim=F),
			modpars(es.lmabcs$estimate,'lm',nil=nil,cni=c(F,F,F,F),trim=F)
		);
	npars<-c(2,4,3,3,2,4,3,3,6,5,5,5,6,5,5,5,8,8,8,7,7,7,7,4);
	models<-c('e','w','wa','wb','wab','g','ga','gb','gm','gma','gmb','gmc','l','la','lb','ls',
			'gmlm','llm','lm','lma','lmb','lmc','lms','lmabcs');
	aic<- 2*npars - 2*maxs;
	bic<- log(length(c(x,y)))*npars - 2*maxs; 
	out<-data.frame(maxs,iters,difs,pars,models,id,nil,lrs,aic,bic,pchi,issig);
	rownames(out)<-c('Exponential','Weibull','Weibull different a','Weibull different b','Weibull different any',
			 'Gompertz','Gompertz different a','Gompertz different b',
			 'Gompertz-Makeham','Gompertz-Makeham different a',
			 'Gompertz-Makeham different b','Gompertz-Makeham different c',
			 'Logistic','Logistic different a','Logistic different b',
			 'Logistic different s',
			 'Gompertz-Makeham -> Logistic-Makeham',
			 'Logistic -> Logistic-Makeham',
			 'Logistic-Makeham',
			 'Logistic-Makeham different a','Logistic-Makeham different b',
			 'Logistic-Makeham different c','Logistic-Makeham different s',
			 'Logistic-Makeham different any');
	colnames(out)<-c('MLE','# cycles','OK?','a1','b1','c1','s1',
			 'a2','b2','c2','s2','model','id','nil','LR','AIC','BIC','p (chi squared)','sig?');
	# The following line drops the unused models
	out<-(out[is.finite(out$MLE),]);
	if(!is.null(y)){
		unc.max0<-c(es.e$group0$max,es.w$group0$max,es.g$group0$max,es.gm$group0$max,
			    es.l$group0$max,es.lm$group0$max);
		unc.max1<-c(es.e$group1$max,es.w$group1$max,es.g$group1$max,es.gm$group1$max,
			    es.l$group1$max,es.lm$group1$max);
		unc.models<-c('e','w','g','gm','l','lm');
		unc.lrs0<-c(0,es.w$group0$max-es.e$group0$max,es.g$group0$max-es.e$group0$max,
			es.gm$group0$max-es.g$group0$max,es.l$group0$max-es.g$group0$max,
			es.lm$group0$max-max(es.gm$group0$max,es.l$group0$max));
		unc.lrs1<-c(0,es.w$group1$max-es.e$group1$max,es.g$group1$max-es.e$group1$max,
			es.gm$group1$max-es.g$group1$max,es.l$group1$max-es.g$group1$max,
			es.lm$group1$max-max(es.gm$group1$max,es.l$group1$max));
		unc.pchi0<-pchisq(2*unc.lrs0,df=1,low=F);
		unc.pchi1<-pchisq(2*unc.lrs1,df=1,low=F);
		out.unc<-data.frame(unc.max0,unc.max1,unc.models,unc.lrs0,unc.lrs1,unc.pchi0,
				    unc.pchi1,id,nil);	
	}

	if(!is.null(label)){
		mainname=paste(label,'.txt',sep='');
		if(file.exists(mainname)){
			write.table(out,file=mainname,sep='\t',append=T,col.names=F,row.names=F);
		} else {write.table(out,file=mainname,sep='\t',col.names=T,row.names=F);}
		if(writeunc&exists('out.unc')){
			uncname=paste(label,'_unc.txt',sep='');
			if(file.exists(uncname)){
				write.table(out.unc,file=uncname,sep='\t',append=T,
				col.names=F,row.names=F); } else {
				write.table(out.unc,file=uncname,sep='\t',col.names=T,row.names=F);}
		}
		if(writegml){
		gmlname=paste(label,'_gml.txt',sep='');
		if(file.exists(gmlname)){write.table(t(c(es.gmlm$max,es.gmlm$titer,es.gmlm$runtime,
				modpars(es.gmlm$estimate,'lm',trim=F),id,'gmlm',nil)),
			    file=gmlname,sep='\t',append=T,col.names=F,row.names=F);
		} else {write.table(t(c(es.gmlm$max,es.gmlm$titer,es.gmlm$runtime,
				modpars(es.gmlm$estimate,'lm',trim=F),id,'gmlm',nil)),
			    file=gmlname,sep='\t',col.names=T,row.names=F);}
		write.table(t(c(es.llm$max,es.llm$titer,es.llm$runtime,
				modpars(es.llm$estimate,'lm',trim=F),id,'llm',nil)),
			    file=gmlname,sep='\t',append=T,col.names=F,row.names=F);
		}
	}

	if(summary){
		cat('\n');
		print(out,digits=digits);#print(ests);
		print(sum(out[,2]));#invisible(list(results=out));
	}

	# somewhere here, test which results are relevant and keep only those
	invisible(out);
}
