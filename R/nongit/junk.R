`startpars.old`<-
function(x,y,label,perms=10){
	d<-int2tab(x,y); lxy<-dim(d)[1]; pr<-matrix(nrow=lxy,ncol=perms);
	for(i in 1:perms){pr[,i]<-sample(d[,2],lxy);}
	fits<-list(lm=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			   b=list(cons=c(1L,0L,1L,1L),out=c()), c=list(cons=c(1L,1L,0L,1L),out=c()),
			   s=list(cons=c(1L,1L,1L,0L),out=c())),
		   l=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			  b=list(cons=c(1L,0L,1L,1L),out=c()), s=list(cons=c(1L,1L,1L,0L),out=c())),
		   gm=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			   b=list(cons=c(1L,0L,1L,1L),out=c()), c=list(cons=c(1L,1L,0L,1L),out=c())),
		   g=list(u=list(cons=c(1L,1L,1L,1L),out=c()), a=list(cons=c(0L,1L,1L,1L),out=c()),
			  b=list(cons=c(1L,0L,1L,1L),out=c())));	

	cat('\nCalculating unconstrained fits from default starting parameters.\n');
	for(i in 1:perms){
		for(j in 1:length(fits)){
			jout<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],names(fits)[j]);
			fits[[j]]$u$out<-rbind(fits[[j]]$u$out,
					       c(jout$maximum,jout$titer,jout$runtime,jout$estimate));
		}
	}

	cat('\nCalculating unconstrained fits from candidate optimal parameters.\n');
	for(j in 1:length(fits)){
		#browser();
		jcols<-4:dim(fits[[j]]$u$out)[2];
		jin<-fits[[j]]$u$out[-1,jcols];
		jmaxs<-numeric(4); jpars<-list();
		jmaxs[1]<-fits[[j]]$u$out[1,1];

		jpars[[1]]<-fits[[j]]$u$out[1,jcols];
		jpars[[2]]<-apply(jin,2,mean);
		jpars[[3]]<-apply(jin,2,median);
		jpars[[4]]<-apply(jin,2,gmean);

		for(l in 2:4){
			jmaxs[l]<-opsurv(x,y,names(fits)[j],par=jpars[[l]])$maximum;
		}
		fits[[j]]$u$label<-paste(label,'_',names(fits)[j],'u',sep='');
		jmaxmax<-which.max(jmaxs); if(jmaxmax>1){
			cat('\nBetter pars found for',fits[[j]]$u$label,':',jmaxmax,'\n');
		}
		fits[[j]]$u$par<-jpars[[jmaxmax]];
		#cat('\nmaxima:',jmaxs,'\n');
		write.table(cbind(jmaxs,t(as.data.frame(jpars))),
			    file=paste(fits[[j]]$u$label,'_opars.txt',sep=''),
			    col.names=F,row.names=F,sep='\t');
	}
	
	cat('\nCalculating constrained fits from unconstrained optimal parameters.\n');
	for(i in 1:perms){
		for(j in 1:length(fits)){
			for(k in 2:length(fits[[j]])){
				jkout<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],
					      names(fits)[j],
					      cons=fits[[j]][[k]]$cons,
					      par=fits[[j]]$u$par);
				fits[[j]][[k]]$out<-rbind(fits[[j]][[k]]$out,
							  c(jkout$maximum,jkout$titer,jkout$runtime,	
							    jkout$estimate));
			}
		}
	}

	cat('\nCalculating constrained from candidate optimal parameters.\n');
	for(j in 1:length(fits)){
		for(k in 2:length(fits[[j]])){
			jkcols<-4:dim(fits[[j]][[k]]$out)[2];
			jkin<-fits[[j]][[k]]$out[-1,jkcols];
			jkmaxs<-numeric(4); jkpars<-list();
			jkmaxs[1]<-fits[[j]][[k]]$out[1,1];
			jkpars[[1]]<-fits[[j]]$u$out[1,jkcols];
			jkpars[[2]]<-apply(jkin,2,mean);
			jkpars[[3]]<-apply(jkin,2,median);
			jkpars[[4]]<-apply(jkin,2,gmean);
	
			for(l in 2:4){
				jkmaxs[l]<-opsurv(x,y,names(fits)[j],cons=fits[[j]][[k]]$cons,
				                  par=jkpars[[l]])$maximum;
			}
			fits[[j]][[k]]$label<-paste(label,'_',names(fits)[j],names(fits[[j]])[k],sep='');
			jkmaxmax<-which.max(jkmaxs); if(jkmaxmax>1){
				cat('\nBetter pars found for',fits[[j]][[k]]$label,':',jkmaxmax,'\n');
			}
			fits[[j]][[k]]$par<-jkpars[[jkmaxmax]];
			write.table(cbind(jkmaxs,t(as.data.frame(jkpars))),
				    file=paste(fits[[j]][[k]]$label,'_opars.txt', sep=''), 
				    col.names=F,row.names=F,sep='\t');
		}
	}
	outname<-paste(label,'pars',sep='');
	assign(outname,fits,envir=.GlobalEnv);
	save(list=outname,file=paste(outname,'.rdata',sep=''));
	invisible(fits);
}

`ressurv.old` <-
function(x,y,label,resam=1000,pars=NULL,...){
	lx<-length(x);ly<-length(y);
	rsx<-matrix(nrow=lx,ncol=resam);
	rsy<-matrix(nrow=ly,ncol=resam);
	for(i in 1:resam){
		rsx[,i]<-sample(x,lx,replace=T);
		rsy[,i]<-sample(y,ly,replace=T);
	}
	rsx<-cbind(x,rsx); rsy<-cbind(y,rsy);
	if(is.character(pars)){parslabel<-pars;pars<-NULL;}
	else{parslabel<-paste(label,'pars',sep='');}
	# need to write some kind of validation function for fits
	if(is.null(pars)|length(pars)!=4){
		if(exists(parslabel)){pars<-get(parslabel);}
		else{pars<-startpars.old(x,y,label,perms=10);}
	}
	for(i in 1:(resam+1)){
		for(j in names(pars)){
			for(k in names(pars[[j]])){
				out<-opsurv(rsx[,i],rsy[,i],model=j,cons=pars[[j]][[k]]$cons,
					    par=pars[[j]][[k]]$par);
				outtab<-c(out$maximum,out$titer,out$runtime,out$estimate,
					  out$group0$estimate,out$group1$estimate);
				write.table(t(outtab),paste(label,'_resam_',j,k,'.txt',sep=''),sep='\t',
					    append=T,col.names=F,row.names=F);
			}
		}
	}
}

`permsurv.old` <-
function(x,y,label,perms=5000,pars=NULL){
	d<-int2tab(x,y); lxy<-dim(d)[1]; pr<-matrix(nrow=lxy,ncol=perms);
	for(i in 1:perms){pr[,i]<-sample(d[,2],lxy);}
	pr<-cbind(d[,2],pr);
	if(is.character(pars)){parslabel<-pars;pars<-NULL;}
	else{parslabel<-paste(label,'pars',sep='');cat('\nLooking for',parslabel,'\n');}
	# need to write some kind of validation function for fits
	if(is.null(pars)|length(pars)!=4){
		if(exists(parslabel)){cat('\nFound',parslabel,'\n');pars<-get(parslabel);}
		else{pars<-startpars.old(x,y,label,perms=10);}
	}
	cat('\nCommencing permutations.\n');
	for(i in 1:(perms+1)){
		for(j in 1:length(pars)){
			for(k in 1:length(pars[[j]])){
				jkts<-as.numeric(Sys.time());
				jklabel<-paste(pars[[j]][[k]]$label,jkts,sep='');
				jkfit<-opsurv(d[pr[,i]==0,1],d[pr[,i]==1,1],names(pars)[j],
						cons=pars[[j]][[k]]$cons, par=pars[[j]][[k]]$par);
				jktab<-c(jkfit$maximum,jkfit$titer,jkfit$runtime,jkfit$estimate,
						jkfit$group0$maximum,jkfit$group1$maximum);
				assign(jklabel,jkfit);
				dump(jklabel, paste(pars[[j]][[k]]$label,'.R',sep=''),
				     append=T,control='S_compatible');
				#print('');print(jktab);
				write.table(t(jktab),paste(pars[[j]][[k]]$label,'.txt',sep=''),sep='\t',
					    append=T,col.names=F,row.names=F);
			}
		}
	}
	save(list=ls(),file=paste(label,'debug.rdata',sep=''));
}


`hazwrapper.old` <-
function(data,foo,outlog='hazlog.R',debuglog='hazdebug.R',debugvars=c('estimate','maximum','titers','runtime'),...){
	ts<-as.numeric(Sys.time()); data[,1]<-data[,1][foo]; 
	outname<-paste('out',ts,sep='.');
	assign(outname,opsurv('lm',data[data[,2]==0,1],data[data[,2]==1,1],...));
	if(!is.null(outlog)){
		dump(outname,outlog,append=T,control='S_compatible')
	}
	get(outname)$maximum;
}

`spyhaz.old`<-
function(label,show=c('models','lm','l','gm','g','pars')){
	fit<-input<-list(lm=list(u=list(),a=list(),b=list(),c=list(),s=list()), 
			 l=list(u=list(),a=list(),b=list(),s=list()),
			 gm=list(u=list(),a=list(),b=list(),c=list()),
			 g=list(u=list(),a=list(),b=list()));
	models<-list(g.l=c(),g.gm=c(),l.lm=c(),gm.lm=c());
	class(fit)<-'spyfit';class(input)<-'spyinput';class(models)<-'spymodels';
	for(i in 1:length(input) ){
		for(j in 1:length(input[[i]]) ){
			input[[i]][[j]]<-read.table(paste(label,'_',names(input)[i],names(input[[i]])[j],
							  '.txt',sep=''));
			ijcifile<-paste(label,'_resam_',names(input)[i],names(input[[i]])[j],
					'.txt',sep='');
			if(file.exists(ijcifile)){
				ijci<-read.table(ijcifile);
				#attr(input[[i]][[j]],'ci')<-apply(ijci,2,quantile,c(.05,.95));
				attr(input[[i]][[j]],'ci')<-ijci;
			} else ijci<-NULL;
			if(j>1){ 
				fit[[i]][[j]]<-2*(input[[i]]$u[,1] - input[[i]][[j]][,1]);
				if(!is.null(ijci) & exists('iuci')){
# 					attr(fit[[i]][[j]],'ci')<-quantile(2*(iuci[,1]-ijci[,1]),
# 									   c(.05,.95));
					attr(fit[[i]][[j]],'ci')<-2*(iuci[,1]-ijci[,1]);
				}
			} else { iuci<-ijci; }
		}
	}
	models$g.l<-2*(input$l$u[,1] - input$g$u[,1]);
	attr(models$g.l,'ci')<-2*(attr(input$l$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);
	models$g.gm<-2*(input$gm$u[,1] - input$g$u[,1]);
	attr(models$g.gm,'ci')<-2*(attr(input$gm$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);
	models$l.lm<-2*(input$lm$u[,1] - input$l$u[,1]);
	attr(models$l.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$l$u,'ci')[,1]);
	models$gm.lm<-2*(input$lm$u[,1] - input$gm$u[,1]);
	attr(models$gm.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$gm$u,'ci')[,1]);
	models$g.lm<-2*(input$lm$u[,1] - input$g$u[,1]);
	attr(models$g.lm,'ci')<-2*(attr(input$lm$u,'ci')[,1]-attr(input$g$u,'ci')[,1]);

	assign(paste(label,'fit',sep=''),fit,envir=.GlobalEnv);
	assign(paste(label,'input',sep=''),input,envir=.GlobalEnv);
	assign(paste(label,'models',sep=''),models,envir=.GlobalEnv);
}

`spyhaz2`<-
function(logfile,realval=NULL){
	env<-new.env();sys.source(logfile,env);logtemp<-as.list(env);
	out<-list(); out$log<-logtemp;
	out$loglen<-loglen<-length(logtemp);
	if(is.null(realval)){out$realval<-realval<-logtemp[[loglen]];}
	out$realmax<-realmax<-realval$maximum;
	ests<-maxs<-iters<-rts<-c();
	for(i in logtemp){
		ests<-rbind(ests,i$estimate);
		maxs<-c(maxs,i$maximum);iters<-c(iters,i$titer);rts<-c(rts,i$runtime);
	};
	out$ests<-ests; out$maxs<-maxs; out$iters<-iters; out$rts<-rts;
	lests<-dim(ests)[2];
	hist(maxs,breaks=1000,xlim=c(min(maxs),max(maxs)),freq=F);
	abline(v=realmax,col='red');myecdf<-ecdf(maxs);
	out$permp<-1-myecdf(realmax); out$avgrt<-sum(rts)/loglen; out$loglen<-loglen;
	maxests<-minests<-means<-gmeans<-medians<-numeric(lests);
	for(i in 1:lests){
		means[i]<-mean(ests[,i]);medians[i]<-median(ests[,i]);gmeans[i]<-gmean(ests[,i]);
		maxests[i]<-max(ests[,i]);minests[i]<-min(ests[,i]);
	}
	avpars<-rbind(means,medians,gmeans,maxests,minests);
	out$avpars<-avpars;
	class(out)<-'spyhaz';
	print(out);
	invisible(out);
}

`print.spyhaz2`<-
function(x,...){
	cat('permuted p:',x$permp,'\n');
	cat('avg runtime:',x$avgrt,'\n');
	cat('# perms:',x$loglen,'\n');
	print(x$avpars);
}

# This is the old version, which calls the interpreted versions of the functions and will be removed once the new compiled version is verified
`opsurv2` <-
function(model,x,y=NULL,par=rep(1e-6,4),cons=1,old=F,usegr=T,usehs=F,debug=F,
	 lb=c(1e-15,1e-10,1e-10,1e-8),ub=c(.5,.5,.5,3),called=NULL,method='nm',...){
	if(!exists("rm.debug")){rm.debug<<-list(0);} 
	call<-sys.call(); tstart<-proc.time()[3];
	# below needs to have the function names replaced with chars: "ofg","ofgm","ofl","oflm"
	switch(model,
		g={np<-2;ex<-exg;hs<-hsg;keep<-c(1,2,5,6);},
		gm={np<-3;ex<-exgm;hs<-hsgm;keep<-c(1:3,5:7);},
		l={np<-3;ex<-exl;hs<-hsl;keep<-c(1,2,4,5,6,8);},
		lm={np<-4;ex<-exlm;hs<-hslm;keep<-1:8;},
		stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
	);
	rm.temp<<-list(cons=checkcons(cons,np,keep),cons2=rep(cons,le=4),x=x,y=y,
			keep=keep[1:np],old=old,debug=debug,called=called,np=np);
	par<-fixpar(par,np,keep);
	if(sum(is.na(par)>0)){
	stop("At least one of the starting parameters (par) you 
specified is an NA or NaN. Here are the parameters:\n",
paste(par,collapse=' '),
"\nPlease specify a different set of parameters or just omit them 
and use the defaults.");}
	lb<-fixpar(lb,np,keep);ub<-fixpar(ub,np,keep); rm.temp$lb<<-lb;rm.temp$ub<<-ub;
	fn<-function(par){
		objf2(par,lb=lb,ub=ub,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],ex=ex);
	}
	gr<-function(par){grf2(par,cons=rep(cons,le=4),x=x,y=y,np=np,keep=keep[1:np],model=model);}
	if(!usegr) gr<-NULL; if(!usehs) hs<-NULL;
	change<-1; out<-list(par); totaliter<-0;
        ctrl$parscale<-rep.int(1, length(par)); ctrl$ndeps<-rep.int(0.001, length(par));
	while(max(change)>0){
		out.old<-out;
		cat(".");
		options(show.error.messages=F);
		out<-try(.Internal(optim(out.old[[1]],fn,gr,"Nelder-Mead",ctrl,lb,ub)))
		options(show.error.messages=T);
		if(class(out)[1]!="try-error"){
			totaliter<-out[[3]][1]+totaliter;
			change<-abs(out.old[[1]]-out[[1]]);
		} else {
			p<-out.old[[1]]; lk<-NA;
			# Bad starting params, huh? We'll see about that.
			while(is.na(lk)|is.infinite(lk)){
			 p<-p/2; lbp<-p>lb;
			 p<-p*lbp+1.1*lb*(1-lbp);
			 if(max(p>1.1*lb)==0){cat('X\n');print(p);break();}
			 lk<-objf2(p,lb,ub,rep(cons,le=4),
			          max(x),max(y),np=np,keep=keep[1:np],ex=ex);
			 # Damn starting params! You have bested me, I surrender.
			 cat('!');
			}
			if(!is.na(lk)){out<-list(p);change<-1;} else {break;}
		} 
	}
	cat(model," ");
	if(class(out)[1]!="try-error"){
		names(out)<-c('estimate','maximum','iterations','code','message');
		if(usegr) { out$gradient<-gr(out$estimate); }
		out$call<-call;
		out$titer<-totaliter; 
		out$runtime<-proc.time()[3]-tstart;
	}
	if(debug){cat("Runtime:",proc.time()[3],"\n");}
	#out$totaliter<-totaliter;
	out;
}

# This is the old, interpreted version, which will be removed once the compiled version is verified
`objf2` <-
function(par,lb=0,ub=1,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,ex=exg){
	if(max(par<=lb)!=0|max(par>=ub)!=0){return(NA);}
	par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np];
	par2[keep][cons[keep]>0]<-par[(np+1):length(par)];
	env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
	env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
	out<-sum(eval(ex,envir=env1))+sum(eval(ex,envir=env2));
	if(!is.na(out) & is.finite(out) & out < -50){return(out);}else{return(NA);}
}

# simulate populations using models
`simsurv.old`<-
function(n,type='g',p=c(2.433083e-05,0.005,3e-11,0.0015),maxt=1e5){
	if(length(p)!=4){stop('You must specify four parameters.
Fill them out with zeroes if necessary.');}
	tt<-1:maxt;out<-integer(n);
	switch(type,
		g={sf<-function(x,a,b,c,s=0){exp(-a*(exp(b*x)-1)/b -c*x);}},
		l={sf<-function(x,a,b,c,s){exp(-c*x)*(1+s*a*(exp(b*x)-1)/b)^(-1/s);}},
		w={sf<-function(x,a,b,c=0,s=0){exp(-(a*x)^b);}});
	risk<-1-sf(tt+1,p[1],p[2],p[3],p[4])/sf(tt,p[1],p[2],p[3],p[4]);
	risk<-risk[!is.na(risk)]
	for(i in 1:n){
		out[i]<-match(T,runif(length(risk))<risk);
	}
	attr(out,'type')<-type;attr(out,'pars')<-p;
	return(out);
}

`scramble` <-
function(x,y,n,srt=T){
	lx<-length(x);ly<-length(y);
	xy <- c(x,y); lxy<-length(xy);
	outx<-c(); outy<-c();
	for(i in 1:n){
		ixy<-sample(xy,lxy)
		ix<-ixy[1:lx]; iy<-ixy[(lx+1):(lx+ly)];
		if(srt){ix<-sort(ix); iy<-sort(iy);}
		outx<-cbind(outx,ix); outy<-cbind(outy,iy);
	}
	return(list(x=outx,y=outy));
}

# some scratch functions that are not actually used
`permstat` <-
function(x,y,n,fn,srt=F,statname='statistic',...){
	fnout<-fn(x,y,...);
	perms<-scramble(x,y,n,srt);
	sout<-list(); length(sout)<-length(statname);
	names(sout)<-statname;
	for(i in 1:n){
		ifnout<-fn(perms$x[,i],perms$y[,i],...);
		for(j in statname){ sout[[j]]<-c(sout[[j]],ifnout[[j]]);}
	}
	for(j in statname){
		jecdf<-ecdf(sout[[j]]);
		cat(j,':',fnout[[j]],'permuted p:',1-jecdf(fnout[[j]]),'\n');
		cat('empirical quantiles:',quantile(sout[[j]],c(.05,.1,.5,.9,.95)),'\n');
	}
	invisible(list(f=fnout,perms=sout));
}

`twrapper` <-
function(data,foo,...){
	foox<-data[,1][foo[data[,2]==0]]; fooy=data[,1][foo[data[,2]==1]];
	out<-t.test(x=foox,y=fooy,alternative='two.sided',mu=0);
	out$statistic;
}




`compsurvs`<-
function(a,b){
	out1<-cbind(c(a$maximum,a$titer),c(b$maximum,b$titer));
	rownames(out1)<-c('maximum','titer');
	out2<-rbind(a$estimate,b$estimate);
	print(out1,digits=22);print(out2);
	invisible(list(out1,out2));
}

