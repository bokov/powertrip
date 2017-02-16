`gmean`<-
function(x){exp(sum(log(x))/length(x));}

`dropn`<-
function(x,n,which='smallest'){
	if(match('largest',which,nomatch=0)){
		for(i in 1:n){
			x<-x[-which.max(x)];
		}
	}
	if(match('largest',which,nomatch=0)){
		for(i in 1:n){
			x<-x[-which.min(x)];
		}
	}
}

`spyhaz`<-
function(label,path='./',extramain=0){
	#main<-paste(label,'main',sep='.');
	main<-label;
	gml<-paste(label,'gml',sep='.');
	unc<-paste(label,'unc',sep='.');
	gmlcols<-maincols<-c("max", "iters", "rt", "a1", "b1", "c1", "s1", "a2", "b2", "c2",
			      "s2", "id", "fit", "nilval");
	unccols<-c('max1','max2','id','fit','nilval');
	length(maincols)<-length(maincols)+extramain;
	assign(main,read.table(paste(path,main,'.txt',sep=''),col.names=maincols),envir=.GlobalEnv);
	assign(gml,read.table(paste(path,gml,'.txt',sep=''),col.names=gmlcols),envir=.GlobalEnv);
	assign(unc,read.table(paste(path,unc,'.txt',sep=''),col.names=unccols),envir=.GlobalEnv);
	cat('Output collected from',max(get(gml)$id),'iterations.\n');
}


`templatelr`<-
function(perm,perm2,which,...){
	switch(which,
	lm={plotmultlrs(perm,perm2,c('lm','lm','lm','lm'),c('lma','lmb','lmc','lms'),
		...)},
	unc={plotmultlrs(perm,perm2,c('gm','l','lm','lm','lm'),c('g','g','l','gm','g'),
		...)}
	)
}

`rowcols`<-
function(n){
	if(n<5){return(c(1,n));}
	if(n==7|n==8){return(c(2,4));}
	if(12<n&n<=15){return(c(3,5));}
	out<-sqrt(n); rout<-round(out);
	if(out<=rout){return(c(rout,rout))
	}else{return(c(rout,rout+1))};
}

`plotmultestdens`<-
function(perms,pars1=c('a1','b1','c1','s1'),pars2=c('a2','b2','c2','s2'),cex=.75,...){
	len<-length(pars1);
	# maybe some input validation for pars1 and pars2 length later
	close.screen(all=T);
	par(cex=cex);par(bg='white');scrs<-split.screen(rowcols(len));
	for(i in 1:len){
		screen(scrs[i]);
		plotestdens(perms,,pars1[i],pars2[i],...);
	}
	invisible(scrs);
}

`plotestdens`<-
function(perms,fits=c('w','wa','wb','g','ga','gb','gm','gma','gmb','gmc','l','la','lb','ls',
		      'lm','lma','lmb','lmc','lms'),
	 par1='a1',par2='a2',log=F,ymax=0.02){
	if(length(grep('c',c(par1,par2)))>0){fits<-fits[grep('m',fits)];} else {
		if(length(grep('s',c(par1,par2)))>0){fits<-fits[grep('l',fits)];}
	}
	temp<-list();overall<-c(); lperms<-length(perms); lfits<-length(fits);
	num<-lperms*lfits;
	cols1<-matrix(rainbow(num,start=.5,end=1),nrow=lperms,ncol=lfits);
	cols2<-matrix(rainbow(num,start=.5,end=1,s=.75,v=1),nrow=lperms,ncol=lfits);
	for(i in perms){
		iperm<-get(i);
		for(j in fits){
			jperm1<-subset(iperm,fit==j)[[par1]];
			jperm2<-subset(iperm,fit==j)[[par2]];
			if(log){jperm1<-log(jperm1);jperm2<-log(jperm2);}
			if(sum(!is.na(jperm1))>0){
				temp[[i]][[j]][[1]]<-density(jperm1,na.rm=T);
				temp[[i]][[j]][[1]]$y<-
					temp[[i]][[j]][[1]]$y/sum(temp[[i]][[j]][[1]]$y);
			}
			if(sum(!is.na(jperm1))>0){
				temp[[i]][[j]][[2]]<-density(jperm2,na.rm=T);
				temp[[i]][[j]][[2]]$y<-
					temp[[i]][[j]][[2]]$y/sum(temp[[i]][[j]][[2]]$y);
			}
			overall<-c(overall,jperm1,jperm2);
		}
	}
	xlim<-range(overall,finite=T);
	plot(temp[[1]][[1]][[1]],xlim=xlim,ylim=c(0,ymax),
	     main=paste('Distribution of Estimates \nfor Parameters',par1,'&',par2),
	     sub='');
	for(i in 1:lperms){
		for(j in 1:lfits){
			lines(temp[[i]][[j]][[1]],col=cols1[i,j]);
			lines(temp[[i]][[j]][[2]],col=cols2[i,j]);
		}
	}
}

`plotmultlrs`<-
function(perm,perm2,fits1=c('g','g','gm','gm','gm','l','l','l','lm','lm','lm','lm'),
	 fits2=c('ga','gb','gma','gmb','gmc','la','lb','ls','lma','lmb','lmc','lms'),
	 df=1,breaks=100,abcol='red',cex=.75,extraab=NULL,...){
	close.screen(all=T);
	len<-min(length(fits1),length(fits2));breaks<-rep(breaks,len=len);
	abcol<-rep(abcol,len=len); df<-rep(df,len=len);
	length(fits1)<-length(fits2)<-len;
	par(cex=cex);par(bg='white');scrs<-split.screen(rowcols(len));
	for(i in 1:len){
		screen(scrs[i]); #par(mar=rep(3,4));
		plotlr(perm,perm2,fits1[i],fits2[i],breaks=breaks[i],abcol=abcol[i],
		       df=df[i],extraab=extraab[i],...);
	}
	invisible(scrs);
}

`plot.mainhaz`<-
function(perm,breaks=100,add=F,col=1,abcol='red',onescreen=T,what='cons'){
	if(what=='cons'){
		if(onescreen){layout(matrix(1:12,4,3,byrow=T));}
		plotlr(perm,'g','ga',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'g','gb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gma',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gmb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','gmc',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','la',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','lb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','ls',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lma',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lmb',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lmc',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','lms',breaks=breaks,col=col,abcol=abcol,add=add);
	}
	if(what=='model'){
		if(onescreen){layout(matrix(1:6,3,2,byrow=T));}
		plotlr(perm,'g','e',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'gm','g',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'l','g',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','gm',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','l',breaks=breaks,col=col,abcol=abcol,add=add);
		plotlr(perm,'lm','g',breaks=breaks,col=col,abcol=abcol,add=add,df=2);
	}
	if(onescreen){layout(matrix(1,1,1));}
}

`survpowab`<-
function(x,label,model='ga',reps=6,nil=0,bnil=0,ns=NULL,add=0,mult=1.1,lrcol='LR',sig=.05,pow=.8,parmax=1){
	# the if is.null idiom is to keep the function declaration uncluttered
	# if(is.null(ns)){ns<-c(5,10,15,20,25,30,35,40,45,50,100,200,500,1000);}
	if(is.null(ns)){ns<-c(20,50,100,200,500);}
	chisig<-qchisq(sig,df=1,lower=F);
	lx<-length(x); basemodel<-sub('a|b|c|s','',model); testparam<-sub('w|g|gm|l|lm','',model);
	realfit<-findpars(x,nil=nil,bnil=bnil,models=basemodel);
	# the hardcoded g is not compatible with 'c' and 's' parameters
	nullpars<-realfit[realfit[['model']]==basemodel,c('a1','b1','c1','s1')];
	names(nullpars)<-c('a','b','c','s'); 
	nullpars[is.na(nullpars)]<-nil;
	if(nullpars[testparam]==0){mult<-1;if(add==0){add<-0.001;}}
	for(i in ns){
		cat('\nn:',i,' ');
		# simsurv automatically downgrades to 'g' family of models if s ~ 0
		ix<-simsurv(i,'l',nullpars);
		pars<-nullpars; plr<-0; pars[testparam]<-pars[testparam]*mult+add; out<-c();
		while(plr<pow & pars[testparam]<=parmax & pars[testparam]>0){
			cat('plr:',plr,' ');
			lr<-c();
			for(j in 1:reps){
				redo=T;
				while(redo){
					jx<-simsurv(i,'l',pars);
					jfit<-try(findpars(ix,jx,nil=nil,bnil=bnil,models=model));
					if(class(jfit)[1]!="try-error"){redo=F;}else{print(jx);}
				}
				jlr<-jfit[jfit[['model']]==model,lrcol]; lr<-c(lr,jlr);
				jfit[['rep']]<-j; jfit[['par']]<-c('a','b','c','s')[testparam];
				jfit[['N']]<-i; jfit[['mean']]<-mean(jx); jfit[['sd']]<-sd(jx);
				jfit[['q50']]<-quantile(jx,.5);jfit[['q75']]<-quantile(jx,.75);
				jfit[['q90']]<-quantile(jx,.9);jfit[['q95']]<-quantile(jx,.95);
				jfit<-cbind(jfit,pars);jfit[['pow']]<-NA;
				jfit[['pLR']]<-NA;jfit[['LR']]<-jlr;
				out<-rbind(out,jfit);
			}
			pars[testparam]<-pars[testparam]*mult+add;
			plr<-sum(lr>chisig)/length(lr);print(lr);print(pars);print(add);
			out[is.na(out$pow),]$pow<-plr;
			out[is.na(out$pLR),]$pLR<-quantile(lr,1-pow);
		}
		if(file.exists(label)){
			write.table(out,file=label,sep='\t',col.names=F,row.names=F,append=T);
		} else {write.table(out,file=label,sep='\t',col.names=T,row.names=F,append=T);}
	}
}

`plotlr`<-
function(x,realid=1,model='g',lrcol=4,modcol=3,idcol=8,k=2,breaks=100,barcol='black',linecol='red',draw=T){
	# x is the dataframe obtained by reading the output file from findpars
	# k is two unless the log ratio has already been doubled in which case
	# is should be set to 1. realid is the id corresponding to the real, 
	# non-bootstrapped comparison. Set it to null in order to not draw the
	# line for that comparison
	# most of the other parameters are the corresponding columns where to
	# look for the respective values
	if(draw){hist(k*x[x[,modcol]==model,lrcol],freq=F,border=barcol,breaks=breaks);}
	out<-c();
	if(!is.null(realid)){
		reallr<-k*x[x[,idcol]==1&x[,modcol]==model,lrcol];
		if(draw){abline(v=reallr,col=linecol);}
		myecdf<-ecdf(k*x[x[,modcol]==model,lrcol]);
		out<-1-myecdf(reallr);
	}
	return(out);
}

`plotlr.old`<-
function(perm,perm2,fit1='gm',fit2='g',breaks=100,df=1,col=c(1,'green'),
	 abcol='red',extraab=NULL,extracol='blue',add=F,log=F,ylim=0.2,...){
	lrs<-2*(subset(perm,fit==fit1)$max-subset(perm,fit==fit2)$max);
	empp<-1-ecdf(lrs)(lrs[1]);
	chip<-pchisq(lrs[1],df,lower=F);
	badlrs<-lrs<0; sumbad<-sum(badlrs);
	if(exists('perm2')){
		lrs2<-2*(subset(perm2,fit==fit1)$max-subset(perm2,fit==fit2)$max);
		badlrs2<-lrs<0; sumbad<-sum(badlrs2);
		if(log){lrs2<-log(lrs2);}
	}
	if(log){lrs<-log(lrs);} 
	if(exists('perm2')){xlim<-range(lrs,lrs2,finite=T);
	} else {xlim<-range(lrs,finite=T);}
	breaks=seq(xlim[1],xlim[2],len=breaks); 

	#cat(fit1,'vs',fit2,'fraction <0:',sumbad/length(lrs),'\n');
	#print(summary(lrs[badlrs]));
	#if(sumbad<10){cat(lrs[badlrs],'\n');} else {cat(lrs[badlrs][1:10],'...\n');}
	bars<-hist(lrs,breaks=breaks,plot=F); 
	bars$counts<-bars$counts/sum(bars$counts);
	plot(bars,col=NULL, border=col[1],main=paste(fit1,'vs',fit2),
		 sub='',add=add,ylim=c(0,ylim),xlim=xlim,
		 xlab=paste('Emp. p:',format(empp,digits=3),
			    '  Chi-sq p:',format(chip,digits=3)));
	abline(v=lrs[1],col=abcol);
	if(exists('perm2')){
		bgbak<-par('bg');par(bg='transparent');
		bars<-hist(lrs2,breaks=breaks,plot=F);
		bars$counts<-bars$counts/sum(bars$counts);
		plot(bars,col=NULL,border=col[2],main='',sub='',xlab='',add=T,
		     ylim=c(0,ylim),xlim=xlim,lty=4);
		abline(v=lrs2[1],col=abcol);
		par(bg=bgbak);
	}
	if(!is.null(extraab)){abline(v=extraab,col=extracol);}
}


`plotest`<-
function(perm,val='gm',par='a1',breaks=100,add=F,col=1){
	ests<-subset(perm,fit==val)[[par]];
	hist(ests,breaks=breaks,add=add,main=paste(par,'parameter for the',val,'fit'),border=col);
	abline(v=ests[1],col='red');
}

`parnames`<-
function(i,j){
	model<-constraint<-parameters<-c();
	i<-toupper(i);
	switch(i,
		{model<-'LM'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','c','s','a','b','c','s');},
				{constraint<-'a'; parameters<-c('a','b','c','s','b','c','s');},
				{constraint<-'b'; parameters<-c('a','b','c','s','a','c','s');},
				{constraint<-'c'; parameters<-c('a','b','c','s','a','b','s');},
				{constraint<-'s'; parameters<-c('a','b','c','s','a','b','c');}
			);
		},
		{model<-'L'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','s','a','b','s');},
				{constraint<-'a'; parameters<-c('a','b','s','b','s');},
				{constraint<-'b'; parameters<-c('a','b','s','a','s');},
				{constraint<-'s'; parameters<-c('a','b','s','a','b');}
			);
		},
		{model<-'GM'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','c','a','b','c');},
				{constraint<-'a'; parameters<-c('a','b','c','b','c');},
				{constraint<-'b'; parameters<-c('a','b','c','a','c');},
				{constraint<-'c'; parameters<-c('a','b','c','a','b');}
			);
		},
		{model<-'G'; switch(j,
				{constraint<-'u'; parameters<-c('a','b','a','b');},
				{constraint<-'a'; parameters<-c('a','b','b');},
				{constraint<-'b'; parameters<-c('a','b','a');}
			);
		}
	)
	list(model=model,constraint=constraint,parameters=parameters);
}


`print.spyfit`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(1:12,4,3,byrow=T));
		for(i in 4:1){
				switch(i,{model='Logistic-Makeham';pars=1:4;},
					 {model='Logistic';pars=c(1,2,4);},
					 {model='Gompertz-Makeham';pars=1:3;},
					 {model='Gompertz';pars=1:2;},);
				for(j in 2:length(x[[i]])){
					param<-switch(pars[j-1],'a','b','c','s');
					ijecdf<-ecdf(x[[i]][[j]]);
					ijp<-1-ijecdf(x[[i]][[j]][1]);
					ijxlab<-paste('parameter',param);
					ijsub<-paste('p =',ijp);
					hist(x[[i]][[j]],breaks=100,xlab=ijxlab,ylab='',
					     freq=T,main=model,sub=ijsub,...);
					abline(v=x[[i]][[j]][1],col='red');
					abline(v=attr(x[[i]][[j]],'ci')[1],col='green',lty=2);
					hist(attr(x[[i]][[j]],'ci'),breaks=100,freq=T,border='red',
					     lty=2,add=T);
				}
			}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		rnames<-out<-c();
		for(i in 1:length(x)){
			for(j in 2:length(x[[i]])){
				out<-rbind(out,c(quantile(x[[i]][[j]],.05),median(x[[i]][[j]]),
						 x[[i]][[j]][1],mean(x[[i]][[j]]), 
						 quantile(x[[i]][[j]],.95),
						 1-ecdf(x[[i]][[j]])(x[[i]][[j]][1]),
						 pchisq(x[[i]][[j]][i],1,lower=F)));
				out<-rbind(out,c(quantile(attr(x[[i]][[j]],'ci'),.05),
						 median(attr(x[[i]][[j]],'ci')),
						 attr(x[[i]][[j]],'ci')[1],mean(attr(x[[i]][[j]],'ci')),
						 quantile(attr(x[[i]][[j]],'ci'),.95),NA,NA));
				ijlabel<-paste(toupper(names(x))[i],',',names(x[[i]])[j],'constrained');
				rnames<-c(rnames,paste('Permuted',ijlabel),paste('Resampled',ijlabel));
			}
		}
		colnames(out)<-c('5%','median','actual','mean','95%','emp p','chisq p');
		rownames(out)<-rnames;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}

`print.spyinput`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(1:12,4,3,byrow=T));
		for(i in 4:1){
				switch(i,{model='Logistic-Makeham';pars=1:4;},
					 {model='Gompertz-Makeham';pars=1:3;},
					 {model='Logistic';pars=c(1,2,4);},
					 {model='Gompertz';pars=1:2;},);
				for(j in 1:length(pars)){
					param<-switch(pars[j],'a','b','c','s');
					ijxlab<-paste('parameter',param);
					ijx0<-x[[i]]$u[,3+j];
					ijx1<-x[[i]]$u[,3+j+length(pars)];
					hist(ijx0[ijx0>0],breaks=100,xlab=ijxlab,freq=F,
					     ylab='',main=model,...);
					hist(ijx1[ijx1>0],breaks=100,lty=2,add=T,freq=F);
					hist(attr(x[[i]]$u,'ci')[,3+j],breaks=100,border='red',freq=F,
					     lty=2,add=T);
					hist(attr(x[[i]]$u,'ci')[,3+j+length(pars)],breaks=100,
					     border='red',freq=F,lty=3,add=T);
				}
			}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		rnames<-out<-c();
		for(i in 1:length(x)){
			for(j in 1:length(x[[i]])){
				ijnames<-parnames(i,j); 
				ijnames$parameters<-c('MLE','iter','rt',ijnames$parameters);
				for(k in 1:length(x[[i]][[j]])){
					out<-rbind(out,c(quantile(x[[i]][[j]][[k]],.05),
						   median(x[[i]][[j]][[k]]),x[[i]][[j]][[k]][1],
						   mean(x[[i]][[j]][[k]]),
						   quantile(x[[i]][[j]][[k]],.95)));
					out<-rbind(out,c(quantile(attr(x[[i]][[j]],'ci')[[k]],.05),
						   median(attr(x[[i]][[j]],'ci')[[k]]),
						   attr(x[[i]][[j]],'ci')[[k]][1],
						   mean(attr(x[[i]][[j]],'ci')[[k]]),
						   quantile(attr(x[[i]][[j]],'ci')[[k]],.95)));
					ijkname<-paste(ijnames$model,ijnames$constraint,
						       ijnames$parameters[k]);
					rnames<-c(rnames,paste(c('Permuted','Resampled'),ijkname))
				}
			}
		}
		colnames(out)<-c('5%','median','actual','mean','95%');
		rownames(out)<-rnames;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}

`print.spymodels`<-
function(x,type=c('plot','table'),...){
	if(match('plot',type,nomatch=0)){
		layout(matrix(c(1:5,0),2,3,byrow=T));
		for(i in 1:5){
			model<-switch(i,'Gompertz vs. Logistic',
					'Gompertz vs. Gompertz-Makeham',
					'Logistic vs. Logistic-Makeham',
					'Gompertz-Makeham vs. Logistic-Makeham',
					'Gompertz vs. Logistic-Makeham');
			iecdf<-ecdf(x[[i]]);
			ip<-1-iecdf(x[[i]][1]);
			isub<-paste('p =', ip);
			hist(x[[i]],breaks=100,xlab='',ylab='',main=model,freq=T,sub=isub,...)
			abline(v=x[[i]][1],col='red');
			abline(v=attr(x[[i]],'ci')[1],col='green',lty=2);
			hist(attr(x[[i]],'ci'),breaks=100,freq=T,border='red',lty=2,add=T);
		}
		layout(matrix(1,1,1));
	}
	if(match('table',type,nomatch=0)){
		out<-c();
		compars<-c('G vs. L', 'G vs. GM', 'L vs. LM', 'GM vs. LM', 'G vs. LM');
		compars<-expand.grid(c('Permuted','Resampled'),compars);
		compars<-paste(compars[,1],compars[,2]);
		for(i in 1:length(x)){
			if(i!=5){idf<-1;}else{idf<-2;}
			out<-rbind(out,c(quantile(x[[i]],.05),median(x[[i]]),x[[i]][1],mean(x[[i]]),
					 quantile(x[[i]],.95),1-ecdf(x[[i]])(x[[i]][1]),
					 pchisq(x[[i]][i],idf,lower=F)));
			out<-rbind(out,c(quantile(attr(x[[i]],'ci'),.05),median(attr(x[[i]],'ci')),
					 attr(x[[i]],'ci')[1],mean(attr(x[[i]],'ci')),
					 quantile(attr(x[[i]],'ci'),.95),NA,NA));
		}
		colnames(out)<-c('5%','median','actual','mean','95%','emp p','chisq p');
		rownames(out)<-compars;
		oldwidth<-options('width')[[1]];options('width'=130);
		print(out);options('width'=oldwidth);invisible(out);
	}
}


`raw2wm`<-
function(t){
	if(!is.list(t)){t<-list(t);}
	out<-c();
	for(i in 1:length(t)){
		irle<-rle(sort(t[[i]]));
		iwm<-cbind(irle$values,irle$lengths,i,1);
		out<-rbind(out,iwm);
	}
	out<-rbind(out,c(-1,NA,NA,NA));
	colnames(out)<-NULL;
	out;
}

`raw2tab`<-
function(x,y,labels=c(0,1)){cbind(c(x,y),c(rep(labels[1],length(x)),rep(labels[2],length(y))));}

`int2tab`<-
function(x,y,labels=c(0L,1L))
	{cbind(as.integer(c(x,y)),c(rep(labels[1],length(x)),rep(labels[2],length(y))));}

# this is to pull out just the significant model and parameters from the output of findpars
`hztblsumm`<-
function(x,fcols=c('MLE','a1','b1','c1','s1','a2','b2','c2','s2','model','LR','p (chi squared)'),
	 frows='g$|gm$|l$|lm$', modcol=fcols[10], pcol=fcols[12], cutoff=.05){
	# fcols are the current names of the columns of interest, frows is a grep
	# search term for pulling out the currently implemented models (the $ avoids
	# pulling out the constrained versions of those models). modcol and pcol
	# identify the name of the input columns containing the model name and the p value,
	# respectively
	out<-x[grep(frows,x[,modcol]),fcols];
	issig=T; i=2; n.models<-dim(out)[1];
	# check each unconstrained model to see if it passes the cutoff. If it doesn't, stop.
	# NOTE: Need to hardcode an exception for the LM model, such that only the one with the
	# lower p-value is included in the table by the time this calculation is run
	while(issig & i < n.models){if(out[i,pcol]>cutoff){issig<-F;
						}else{bestmdl<-out[i,modcol];i<-i+1;}}
	# bestmdlc is assigned the grep search term for all the possible constrained versions
	bestmdlc<-paste(bestmdl,'[a|b|c|s]',sep='');
	# now we pull out the constrained fits from the input that correspond to the
	# best model 
	out<-rbind(x[grep('w',x[,modcol]),fcols],out,x[grep(bestmdlc,x[,modcol]),fcols]);
	# now we add and populate a 'sig' column for quick visual reference of which models
	# are good fits; note that for the constraints the fit criteria are REVERSED and
	# might not be the right way to do it.
	out$sig<-''; out[match(bestmdl,out[,modcol]),'sig']<-'*';
	out[out[,pcol]>.05,'sig']<-'*'; attr(out,'n.models')<-n.models;
	return(out);
}

exdens<-function(t,a){exp(-t/a)/a}
wdens<-function(t,a,b){exp((b-1)*log(t)+(a^b)*(-t^b)+log(b)+log(a)*b)}
gdens<-function(t,a,b){a*exp(b*t-a*(exp(b*t)-1)/b)}
gmdens<-function(t,a,b,c){(c+a*exp(b*t))*exp(-a*(exp(b*t)-1)/b -c*t)}
ldens<-function(t,a,b,s){a*exp(b*t)*(1+s*a*(exp(b*t)-1)/b)^(-(s+1)/s)}
lmdens<-function(t,a,b,c,s){(c+a*exp(b*t)/(1+s*a*(exp(b*t)-1)/b))*(exp(-c*t)*(1+s*a*(exp(b*t)-1)/b)^(-(s+1)/s))}