mh<-function(par0,par1,l0,l1){
	if(is.na(l1)){return(par0);}
	if(l1>l0){return(par1);}
	return(ifelse(rbinom(1,1,exp(l1-l0)),par1,par0));
}

sampler<-function(x,y=NULL,pars,cx=1,cy=1,cycles=10000,sds=c(.1,.1,.1,.1),model='g',abs=F){
	nx<-length(x); ny<-length(y); 
	if(length(cx)==1){cx<-rep(1,nx);}
	if(length(cy)==1){cy<-rep(1,ny);}
	lp<-length(pars); ex<-paste('of',model,sep='');
	out<-matrix(NA,nrow=cycles,ncol=lp);
	switch(model,
		w={np<-2;keep<-c(1,2,5,6);ub<-c(.5,10,0,0);},
		g={np<-2;keep<-c(1,2,5,6);ub<-10;},
		gm={np<-3;keep<-c(1:3,5:7);ub<-10;},
		l={np<-3;keep<-c(1,2,4,5,6,8);ub<-10;},
		lm={np<-4;keep<-1:8;ub<-10;},
		stop("The model you specified, ",model,", does not exist.
The only models supported in this version are: 
g, gm, l, and lm.")
	);
	sds<-c(sds,sds)[keep];
	if(is.null(y)){sds<-sds[1:(length(sds)/2)];}
	ipars<-pars; bar<-txtProgressBar(1,cycles,initial=1,title='Sampling...',style=3);
	for(i in 1:cycles){
		out[i,]<-ipars;
		# is taking an absolute value problematic?
		if(abs){randpars<-abs(pars+rnorm(lp,sd=sds));}else{randpars<-pars+rnorm(lp,sd=sds);}
		for(j in 1:lp){
			jpars<-ipars; jpars[j]<-randpars[j];
			l0<-objf(ipars,x=x,y=y,nx=nx,ny=ny,cx=cx,cy=cy,ub=ub,keep=keep,np=np,ex=ex);
			l1<-objf(jpars,x=x,y=y,nx=nx,ny=ny,cx=cx,cy=cy,ub=ub,keep=keep,np=np,ex=ex);
			#if(is.na(l0)|is.na(l1)){browser();}
			ipars[j]<-mh(ipars[j],randpars[j],l0,l1);
		}
		setTxtProgressBar(bar,i);
	}
	close(bar);
	return(out);
}