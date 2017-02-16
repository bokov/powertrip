# mledrv<-function(x,y=NULL,model='g',pars,what='hess'){
# 	mydrv<-get(paste('d',model,sep=''));
# 	pars<-as.numeric(pars);
# 	nullout<-rep(NA,length(pars));
# 	if(!is.null(y)){return(nullout);}
# 	out<-mydrv(pars,x);
# 	switch(what,
# 		hess = return(tryCatch(apply(attributes(out)$hessian,c(2,3),sum),
# 			      error=function(e){print(e);return(nullout);})),
# 		sd = return(tryCatch(sqrt(abs(diag(solve(apply(attributes(out)$hessian,c(2,3),sum))))),
# 			      error=function(e){print(e);return(nullout);})),
# 		grad = return(tryCatch(apply(attributes(out)$gradient,2,sum),
# 			      error=function(e){print(e);return(nullout);}))
# 	)
# }

# when one of the parameters is zero, gives NA answers; need to come up with a clean way to catch that
# if a simpler model gets substituted, sd and grad should get run through modpars as follows:
# modpars(sd,model,origmodel,nil=0,onegrp=T)
# no idea what, if anything, to do with hess
mledrv<-function(x,y=NULL,model='g',pars,what='hess'){
	pars<-as.numeric(pars); origmodel<-model;
	nullout<-rep(NA,length(pars));
	if(any(pars==0)){
		if(model=='lm'){if(pars[4]==0){pars<-pars[-4];model<-'gm';} else {
			if(pars[3]==0){pars<-pars[-3];model<-'l';}}}
		if(model=='l'|model=='gm'){if(pars[3]==0){pars<-pars[-3];model<-'g';}}
		if((model=='g'|model=='w')&any(pars==0)){return(nullout);}
	}
	hz<-get(paste('d',model,'h',sep=''));
	sv<-get(paste('d',model,'s',sep=''));
	if(!is.null(y)){return(nullout);}
	hout<-hz(pars,x);sout<-sv(pars,x);
	switch(what,
		hess = return(tryCatch(apply(attributes(hout)$hessian+attributes(sout)$hessian,c(2,3),sum),
			      error=function(e){print(e);return(nullout);})),
		sd = return(tryCatch(sqrt(abs(diag(solve(apply(attributes(hout)$hessian+attributes(sout)$hessian,
							 c(2,3),sum))))),
			      error=function(e){print(e);return(nullout);})),
		grad = return(tryCatch(apply(attributes(hout)$gradient+attributes(sout)$gradient,2,sum),
			      error=function(e){print(e);return(nullout);}))
	)
}