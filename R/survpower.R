survpower<-function(x,n=50,cx=NULL,xpars=NULL,pow=0.8,sig=0.05,afac=0.5,bfac=0.5,
                    abounds=NULL,bbounds=NULL,tests=c(50,100,500),
                    #fun=function(z,x) surv2.logrank(Surv(c(z,x),c(sign(z),cx)),rep(1:2,c(length(z),length(x))))$pval,
                    #fun=function(z,x) surv2.logrank(Surv(c(z,x)),rep(1:2,c(length(z),length(x))))$pval,
		    fun=function(z,x) summary(coxph(Surv(c(z,x))~rep(1:2,c(length(z),length(x)))))$sctest['pvalue'],
                    weibull=F,xmodel=NULL,tempname='survpower'){

#n<-sort(n,T); tests<-sort(tests);
if(is.null(cx)) cx<-sign(x);
# if parameters aren't known for the sample population estimate them
if(is.null(xpars)){
  xpars<-findpars(x);
  # if model has not been chosen, choose one
  if(is.null(xmodel)){
    if(!weibull) xpars<-subset(xpars,model!='w');
    xmodel<-choosemodel(xpars[xpars[,'sig?'],'model']);
  }
  xpars<-as.numeric(xpars[xpars$model==xmodel,c('a1','b1','c1','s1')]);
}

out<-list();
# for each n
for(ni in n){
  if(afac > 0){
    mypars<-xpars;
    if(is.null(abounds)){
      abounds<-c(xpars[1],0);
      mypars[1]<-mypars[1]*afac;
    } else mypars[1]<-abounds[2];
    aout<-findbounds(ni,sample(cx,n,rep=T),xpars,mypars,xmodel,pow,sig,
                   matrix(c(afac,abounds,NA,NA,NA),nrow=3),tests,NA,NA,fun);
  # if tempname save aout to a temp file
  }

  if(bfac > 0){
    mypars<-xpars;
    if(is.null(bbounds)){
      bbounds<-c(xpars[2],0);
      mypars[2]<-mypars[2]*bfac;
    } else mypars[2]<-bbounds[2];
    bout<-findbounds(ni,sample(cx,n,rep=T),xpars,mypars,xmodel,pow,sig,
                     matrix(c(NA,NA,NA,bfac,bbounds),nrow=3),tests,NA,NA,fun);
  # if tempname save aout to a temp file
  }

  if(afac > 0 & bfac > 0){
  # similar to a and b but now take the final values for a and b as starting values
    mypars<-c(aout$mypars[1],bout$mypars[2],xpars[3:4]);
    about<-findbounds(ni,sample(cx,n,rep=T),xpars,mypars,xmodel,pow,sig,
                      matrix(c(afac,xpars[1],0,bfac,xpars[2],0),nrow=3),tests,NA,NA,fun);
  # these will probably be a detectable difference, so we step them both down by averaging
  # with original xpars or up by multiplying by their respective factors
  }
  out$ni<-list(aout=aout,bout=bout,about=about);
}

return(out);
# return, for each n...
# call
# the upper and lower bounds on a, b, a&b
# summary of sample means for a, b, a&b
# summary of sample sds for a, b, a&b
# summary of sample medians for a, b, a&b
# summary of sample 90th percentiles for a, b, a&b
}


