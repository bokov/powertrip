smartpars<-function(x,y,cx=NULL,cy=NULL,tlog=F,digits=22,sig=0.05,tidy=T){
	xp<-findpars(x,cx=cx);yp<-findpars(y,cy=cy);
	xm<-as.character(tidyfits(xp,rm.rejected.models=T)$model);
	ym<-as.character(tidyfits(yp,rm.rejected.models=T)$model);
	ms<-c(xm,ym);
	# find the most complicated justified joint model for x and y
	#browser();
	models<-c();
	if(any(grep('^w$',ms))){models<-c(models,'wa','wb')}
	if(any(grep('^lm$',ms))|(any(grep('^gm$',ms))&any(grep('^l$',ms)))){
		models<-c(models,'lma','lmb','lmc','lms');} else {
		if(any(grep('^gm$',ms))){models<-c(models,'gma','gmb','gmc');} else {
			if(any(grep('^l$',ms))){models<-c(models,'la','lb','ls');}else{
				if(any(grep('^g$',ms))){models<-c(models,'ga','gb');}
			}
		}
	}
	out<-findpars(x,y,cx=cx,cy=cy,tlog=tlog,digits=digits,sig=sig,models=models);
	out<-out[,grep('MLE|a1|b1|c1|s1|a2|b2|c2|s2|model|LR|AIC|BIC|p .chi squared.|sig?',names(out))]
	cat('\n');
	return(tidyfits(out,rm.rejected.models=tidy));
}