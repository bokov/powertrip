`constrshow`<-
function(d,x=NULL,k=2,datacol='LR',idcol='id',idobs=1,modelcol='model',nullmodelcol='null_model',
	 chicol='p (chi squared)',figure=T,models=c('wa','wb','ga','gb'),breaks=200,hcol='black',
	 vcol='red',chidf=1,rc=NULL){
	if(is.null(x)){x<-d[d[[idcol]]==idobs,];}; if(is.null(rc)){rc<-rowcols(length(models));}
	if(match('emp.p',names(x),nomatch=0)==0){x<-cbind(x,'emp.p'=NA);updatex=T;}
	close.screen(all=T); scrs<-split.screen(rc);
	for(i in 1:length(models)){
		screen(i);
		reallr<-x[x[[modelcol]]==models[i],datacol];
		if(updatex){
			iecdf<-ecdf(k*d[d[[modelcol]]==models[i],datacol]);
			x[x[[modelcol]]==models[i],'emp.p']<-1-iecdf(reallr);
		}
		if(figure){
			hist(k*d[d[[modelcol]]==models[i],datacol],breaks=breaks,freq=F,border=hcol,
			     main=paste('MLE Log Ratio:',models[i]),xlab='');
			title(sub=sprintf('p chisq = %1.4G,   p emp = %1.4G', x[x[[modelcol]]==models[i],chicol],x[x[[modelcol]]==models[i],'emp.p']));
			if(!is.null(chidf)){
				lines(seq(0,15,l=500),dchisq(seq(0,15,l=500),df=chidf),col=vcol);
			}
			abline(v=reallr,col=vcol);
		}
	}
	invisible(x);	
}