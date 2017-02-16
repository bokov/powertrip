`modelshow`<-
function(label,x=NULL,reload=F,k=2,datacol='LR',modelcol='model',nullmodelcol='null_model',
	 chicol='p (chi squared)',figure=T, models=c('w','g','gm','l','lm'),breaks=200,hcol='black',vcol='red',
	 chidf=1,rc=NULL,namesep='_',fpsuffix='.fp',header=1,sep='\t'){
	# x is output from findpars; the names of data.frames from which to make histograms
	# are constructed from x and then looked for in the parent environment
	if(is.null(x)){
		input<-paste(label,fpsuffix,sep=''); if(!reload & exists(input)){x<-get(input);} else {
			rinput<-paste(label,'.rdata',sep='');
			if(file.exists(rinput)){
				finput<-load(rinput);x<-get(finput);assign(input,x,envir=.GlobalEnv);
			}
		}
	}
	ms<-c();
	for(m in models){
		nm<-.models[[m]]$nests;
		for(nmi in nm){				# sep used to be '', changed it to '_'
			if(length(nm)>1){pairname<-paste(nmi,m,sep='');}else{pairname<-m;}
			ms<-rbind(ms,data.frame(nmod=nmi,mod=m,pairname=pairname,
				  dname=paste(label,m,sep=namesep),stringsAsFactors=F));
		}
	}
	ml<-dim(ms)[1];updatex=F;
	if(!is.null(x)&match('emp.p',names(x),nomatch=0)==0){x<-cbind(x,'emp.p'=NA);updatex=T;}
	if(is.null(rc)){rc<-rowcols(ml);}
	close.screen(all=T); scrs<-split.screen(rc);
	for(i in 1:ml){
		screen(i);
		fname<-ms[i,'dname'];
		if(!reload & exists(fname)){d<-get(fname);} else {
			if(file.exists(fname)){
				d<-read.table(fname,header=header,sep=sep); assign(fname,d,envir=.GlobalEnv);
			} else rm(d);
		}
		if(exists('d')){
			if(figure){
				hist(k*d[d[[modelcol]]==ms[i,'pairname'] & 
						d[[nullmodelcol]]==ms[i,'nmod'],datacol],
					breaks=breaks,freq=F,border=hcol,
				main=paste('MLE Log Ratio:',ms[i,'nmod'],'vs',ms[i,'mod']),xlab='');
			}
			if(is.null(x)){x<-d[d$id==0,];}
			if(!is.null(x)){
				reallr<-x[x[[modelcol]]==ms[i,'pairname'],datacol];
				if(figure){abline(v=reallr,col=vcol);}
				if(updatex){
					iecdf<-ecdf(k*d[d[[modelcol]]==ms[i,'pairname'] & 
						    d[[nullmodelcol]]==ms[i,'nmod'],datacol]);
					x[x[[modelcol]]==ms[i,'pairname'],'emp.p']<-1-iecdf(reallr);
				}
			}
			if(!is.null(chidf)&figure){
				lines(seq(0,15,l=500),dchisq(seq(0,15,l=500),df=chidf),col=vcol);
			}
			if(figure){
				title(sub=sprintf('p chisq = %1.4G,   p emp = %1.4G',
					x[x[[modelcol]]==ms[i,'pairname'],chicol],x[x[[modelcol]]==ms[i,'pairname'],'emp.p']));
			}
		} else if(figure){hist(0,main=paste('Cannot find data.frame',fname));}
	}
	close.screen(all=T);
	if(!exists(input)){assign(input,x,envir=.GlobalEnv)}
	invisible(x);
}

