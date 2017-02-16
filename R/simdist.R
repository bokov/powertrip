`simdist`<-
function(x,label,cx=NULL,nil=0,bnil=0,wbnil=1,pf=mean,rounds=5000,models=c('w','g','gm','l','lm'),
	 dropcols=c('a2','b2','c2','s2'),sims=NULL,pars=NULL){
	# if sims and pars non-null, should write some kind of validation to make sure they match
	if(is.null(pars)){pars<-list();}
	lx<-length(x);
	ihaz<-findpars(x,cx=cx,nil=nil,bnil=bnil,wbnil=wbnil,models=models); 
	if(is.null(sims)) {
		sims<-list();
		msize<-lx*rounds;
		for(m in models){
			#.models is a data object with information about the supported models
			# that gets loaded when the library is loaded
			#nm<-.models[[m]]$nests;
			for(nmi in modelinfo(m)){
				pairname=paste(nmi,m,sep='');
				npars<-ihaz[ihaz$model==nmi,c('a1','b1','c1','s1')];
				npars[is.na(npars)]<-0;
				pars[[pairname]]<-npars;
				if(match(nmi,names(sims),nomatch=0)==0){
					sims[[nmi]]<-matrix(simsurv(msize,nmi,as.numeric(npars)),nrow=lx,ncol=rounds);
					sims[[nmi]]<-cbind(x,sims[[nmi]]);
					cat('Done simulating',nmi,'\n');
				}
			}
		}
		save(list=c('pars','sims'),file=paste(label,'.sims.rdata',sep=''));
	} else {rounds<-min(unlist(lapply(sims,function(x){dim(x)[2]})));}
	save(ihaz,file=paste(label,'.rdata',sep=''));
	print('Commencing parameter fitting');
	# used to be 1:rounds, but that may cause errors if sims argument is not null
	for(i in 1:(rounds+1)){
		for(m in models){
			#nm<-.models[[m]]$nests;
			for(nmi in modelinfo(m)){
				pairname=paste(nmi,m,sep='');
				#if(m=='w'){cat('\n');print(sims[[nmi]][,i]);browser();}
				out<-findpars(sims[[nmi]][,i],cx=cx,nil=nil,bnil=bnil,wbnil=wbnil,
					      models=m,id=i-1);
				out<-out[,setdiff(names(out),dropcols)];
				out<-cbind(out,null_pars=pars[[pairname]],null_model=nmi,target_model=m);
				mlabel<-paste(label,m,sep='_');
				if(file.exists(mlabel)){cn=FALSE;} else {cn=TRUE;}
				write.table(out,file=mlabel,sep='\t',col.names=cn,row.names=F,append=T);
			}
		}
	}
}
