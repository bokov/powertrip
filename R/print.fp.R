`print.fp`<-function(d,what=c('fitmodels','sigparams'),verbose=1,showsigonly=T,digits=3, columns=c('Model','modlo','MLE','df','LR','a','b','c','s','p','sig'),print.gap=1,filename=NULL){
    what<-match.arg(what,several.ok=T); out<-list();
    outputwidth=80; # happens to be the
    if(!is.null(filename)) write.csv('',file=filename);
    if('fitmodels' %in% what){
	cat('\nCandidate Model Fits\n');
	cat(paste(rep('-',outputwidth),collapse=''),'\n');
	mf<-function(a,b,name,columns){
	    a<-data.frame(a);a$Model<-rownames(a);
	    a<-merge(b,a,by='Model',sort=F);
	    levels(a$Model)<-modelinfo(levels(a$Model),'full');
	    levels(a$modlo)<-modelinfo(levels(a$modlo),'full');
	    
	    a<-a[,columns];
	    columns<-gsub('modlo','Compared to',gsub('LR','Log ratio',columns));
	    names(a)<-columns;
	    cat('\n',name,'\n');
	    print(a,row.names=F,digits=digits,print.gap=print.gap);
	    return(a);
	}
	
	out$xmodels<-mf(d$unc[,,'x'],d$lrsx,d$groupnames[1],columns);
	out$ymodels<-mf(d$unc[,,'y'],d$lrsy,d$groupnames[2],columns);
	out$modsig<-d$modsig; print(d$modsig,print.gap=print.gap);
        if(!is.null(filename)){
          write(d$groupnames[1],file=filename,append=T);
          write.csv(out$xmodels,file=filename,append=T);
          write('',file=filename,append=T);
          write(d$groupnames[2],file=filename,append=T);
          write.csv(out$ymodels,file=filename,append=T);
        }
    }
    out$jointmodels<-d$jointmodels; print(d$jointmodels,print.gap=print.gap);
    if(!is.null(filename)) write(out$jointmodels,file=filename,append=T);
    if('sigparams' %in% what){
	cat('\nParameter differences between',d$groupnames[1],'and',d$groupnames[2],'\n');
	cat(paste(rep('-',outputwidth),collapse=''),'\n');
	#out<-data.frame(Model=NA,'Parameter/s'=NA,MLE=NA,'Log ratio'=NA,df=NA,'p (chi squared)'=NA,'p (empirical'=NA);
	out$paramdiff<-data.frame(NA); n<-1;
	for(i in names(d$cons)) for(j in 2:nrow(d$cons[[i]])){
	    ijcons<-1+d$cons[[i]][j,];
	    out$paramdiff[n,'Model']<-i;
	    out$paramdiff[n,'Parameter/s']<-paste(colnames(d$cons[[i]])[!as.logical(d$cons[[i]][j,])],collapse=',');
	    out$paramdiff[n,'MLE']<-d$pers[i,ijcons[1],ijcons[2],ijcons[3],ijcons[4],'MLE',1];
	    out$paramdiff[n,'Log ratio']<-d$pers[i,ijcons[1],ijcons[2],ijcons[3],ijcons[4],'LR',1];
	    out$paramdiff[n,'df']<-d$pers[i,ijcons[1],ijcons[2],ijcons[3],ijcons[4],'df',1];
	    ijecdf<-ecdf(d$pers[i,ijcons[1],ijcons[2],ijcons[3],ijcons[4],'LR',][-1]);
	    out$paramdiff[n,'p (chi squared)']<-pchisq(2*out$paramdiff[n,'Log ratio'],out$paramdiff[n,'df'],lower=F);
	    out$paramdiff[n,'p (empirical)']<-1-ijecdf(out$paramdiff[n,'Log ratio']);
	    n<-n+1; out$paramdiff<-rbind(out$paramdiff,NA);
	}
	print(out$paramdiff[-n,-1],row.names=F,digits=digits,print.gap=print.gap);
        if(!is.null(filename)){
          write('',file=filename,append=T);
          write('Parameter differences',file=filename,append=T);
          write.csv(out$paramdiff[-n,-1],file=filename,append=T);
        }
    }
    class(out)<-c('fp.summary','list');
    invisible(out);
}
