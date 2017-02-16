`go` <-
function(x,y,xynames=c(),save=T,path,prompt=2,xlim=1350,slwd=4,
scol=c('darkred',1),spch=24:25,spcex=0,splwd=slwd,spbg=scol,sleg=T,sxcex=1.5,sxlwd=3,qstart=1,demint=1, qrse='boot',dohz=1){
	library(tcltk);
	options(warn=-1);options(show.error.messages=T);
	quantiles<-list(seq(.01,.99,.01),seq(.05,.95,.05),seq(.1,.9,.1),seq(.2,1,.2));
	mycall<-strsplit(as.character(sys.call()[1]),"\\$")[[1]];
	if(length(mycall)>1){
		top=F;d<-get(mycall[1]);
		if(length(mycall)>2){
			for(i in 2:(length(mycall)-1)){
				d<-d[[mycall[i]]];
		}}
		x<-d$x;y<-d$y;xynames<-d$xynames;smry<-d$smry;zsc<-d$zsc;xy<-d$xy;
		group<-d$group;lr<-d$lr;qreg<-d$qreg;qreg.sum<-d$qreg.sum;tt<-d$tt;
		qreg.tab<-d$qreg.tab;mod<-d$mod;demint<-d$demint;
		xd<-d$xd;yd<-d$yd;xysurvfit<-d$xysurvfit;
		sig.tests<-d$sig.tests;tests<-d$tests; report<-d$report; path<-d$path;
		sigqreg<-d$sigqreg;sigzsc<-d$sigzsc;slwd<-d$slwd;
		scol<-d$scol;spch<-d$spch;spcex<-d$spcex;splwd<-d$splwd;spbg<-d$spbg;
		sleg<-d$sleg;sxcex<-d$sxcex;sxlwd<-d$sxlwd; dohz<-d$dohz;
	}else{top=T;}
	if(missing(path)){
		if(top&save){
			message<-"Where would you like to save your results?\n";
			choices<-c("Let me choose using the point and click method.",
				"Let me type in the location from the console.",
				"I don\'t want to save any files. Run analysis without saving.");
			cat(message);
			input<-menu(choices,graphics=T,title=message);
			switch(input,
				{path<-paste(tclvalue(tkchooseDirectory()),"/",sep="");},
				{path<-readline("Please type in the location where you would like to save files: ");
				path<-paste(path,"/",sep="");},
				{save<-F;path=NULL;cat('\nContinuing without saving.');}
			);
			if(save){cat('\nResults will be saved to:',path,'\n');}
		} else {path=NULL;}
	}
	if(missing(x)){
		cat("\nYou need a group to compare. Please copy and paste one column of survival times from a spreadsheet.
They don\'t have to be in order, but make sure there are no empty spaces. Hit the enter key when you 
are done.\n");
		x<-scan(); 
	}
	if(length(xynames)<2){xynames[1]<-readline("Please give a brief name for the first group: ");}
	if(missing(y)){
		cat("\nYou need a second group to compare. Please copy and paste one column of survival times from a 
spreadsheet. They don\'t have to be in order, but make sure there are no empty spaces. Hit 
the enter key when you are done.\n");
		y<-scan(); 
	}
	if(length(xynames)<2){xynames[2]<-readline("Please give a brief name for the second group: ");}
	xynames<-gsub('[[:space:]]+','.',xynames)
	if(save&top){
		assign(xynames[1],x,envir=.GlobalEnv);
		assign(xynames[2],y,envir=.GlobalEnv);
		cat("The survival times for the two groups have been saved as",xynames[1],"and",xynames[2],"inside R 
so you won\'t have to paste them in next time. Instead you can type 
\'go(",xynames[1],",",xynames[2],")\'\n");
	}

	if(top){
		tests<-c("t test"="tt","z score"="zsc","log rank"="lr","quantile regression"="qreg.tab");
		if(dohz){tests<-c(tests,"mortality parameters"="mod");}
		sig.tests<-c();
	}

	cat("\nDoing summary statistics...");
	if(top){smry<-t(cbind(summary(x),summary(y))); rownames(smry)<-xynames;}
	cat('\n'); print(smry);
	cat("\nDoing t-test on log-transformed data...");
	if(top){xt<-x[which(x>0)]; yt<-y[which(y>0)];tt<-t.test(log(xt),log(yt));if(tt$p.value<=.05) sig.tests<-c(sig.tests,1);}
	cat('\n'); print(tt);
	cat("\nDoing z-score test to see which quantiles differ..."); 
	if(top){
		zsc<-c();
		zsc<-ezz(x,y,xynames,quant=seq(0,.95,.05));
		if(length(levels(zsc$sig))>1) sig.tests<-c(sig.tests,2);
	} 
	cat('\n'); print(zsc);
	cat("\nDoing log-rank test with permutations...");
	if(top){
		xy<-c(x,y); group<-c(rep(xynames[1],length(x)),rep(xynames[2],length(y)));
		#print(group);
		lr<-surv2.logrank(Surv(xy),factor(group));
		if(lr$pval<=.05) sig.tests<-c(sig.tests,3);
	}
	print(lr);
	cat("\nDoing quantile regression...");
	if(top){
		i<-qstart;qreg.sum<-0;class(qreg.sum)<-"try-error";
		while(class(qreg.sum)=="try-error"){
			if(i>length(quantiles)){
				qreg.sum<-"Unable to perform quantile regression.";
				next;
			}
			#cat('\nTrying quantiles...',paste(quantiles[[i]],collapse=" "));
			qreg<-rq(log(xy)~factor(group,levels=xynames),tau=quantiles[[i]]);i<-i+1;
			qreg.sum<-try(summary(qreg,se=qrse)); cat('.');
		}
		if(class(qreg.sum)!="character"){
			qreg.tab<-c(); qreg.sig=F;
			for(i in qreg.sum){
				qreg.tab<-rbind(qreg.tab,c(quantile=i$tau,i$coefficients[2,]));
				if(min(qreg.tab[,5])<=.05){qreg.sig=T;}
			}
			if(qreg.sig){sig.tests<-c(sig.tests,4);}
		}else{qreg.tab<-qreg.sum;}
	}
	cat('\n'); print(qreg.tab);

	cat("\nFinding the best mortality models..."); 
	if(top & dohz){xd<-demog(x,demint); yd<-demog(y,demint);mod<-findpars(x,y);}
	if(top){
		sigzsc<-zsc[which(zsc$sig=="*"),c(2,4,8)];
		if(class(qreg.tab)!="character"){
			sigqreg<-matrix(qreg.tab[qreg.tab[,5]<.05&qreg.tab[,3]>0,c(1:2,5)],ncol=3);
	} else {sigqreg<-matrix(NA,ncol=3);}
	if(dim(sigqreg)[1]==0){sigqreg<-matrix(NA,ncol=3);}
	lz<-dim(sigzsc)[1];lqr<-dim(sigqreg)[1];
	report<-c(n=length(x),mean=mean(x),round(quantile(x,c(.5,.9)),0),max=max(x),`log-rank`=lr$stat,NA,
		NA,if(match(NA,sigqreg[,1],no=F)){NA}else{quantile(x,sigqreg[,1])},sigzsc[,1]);
	report<-cbind(report,NA,NA,c(length(y),mean(y),round(quantile(y,c(.5,.9)),0),max(y),NA,NA,NA,if(match(NA,sigqreg[,1],no=F)){NA}else{quantile(y,sigqreg[,1])},sigzsc[,2]),NA,NA);
	report[2,2]<-mean(x)-sd(x)/sqrt(length(x)); report[2,3]<-mean(x)+sd(x)/sqrt(length(x));
	report[2,5]<-mean(y)-sd(y)/sqrt(length(y)); report[2,6]<-mean(y)+sd(y)/sqrt(length(y));
	report[3:4,2:3]<-qci(x,c(.5,.9)); report[3:4,5:6]<-qci(y,c(.5,.9));
	if(match(NA,sigqreg[,1],no=F)==0 & dim(report)[1]>10){
		#browser();
		# these are terrible hacks, should be fixed
		inrange<-9:(8+dim(sigqreg)[1]);
		report[inrange,2:3]<-qci(x,sigqreg[,1])[1:length(inrange),];
		report[inrange,5:6]<-qci(y,sigqreg[,1])[1:length(inrange),];
	}
	rownames(report)[3]<-'median';
	if(match(NA,sigqreg[,1],no=F)==0){
		#rownames(report)[dim(report)[1]:(10+dim(sigqreg)[1])]<-paste("(qr) age, quantile",sigqreg[,1]);
	}
	if(dim(sigzsc)[1]>0&dim(report)[1]>11){
		rownames(report)[(11+dim(sigqreg)[1]):dim(report)[1]]<-paste("(score) # alive,  quantile",rownames(sigzsc));
	}
	#colnames(report)[c(1]<-c(xynames[1],'LBx95','UBx95','xynames[2]','LBy95','UBy95');
	report<-cbind(report,p=
	   c(NA,tt$p.value,zsc[match(c(.5,.9),rownames(zsc)),8],NA,lr$pval,NA,
	   NA,
	   #if(dohz){xymod.tab[,5]} else {NA},
	   NA,sigqreg[,3],sigzsc[,3]));
	}
	bndlbs<-c('LB95','UB95');
	colnames(report)<-c(xynames[1],paste(xynames[1],bndlbs,sep='.'),
			    xynames[2],paste(xynames[2],bndlbs,sep='.'),'p');

	if(top){xysurvfit<-survfit(Surv(xy)~group);}
	tmain<-paste(xynames,collapse=" vs. ");
	cat("\nPlotting survival curves.\n");
	plot(xysurvfit,lwd=slwd,bty="n",xlim=c(0,xlim),col=scol,cex.axis=sxcex,las=1,yscale=100);
	axis(1,lwd=sxlwd,labels=F);axis(2,lwd=sxlwd,labels=F);
	for(i in 1:2){points(xysurvfit[i],pch=spch[i],cex=spcex[i],lwd=splwd[i],bg=spbg[i]);}
	title(main=paste(tmain,"Survival"));
	if(sleg){legend("topright",legend=xynames[2:1],pch=spch,y.intersp=1.5,
		bty="n",pt.cex=spcex,cex=sxcex,lwd=splwd,pt.bg=spbg,col=scol);}
	
	if(save){graphsave(path,xynames,".srv.eps",prompt);}else{readline("Press enter to continue.")}
	if(class(qreg.sum)!="character"){
		cat("\n\nPlotting quantile regression curves.\n");
		plot(qreg.sum,parm=2,ols=F,ylim=c(-1,1),main=paste(tmain,"Quantile Regression"));
		points(sigqreg[,1],sigqreg[,2],col='red',lwd=3);
		if(save){graphsave(path,xynames,".qr.eps",prompt);}else{readline("Press enter to continue.")}
	}

	#if(top & dohz){
	#	xhzpars<-rbind(c(mod$g$group0$estimate,NA),c(mod$gm$group0$estimate));
	#	yhzpars<-rbind(c(mod$g$group1$estimate,NA),c(mod$gm$group1$estimate));}
# 	cat("\n\nPlotting hazard curves.\n");
# 	#lazyhazplots(list(xd1,xd7,xd30),xhzpars,c(1,7,30),xlim=xlim);
# 	browser();
# 	lazyhazplots(list(xd),xhzpars,demint,xlim=xlim/30,cols='black',main=paste(tmain,"Hazard"));
# 	lazyhazplots(list(yd),yhzpars,demint,xlim=xlim/30,cols='darkred',add=T);
# 	if(save){graphsave(path,xynames[2],".hz.eps",prompt);}else{readline("Press enter to continue.")}

	if(top){
		sig.tests<-tests[sig.tests]; 
		if(length(sig.tests)==0){sig.tests<-0;names(sig.tests)<-"None";}
		if(!dohz){xd<-yd<-mod<-c();}
		out<-list(smry=smry,zsc=zsc,x=x,y=y,xy=xy,xynames=xynames,path=path,
			group=group,lr=lr,tt=tt,tests=tests,sig.tests=sig.tests,
			qreg=qreg,qreg.tab=qreg.tab,qreg.sum=qreg.sum,demint=demint,
			xd=xd,yd=yd,mod=mod,report=report,sigqreg=sigqreg,
			sigzsc=sigzsc,xysurvfit=xysurvfit,
			slwd=slwd,scol=scol,spch=spch,spcex=spcex,splwd=splwd,
			spbg=spbg,sleg=sleg,sxcex=sxcex,sxlwd=sxlwd,dohz=dohz)
		class(out)<-"survomatic";
		go<-get(as.character(sys.call()[1]));
		out$go<-go; out$sys<-function(q){print(group);print(sys.status());print(sys.call());}
		if(save){
			rname<-paste(paste(xynames,collapse=""),"surv",sep=".");
			fname<-paste(path,rname,".rdata",sep="");
			assign(rname,out);
			.First<-function(x=rname){
				library(utils);
                                ltry<-try(library(Survomatic));
                                if(class(ltry)=='try-error'){
                                        install.packages('Survomatic',repos=c('http://rforge.net/','http://cran.r-project.org/'));}
                                ltry<-try(library(Survomatic));
                                if(class(ltry)=='try-error'){
                                        cat('\nProblem installing library. Functions may not work but
your data is still safe. See if you have a problem with
your internet connection. If it\'s working properly and
you still get this message, contact
alex.bokov@gmail.com\n');} else {.fn<<-get(x);.fn$go();rm(.fn,pos=parent.frame())}
                        }
			save(list=c(rname,'.First'),file=fname);
			cat("\nYour results have been saved as ",fname);
		}
	}

	cat("\nThe following test/s may have produced significant results:\n");
	print(names(sig.tests));
	cat("\nExecutive Summary:\n");
	print(report[1:5,],na.print="");
	cat('\n-----------------------\nLog-rank: p =',report[6,7],'\n');
	#cat('\n-----------------------\nMortality\n');
	#print(report[8:9,],na.print="");
	cat('\n-----------------------\nQuantile Regression & Z-score\n');
	print(report[9:dim(report)[1],],na.print="");
	if(dohz){
		cat('\n-----------------------\nMortality Models\n');
		if(!mod[mod$model=='gm','sig?']){
			mod[mod$model=='gma','sig?']<-mod[mod$model=='gmb','sig?']<-
			mod[mod$model=='gmc','sig?']<-F;
		}
		if(!mod[mod$model=='l','sig?']){
			mod[mod$model=='la','sig?']<-mod[mod$model=='lb','sig?']<-
			mod[mod$model=='ls','sig?']<-F;
		}
		if(!mod[mod$model=='l','sig?']&!mod[mod$model=='g','sig?']){
			mod[mod$model=='lma','sig?']<-F;
		}
		if(!mod[mod$model=='lm','sig?']){
			mod[mod$model=='lma','sig?']<-mod[mod$model=='lmb','sig?']<-
			mod[mod$model=='lmc','sig?']<-mod[mod$model=='lms','sig?']<-F;
		}
		print(subset(mod,`sig?`==T,select=c('MLE','a1','b1','c1','s1','a2','b2','c2','s2',
			'LR','AIC','BIC','p (chi squared)')),na.print="");
	}
	cat("\n\nMethods:");
	cat("\nMeans of log-transformed survival times were compared using the Student\'s 
t-test. Quantiles, including the median, were compared using a modified 
version of the score test described in Wang et. al. (2004). Quantile 
regression was done as described by Koenker and Geling (2001). The 
log-rank test (Fleming and Harrington, 1991) was used to compare 
distributions of survival times with the p-value adjusted based on 2000 
permutations of the data (Andersen et. al., 1993). The mortality 
parameters were fitted and compared using the maximum likelihood method 
adapted from Pletcher et. al. (2000). The R statistical language (R 
Development Core Team, 2008) was used for all calculations except where 
specified. The R scripts used to perform these calculations are 
available from the authors upon request.\n");
	cat("\nThank you for using Survomatic! "); 
	if(top&save) cat("All the analysis you have done will now be 
saved for later use as an R file entitled
",fname,".")
	cat("\n\n\n");
	options('warn'= 0);options(show.error.messages=T);
	if(top) invisible(out);
}

