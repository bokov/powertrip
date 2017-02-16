`ezz` <-
function(d,d2=NULL,g,c1=NULL,c2=NULL,quant=.9,step=.01,thresh=.05,
	      gnames=NULL,noob=F,thetas=c(.05,.95),debug=F,silent=F){
	if(!silent) tstart<-proc.time()[3];
	if(is.null(d2)){
		groups<-levels(as.factor(g));
		d1<-subset(d,g==groups[1])[,1];
		d2<-subset(d,g==groups[2])[,1];
	} else { d1<-d; d<-c(d1,d2);
	if(is.null(g)){groups<-c('a','b');} else {groups<-g;}}
	#mintheta<-min(1/d1,1/d2);
	n1<-length(d1); n2<-length(d2); n<-n1+n2;
	if(is.null(c1)){c1<-rep(T,n1)} else {c1<-as.logical(c1)}; 
	if(is.null(c2)){c2<-rep(T,n2)} else {c2<-as.logical(c2)};
	d1<-d1[c1]; d2<-d2[c2];
	output<-c(); q<-quantile(d,quant);
	if(debug) cat("\nRunning zptest on x1 and x2.");
	x1<-x2<-c();
	for(i in q){x1<-c(x1,sum(d1>i)); x2<-c(x2,sum(d2>i));}
	output<-zptest(x1,x2,n1,n2,n);
	output<-cbind(q,output);
	if(is.null(gnames)){g1<-groups[1];g2<-groups[2];}
	else{g1<-gnames[1];g2<-gnames[2];}
	surv1<-paste("srv",g1,sep='.');
	surv2<-paste("srv",g2,sep='.');
	surv1n1<-paste("srv",g1,"frc",sep='.');
	surv2n2<-paste("srv",g2,"frc",sep='.');
	colnames(output)<-c('age',surv1,surv1n1,surv2,surv2n2,'survtotal.frc','zp');
	output<-data.frame(output);
	output$zp[is.nan(output$zp)]<-9999; #print(output$zp);

	p<-zptests<-c();
	if(debug) cat("\nRunning zptest on k n1\'s and n2\'s.");
	zptests<-outer(0:n1,0:n2,FUN="zptest",n1=n1,n2=n2,n=n,zponly=T);
	zptests[is.nan(zptests)]<-9999;

	for(i in 1:dim(output)[1]){
		keep<-abs(zptests)>abs(output$zp[i]);
		sup<-c(); if(!silent) cat('.');
		for(j in seq(thetas[1],thetas[2],step)){
			if(debug) cat("\nRunning zpprob on n1\'s and n2\'s.");
			zptemp<-outer(0:n1,0:n2,FUN="zpprob",n1=n1,n2=n2,th=j);
			if(any(is.nan(zptemp))) browser();
			sup<-c(sup,sum(zptemp[keep]));
		}
		p<-c(p,max(sup));
	}
	sig<-p<thresh; sig[sig==T]<-"*"; sig[sig==F]<-"";
	output<-cbind(output,p,sig);
	rownames(output)<-quant;
	if(noob){output<-output[,c(1:5,8,9)];}
	if(!silent) print(proc.time()[3]-tstart);
	output;
}

