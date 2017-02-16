srvsrange <-
function(floors,ceilings=4*floors,fq,cq,lb=c(-20,-8,-20,0),ub=c(-3,-1,-3,10),nn=50,experimental=F){
#function(maxage=2800,earlymort=100,medfloor=365,medmid=500,medceil=2000,lowers=c(0,-16,-16,0),uppers=c(1,-1,-1,10),nn=50){
	lf<-length(floors); lc<-length(ceilings);
	sf<-rep(0,length=lf); sc<-rep(1,length=lc); iters<-0;
	#srvs<-memoize(srvshp);
	# automatically generate equally spaced quantiles if not specified
	if(missing(fq)) fq<-seq(1,0,length=lf+2)[-c(1,lf)]; 
	if(missing(cq)) cq<-seq(1,0,length=lc+2)[-c(1,lc)];
	qq<-sort(unique(c(fq,cq)),decr=T);
	#smax<-1;slow<-0;shi<-1;smid<-0;searly<-0;iqr<-0;iters<-0;pp<-c(a=0,b=0,c=0,s=0);
	while(any(sf<fq)|any(sc>cq)){
	#while(smax>0.01|slow<0.75|shi>0.5|searly<0.95|smid<.25){
		iters<-iters+1;
		#dd<-runif(1,lowers[1],uppers[1]);
                bb<-exp(runif(1,lb[2],ub[2]));
		if(experimental){
			dd<-runif(1,lb[1],ub[1]);
			aa<-exp(-bb*dd);
		} else { aa<-exp(runif(1,lb[1],ub[1])); }
		#aa<-max(bb^(1/dd),.Machine$double.eps);
		cc<-exp(runif(1,lb[3],ub[3]));
		ss<-runif(1,lb[4],ub[4]);
		sf<-srvshp(floors,a=aa,b=bb,c=cc,s=ss,model='lm');
		sc<-srvshp(ceilings,a=aa,b=bb,c=cc,s=ss,model='lm');
		#searly<-srvshp(earlymort,a=aa,b=bb,c=cc,s=ss,model='lm'); if(searly<0.95) next;
		#smax<-srvshp(maxage,a=aa,b=bb,c=cc,s=ss,model='lm'); if(smax>0.01) next;
		#slow<-srvshp(medfloor,a=aa,b=bb,c=cc,s=ss,model='lm'); if(slow<0.75) next;
		#smid<-srvshp(medmid,a=aa,b=bb,c=cc,s=ss,model='lm'); if(smid<0.25) next;
		#shi<-srvshp(medceil,a=aa,b=bb,c=cc,s=ss,model='lm'); #if(shi>0.5) next;
	}
	fullsrv<-sapply(qq,function(ii) which.min(abs(srvshp(1:max(ceilings),a=aa,b=bb,c=cc,s=ss,model='lm')-ii)));
	c(a=aa,b=bb,c=cc,s=ss,iters=iters,fullsrv);
}


