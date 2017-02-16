`zptest` <-
function(x1,x2,n1,n2,n,zponly=F){
	that1<-x1/n1; that2<-x2/n2;
	ttil<-(x1+x2)/n;
	denom<-sqrt(ttil*(1 - ttil)*(1/n1 + 1/n2));
	zp<-(that1 - that2)/denom;
	if(zponly){return(zp);}else{
		output<-cbind(x1,that1,x2,that2,ttil,zp);
		return(output);
	}
}

