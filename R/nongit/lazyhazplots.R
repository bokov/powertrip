`lazyhazplots` <-
function(ts,pars,ints,cols=rainbow(length(ts)),xlim=NULL,add=F,main="",lwd=2){
	mnt<-mnx<-mx<-c(); 
	for(j in ts){
		mx<-c(mx,max(j[j$ux>0&is.finite(j$ux),]$ux,na.rm=T)); 
		mnx<-c(mnx,min(j[j$ux>0&is.finite(j$ux),]$ux,na.rm=T));}
	mnt<-match(T,ts[[which.max(ints)]]$ux>0)
	mx<-max(mx); mnx<-min(mnx); if(is.null(xlim)){xrng<-c(mnt,max(ts[[which.min(ints)]]$time));}else{
		xrng<-c(10,xlim);}
	#print(xrng);
	if(add){points(1+ts[[1]]$time,ts[[1]]$ux,pch=3,col=cols[1],lwd=lwd);} else {
		plot(1+ts[[1]]$time,ts[[1]]$ux,log='y',type='p', pch=3, lwd=lwd,
			xlim=xrng,ylim=c(mnx,mx),xlab="Time",ylab="Hazard Rate",col=cols[1],main=main);
	} 
	lines(1:xrng[2],ints[1]*srvhaz(1:xrng[2],pars[1,1],pars[1,2],i=ints[1]),col=cols[1],lty=2,lwd=lwd);
	lines(1:xrng[2],ints[1]*srvhaz(1:xrng[2],pars[2,1],pars[2,2],pars[2,3],i=ints[1]),col=cols[1],lty=1,lwd=lwd);
	xrng<-xrng[1]:xrng[2];
	if(length(ts)>1){
		for(j in 2:length(ts)){
		points(1+ts[[j]]$time,ts[[j]]$ux,col=cols[j],pch=3); 
			lines(1:xrng[2],ints[j]*srvhaz(1:xrng[2],pars[1,1],pars[1,2],i=ints[j]),col=cols[j],lty=2,lwd=lwd);
			lines(1:xrng[2],ints[j]*srvhaz(1:xrng[2],pars[2,1],pars[2,2],pars[2,3],i=ints[j]),col=cols[j],lty=1,lwd=lwd);
		}
	}
}

