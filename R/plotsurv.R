#`plotsurv.old`<-
#function(x,y,cx=NULL,cy=NULL,col=c(1,'darkred'),lwd=2,legend=NULL,lloc='bottomleft',pcex=1,pch=24:25,ylim=0:1,bg=NULL,add=F,...){
#	lx<-length(x);ly<-length(y);args<-list(...);
#	if(is.null(cx)){cx<-rep(1,lx);}; if(is.null(cy)){cy<-rep(1,ly);};
#	if(is.null(bg)){bg<-col;}
#	sf<-survfit(Surv(c(x,y),event=c(cx,cy))~c(rep(0,lx),rep(1,ly)));
#	plot(sf,col=col,lwd=lwd,ylim=ylim,...);
#	# not a great way to do it, but will do for now
#	if(!is.null(args$fun)){
#		if(args$fun=='event')for(i in 1:2)
#			{points(sf[i]$time,1-sf[i]$surv,pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=col[i])}
#	} else {for(i in 1:2){points(sf[i],pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=col[i]);}}
#	if(!is.null(legend)){legend(x=lloc,col=col,pt.bg=bg,lwd=lwd,pch=pch,legend=legend);}
#	invisible(sf);
#}

`plotsurv`<-
function(x, y, cx=NULL, cy=NULL, lloc='bottomleft', add=F, pcex=1, xscale=1, pch=c(NA,NA), legend=NULL, ...){
# parse arguments
lx<-length(x);ly<-length(y);myargs<-list(...);
# set some defaults if not explicitly set by arguments
if(is.null(cx)){cx<-rep(1,lx);}; if(is.null(cy)){cy<-rep(1,ly);};
if(is.null(myargs$bg)) myargs$bg<-c(1,'darkred');
if(is.null(myargs$col)) myargs$col <- myargs$bg;
if(is.null(myargs$lwd)) myargs$lwd<-2;
#if(is.null(myargs$pcex)) myargs$pcex<-1;
#if(is.null(myargs$pch)) myargs$pch<-24:25;
#if(is.null(myargs$pch)) myargs$pch <- c(NA,NA);
if(is.null(myargs$ylim)) myargs$ylim<-0:1;
# generate data to be plotted
sf<-survfit(Surv(c(x,y),event=c(cx,cy))~c(rep(0,lx),rep(1,ly)));
myargs$x<-sf;
# set and call the plotting function
if(add) {do.call(lines, myargs)} else {
    do.call(plot, myargs);
}
iargs<-myargs;
if(!is.null(iargs$fun)){
	if(iargs$fun=='event')for(i in 1:2)
		with(iargs, points(sf[i]$time, 1-sf[i]$surv,pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=col[i],xscale=xscale))
} else {for(i in 1:2){
	with(iargs,points(sf[i],pch=pch[i],cex=pcex,lwd=lwd,bg=bg[i],col=bg[i],xscale=xscale));
    }
}
if(!is.null(legend)){
    with(myargs,legend(x=lloc,col=col,pt.bg=bg,lwd=lwd,pch=pch,legend=legend));
}
invisible(myargs);
}
