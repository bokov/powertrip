# todo:
#	have all the repeating params get set to default values inside the case statements for maximum flexibility
#	make the abline for the LR plot less ugly
#	add lty
#	user-configurable legend location
#	make all these customization options settable through srvgui

`plot.fp`<-function(d, what=c('srv','haz','den'), real_or_fit=NULL, model=d$jointmodels[!grepl('w',d$jointmodels)], par='a', type='l', col=c(1,'darkred','gray','red'), bg=col, lwd=rep(2,4), pch=c(24:25,NA,NA), ages=0:max(d$xy), leg=list(cex=.85), groupnames=d$groupnames, lrexp=1.2, lrmin=1e-16, lrsmooth=1000, smooth=10, ...){
    #previous=NULL, 
    myargs<-list(...); if(length(col)<4) col<-rep(col,length=4);
    what<-tolower(what); 
    if(any(what %in% c('l','lr','lrs'))){
	cons<- 2-(c('a','b','c','s') %in% par);
	mylrs<-d$pers[model,cons[1],cons[2],cons[3],cons[4],'LR',];
	mydf<-d$pers[model,cons[1],cons[2],cons[3],cons[4],'df',1];
	if(all(is.na(mylrs))){
	    plot.new();text(.5,.5,'Data not available.'); 
	    print(paste('For model',model,'and parameter/s',par,'no log-ratio estimates are available. They might not have been calculated when the input data,',as.character(match.call())[2],', was created. Please check to make sure that you really intend to compare the parameter/s you specified between the two experimental groups. If you do, you will need to manually run the fp() function on your original survival data with the optional argument \'cons\' set to \'list(',model,'=matrix(c(1,1,1,1,',paste(cons-1,sep=','),'),ncol=4,byrow=T,dimnames=list(NULL,par=c(\'a\',\'b\',\'c\',\'s\')))\' ). Please contact the package maintainers and we will consider adding your desired combination of parameters to the defaults settings in a future version of Survomatic.'));
	} else {
	    if(is.null(myargs$xlab)) myargs$xlab<-'Log Ratio';
	    if(is.null(myargs$main)){
		myargs$main<-paste('Observed vs. Permuted Log-Ratios for Parameter:',paste(par,collapse=','));
	    }
	    if(is.null(myargs$xlim)) myargs$xlim<-c(0,max(mylrs)*lrexp);
	    mysmooth<-seq(lrmin,myargs$xlim[2],len=lrsmooth);
	    if(is.null(myargs$breaks)) myargs$breaks<-100;
	    if(is.null(myargs$border)) myargs$border<-col[1];
	    myargs$col<-bg[1]; myargs$lwd<-lwd[1];
	    myargs$x<-mylrs[-1];
	    do.call(hist,myargs);
	    lines(mysmooth,d$perms*dchisq(mysmooth,df=mydf),col=col[2],lwd=lwd[2])
	    abline(v=mylrs[1],col=col[3],lwd=lwd[3]);
	    if(!is.na(leg)&&!is.null(leg)){
		leg$col<-c(leg$col,col[1:3]);
		leg$lwd<-c(leg$lwd,lwd[1:3]);
		leg$legend<-c(leg$legend,c('Empirical distribution of log-ratios','Chi-squared distribution','Actual log-ratio'));
		leg$x<-'top';
		do.call(legend,leg);
	    }
	}
	# print the empirical hazard histograms and exit
    } else {
	what<-match.arg(what);
	myunc<-d$unc; myunc[is.na(myunc)]<-0; mydx<-demog(d$x); mydy<-demog(d$y);
        if(smooth>1){
          smoother<-function(x,s=smooth){
            c(rep(NA,floor((s-1)/2)),
              apply(embed(x,s),1,mean,na.rm=T),
              rep(NA,ceiling((s-1)/2)))
          }
        } else smoother<-function(x,s=smooth) x;
	simmodel<-substr(model[1],1,1);
	if(is.null(myargs$xlab)) myargs$xlab<-'Age';
	if(is.null(myargs$xlim)) myargs$xlim<-c(0,max(ages));
	if(model=='w') hf<-weibhaz else hf<-srvhaz;
	xsrv<-srvshp(ages, a=myunc[model,'a','x'], b=myunc[model,'b','x'], c=myunc[model,'c','x'], s=myunc[model,'s','x'], model=simmodel);
	ysrv<-srvshp(ages, a=myunc[model,'a','y'], b=myunc[model,'b','y'], c=myunc[model,'c','y'], s=myunc[model,'s','y'], model=simmodel);
	xhaz<-hf(ages, a=myunc[model,'a','x'], b=myunc[model,'b','x'], c=myunc[model,'c','x'], s=myunc[model,'s','x']);
	yhaz<-hf(ages, a=myunc[model,'a','y'], b=myunc[model,'b','y'], c=myunc[model,'c','y'], s=myunc[model,'s','y']);
	xden<-xsrv*xhaz*length(d$x);
	yden<-ysrv*yhaz*length(d$y);
	switch(what,
	    'srv'={
		if(is.null(myargs$ylab)) myargs$ylab<-'Fraction Alive';
		if(is.null(myargs$main)) myargs$main<-'Survival';
		if(is.null(real_or_fit)) real_or_fit<-'both';
		if(real_or_fit %in% c('real','both')){
		    myargs$x<-d$x; myargs$y<-d$y; myargs$cx<-d$cx; myargs$cy<-d$cy; 
		    myargs$col<-col[1:2]; myargs$lwd<-lwd[1:2]; myargs$pch<-pch[1:2]; myargs$bg<-bg[1:2];
		    do.call(plotsurv,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[1:2]);
			leg$lwd<-c(leg$lwd,lwd[1:2]);
			leg$pch<-c(leg$pch,pch[1:2]);
			leg$pt.bg<-c(leg$pt.bg,bg[1:2]);
			leg$legend<-c(leg$legend,groupnames);
		    }
		}
		if(real_or_fit!='real'){
		    myargs$x<-ages; myargs$col<-col[3]; myargs$lwd<-lwd[3]; myargs$pch<-pch[3]; myargs$bg<-bg[3];
		    myargs$y<-xsrv;
		    if(real_or_fit=='both') {pf<-lines; myargs$cx<-myargs$cy<-NULL;} else { pf<-plot; myargs$type<-'l';};
		    do.call(pf,myargs); do.call(points,myargs);
		    myargs$type<-NULL; 
		    myargs$col<-col[4]; myargs$lwd<-lwd[4]; myargs$pch<-pch[4]; myargs$bg<-bg[4];
		    myargs$y<-ysrv;
		    do.call(lines,myargs); do.call(points,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[3:4]);
			leg$lwd<-c(leg$lwd,lwd[3:4]);
			leg$pch<-c(leg$pch,pch[3:4]);
			leg$pt.bg<-c(leg$pt.bg,bg[3:4]);
			leg$legend<-c(leg$legend,paste(groupnames,',',toupper(model),'predicted'));
		    }
		 }
		 if(!is.na(leg)&&!is.null(leg)){
		    leg$x<-'bottomleft';
		    do.call(legend,leg);
		 }
	    },
	    'haz'={
		if(is.null(myargs$ylab)) myargs$ylab<-'Death Rate';
		if(is.null(myargs$main)) myargs$main<-'Hazard';
		if(is.null(real_or_fit)) real_or_fit<-'both';
                xx<-smoother(mydx$ux); yy<-smoother(mydy$ux);
		if(is.null(myargs$ylim)) myargs$ylim<-c(0,max(xx,yy,na.rm=T));
		if(real_or_fit %in% c('real','both')){
		    myargs$type<-'l';
		    myargs$col<-col[1]; myargs$lwd<-lwd[1]; myargs$pch<-pch[1]; myargs$bg<-bg[1];
		    myargs$x<-mydx$time; myargs$y<-xx;
		    do.call(plot,myargs); myargs$type<-NULL; do.call(points,myargs);
		    myargs$col<-col[2]; myargs$lwd<-lwd[2]; myargs$pch<-pch[2]; myargs$bg<-bg[2];
		    myargs$x<-mydy$time; myargs$y<-yy;
		    do.call(lines,myargs); do.call(points,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[1:2]);
			leg$lwd<-c(leg$lwd,lwd[1:2]);
			leg$pch<-c(leg$pch,pch[1:2]);
			leg$pt.bg<-c(leg$pt.bg,bg[1:2]);
			leg$legend<-c(leg$legend,groupnames);
		    }
		}
		if(real_or_fit!='real'){
		    myargs$x<-ages;
		    myargs$col<-col[3]; myargs$lwd<-lwd[3]; myargs$pch<-pch[3]; myargs$bg<-bg[3];
		    myargs$y<-xhaz; 
		    if(real_or_fit=='both') pf<-lines else { pf<-plot; myargs$type<-'l'; };
		    do.call(pf,myargs); myargs$type<-NULL; do.call(points,myargs);
		    myargs$col<-col[4]; myargs$lwd<-lwd[4]; myargs$pch<-pch[4]; myargs$bg<-bg[4];
		    myargs$y<-yhaz;
		    do.call(lines,myargs); do.call(points,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[3:4]);
			leg$lwd<-c(leg$lwd,lwd[3:4]);
			leg$pch<-c(leg$pch,pch[3:4]);
			leg$pt.bg<-c(leg$pt.bg,bg[3:4]);
			leg$legend<-c(leg$legend,paste(groupnames,',',toupper(model),'predicted'));
		    }
		 }
		 if(!is.na(leg)&&!is.null(leg)){
		    leg$x<-'topleft';
		    do.call(legend,leg);
		 }
	    },
	    'den'={
		if(is.null(myargs$ylab)) myargs$ylab<-'Number of Deaths';
		if(is.null(myargs$main)) myargs$main<-'Density Function';
		if(is.null(real_or_fit)) real_or_fit<-'both';
                xx<-smoother(mydx$deaths); yy<-smoother(mydy$deaths);
		if(is.null(myargs$ylim)) myargs$ylim<-c(0,max(xx,yy,na.rm=T));
		if(real_or_fit %in% c('real','both')){
		    myargs$type<-'l';
		    myargs$col<-col[1]; myargs$lwd<-lwd[1];
                    myargs$pch<-pch[1]; myargs$bg<-bg[1];
		    myargs$x<-mydx$time; myargs$y<-xx;
		    do.call(plot,myargs); myargs$type<-NULL; do.call(points,myargs);
		    myargs$col<-col[2]; myargs$lwd<-lwd[2];
                    myargs$pch<-pch[2]; myargs$bg<-bg[2];
		    myargs$x<-mydy$time; myargs$y<-yy;
		    do.call(lines,myargs); do.call(points,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[1:2]);
			leg$lwd<-c(leg$lwd,lwd[1:2]);
			leg$pch<-c(leg$pch,pch[1:2]);
			leg$pt.bg<-c(leg$pt.bg,bg[1:2]);
			leg$legend<-c(leg$legend,groupnames);
		    }
		}
		if(real_or_fit!='real'){
		    myargs$x<-ages;
		    myargs$col<-col[3]; myargs$lwd<-lwd[3]; myargs$pch<-pch[3]; myargs$bg<-bg[3];
		    myargs$y<-hf(ages, a=myunc[model,'a','x'], b=myunc[model,'b','x'], c=myunc[model,'c','x'], s=myunc[model,'s','x']); 
		    myargs$y<-xden;
		    if(real_or_fit=='both') pf<-lines else { pf<-plot; myargs$type<-'l'; };
		    do.call(pf,myargs); myargs$type<-NULL; do.call(points,myargs);
		    myargs$col<-col[4]; myargs$lwd<-lwd[4]; myargs$pch<-pch[4]; myargs$bg<-bg[4];
		    myargs$y<-yden;
		    do.call(lines,myargs); do.call(points,myargs);
		    if(!is.na(leg)&&!is.null(leg)){
			leg$col<-c(leg$col,col[3:4]);
			leg$lwd<-c(leg$lwd,lwd[3:4]);
			leg$pch<-c(leg$pch,pch[3:4]);
			leg$pt.bg<-c(leg$pt.bg,bg[3:4]);
			leg$legend<-c(leg$legend,paste(groupnames,',',toupper(model),'predicted'));
		    }
		 }
		 if(!is.na(leg)&&!is.null(leg)){
		    leg$x<-'topleft';
		    do.call(legend,leg);
		 }
	    }
	)
    }
}

# for LR only-- whichpars: 'a','b','c,','s', combos of, c('a','b') by default
# which: empirical, fitted
# type: points, lines
# overplot: T, F
# previous plot params
# color, lwd, lty, pch, pcx, border, bg, labels, main, xaxis, yaxis (left/right), axis names
