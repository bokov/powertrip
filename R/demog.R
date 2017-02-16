`demog` <-
function(t,int=1){
	if(int>1){t<-round(t/int);}
	t<-sort(t); t.rle<-rle(t); d<-px<-lx<-c();
	ts<-1:max(t.rle$values);
	for(i in ts){
		j<-match(i,t.rle$values);
		if(is.na(j)){d<-c(d,0);}else{d<-c(d,t.rle$lengths[j]);}
	}
	for(i in ts) {
		lx<-c(lx,(sum(d)-sum(d[1:i]))/sum(d));
	}
	for(i in ts) {px<-c(px,lx[i+1]/lx[i]);}
	ux<- -log(px); ux[is.infinite(ux)]<-NA; lnux<-log(ux); lnux[is.infinite(lnux)]<-NA;
	out<-data.frame(time=ts,deaths=d,lx=lx,px=px,ux=ux,lnux=lnux);
	return(out);
}

