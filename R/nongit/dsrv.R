dsrv <- function(xx, aa, bb, cc=0, ss=0, int=1, model='g'){
  srvshp(xx,a=aa,b=bb,c=cc,s=ss,i=int,model=model)*srvhaz(xx,a=aa,b=bb,c=cc,s=ss,i=int);
}

qsrv <- function(qq, aa, bb, cc=0, ss=0, int=1, model='g',lower.tail=T,range=0:3000){
  if(!lower.tail) qq <- 1-qq;
  range[which.min(abs(srvshp(range,a=aa,b=bb,c=cc,s=ss,i=int,model=model)-qq))];
}

predict.surv <- function(xx, cx=rep(1,length(xx)), aa, bb, cc=0, ss=0, int=1, model='g', range=0:3000,resid=F){
  if(missing(aa)|missing(bb)) {
    fm<-subset(fitmods(xx,cx,models=model),model==model)[2,c('a','b','c','s')];
    if(missing(aa)) aa<-fm[1]; if(missing(bb)) bb<-fm[2]; if(cc==0&!is.na(fm[3])) cc<-fm[3]; if(ss==0&!is.na(fm[4])) ss<-fm[4];
  }
  oo<-sapply(xx,function(ii) qsrv(mean(c(Inf,0,xx)>ii),a=aa,b=bb, c=cc, s=ss, i=int, model=model, range=range));
  if(resid) return(cx*(xx-oo)) else return(oo);
}
