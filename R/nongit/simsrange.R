simsrange <-
function(nn,maxage=2800,medfloor=355,medceil=2000,lowers=c(-16,-16,-16,0),uppers=c(-1,-1,-1,10)){
 smax<-Inf;smed<-0;iters<-0;pp<-c(a=0,b=0,c=0,s=0);
 while(smax>maxage|smed>medceil|smed<medfloor){
 pp[1:3]<-exp(runif(3,lowers[-4],uppers[-4])); pp[4]<-runif(1,lowers[4],uppers[4]); xx<-simsurv(nn,'lm',pp);
 smax<-max(xx); smed<-median(xx); iters<-iters+1;
 }
 c(pp,iters=iters,summary(xx));
 }


