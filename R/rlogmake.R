rlogmake<-function(n,a,b,c,s){
  if(c<0) stop("Argument c must be positive");
  x2<-rlgst(n,a,b,s);
  if(c==0) return(x2);
  ##if(c==0) return(rlgst(n,a,b,s));
  x1<-rexp(n,rate=c);
  return(pmin(x1,x2))
}
