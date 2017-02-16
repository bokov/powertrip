smoother<-function(x,s=smooth){
  c(rep(NA,floor((s-1)/2)),
    apply(embed(x,s),1,mean,na.rm=T),
    rep(NA,ceiling((s-1)/2)))
}  
