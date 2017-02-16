`plot.survomatic` <- function(xx,what=c('srv','haz','den'),col=c('black','black','darkred','darkred'),lwd=c(2,2,2,2),lty=c(1,2,1,2),xlim=NULL){
  if(is.null(xlim)) xlim <- with(xx,c(0,max(x.d$Days,y.d$Days)));
  if('srv' %in% what){
    with(xx,{plot(survfit(Surv(c(x,y),c(cx,cy))~rep(1:2,c(length(x),length(y)))),col=col[c(1,3)],lwd=lwd[1:2],xlab='Age, Days',ylab='Fraction Alive')});
  }
  if('haz' %in% what){
    with(xx,{plot(`Smoothed Hazard`~Days,data=x.d,type='l',xlim=xlim,lwd=lwd[1],lty=lty[1],col=col[1],
                  ylim=c(0,max(x.d$`Predicted Hazard`,x.d$`Smoothed Hazard`,y.d$`Predicted Hazard`,y.d$`Smoothed Hazard`,na.rm=T)));
             lines(`Predicted Hazard`~Days,data=x.d,lwd=lwd[2],lty=lty[2],col=col[2]);
             lines(`Smoothed Hazard`~Days,data=y.d,lwd=lwd[3],lty=lty[3],col=col[3]);
             lines(`Predicted Hazard`~Days,data=y.d,lwd=lwd[4],lty=lty[4],col=col[4]);});
  }

  if('den' %in% what){
    with(xx,{plot(`Smoothed Deaths`~Days,data=x.d,type='l',xlim=xlim,lwd=lwd[1],lty=lty[1],col=col[1],
                  ylim=c(0,max(x.d$`Predicted Deaths`,x.d$`Smoothed Deaths`,y.d$`Predicted Deaths`,y.d$`Smoothed Deaths`,na.rm=T)));
             lines(`Predicted Deaths`~Days,data=x.d,lwd=lwd[2],lty=lty[2],col=col[2]);
             lines(`Smoothed Deaths`~Days,data=y.d,lwd=lwd[3],lty=lty[3],col=col[3]);
             lines(`Predicted Deaths`~Days,data=y.d,lwd=lwd[4],lty=lty[4],col=col[4]);});
  }
}
