survwrapper <- function(x,y=NULL,cx=rep(1,length(x)),cy=rep(1,length(y)),forcemodel=NULL,
                        ext=F,n=length(c(x,y)),AIC=F,BIC=F,breakties='AIC',use_g_start=T,hse=T,
                        compare.matrix=NULL,constraint.matrix=NULL,thresh=0.05,smooth=7,silent=F,batch=F){
  # x, y are vectors of survival times; cx and cy are corresponding vectors of censoring
  # models is a character vector with any combination of e, w, g, gm, l, and lm
  # ext is whether to use extended comparisons FOR PARAMETER estimation... probably doesn't work
  # n is the total sample size
  # AIC and BIC are whether to calculate the AIC and BIC statistics
  # use_g_start = whether to do a quick estimate of the initial g parameters using least squares
  # breakties is what to use for choosing a model if there are more than one justified by the comparisons
  # compare.matrix is a matrix for specifying a customized comparison algorithm
  # constraint.matrix is a matrix of 1's and 0's for specifying a customized set of parameter constraints to test
  # thresh is the significance cutoff
  models=c('g','gm','l','lm');
  # First, we fit the parameters for each model of interest and get the log likelihood
  if(use_g_start){
    g_startx<-with(survfit(Surv(x,cx)~1),lm(lnux~time,subset(data.frame(lnux=log(-log(c(surv[-1],NA)/surv)),time),!is.infinite(lnux)))$coef);
    if(!is.null(y)) g_starty<-with(survfit(Surv(y,cy)~1),lm(lnux~time,subset(data.frame(lnux=log(-log(c(surv[-1],NA)/surv)),time),!is.infinite(lnux)))$coef);
    g_startx[1]<-exp(g_startx[1]); g_starty[1]<-exp(g_starty[1]);
  } else g_startx <- g_starty <- NULL;
  x.m <- fitmods(x,cx,models,ext,g_start=g_startx,hse=hse); if(!is.null(y)) y.m <- fitmods(y,cy,models,ext,g_start=g_starty,hse=hse) else y.m <- NULL;
  # Then, using likelihood ratios, the models (or joint models if y is non null) are compared
  xy.sm <- sigmods(x.m[,c('model','LL')],y.m[,c('model','LL')],n=n,breakties=breakties,AIC=AIC,BIC=BIC,cmat=compare.matrix,thresh=thresh);
  #browser();
  if(!is.null(forcemodel)) suggested.models <- forcemodel else suggested.models <- unique(xy.sm$complex[xy.sm$chosen]);
  par.differences <- list();
  # Then, for each suggested model, the significance of constraints for the parameters of interest is calculated
  if(!is.null(y)){
    for(i in suggested.models){
      par.differences[[i]] <- cpars(x,y,xpars=x.m[x.m$model==i,c('a','b','c','s')],ypars=y.m[y.m$model==i,c('a','b','c','s')],
                                    cx=cx,cy=cy,LL=xy.sm[xy.sm$complex==i,'LLc'][1],model=i,const=constraint.matrix);
    }
  }
  if(batch) return(list(x.m=x.m,y.m=y.m,xy.sm=xy.sm,suggested.models=suggested.models,par.differences=par.differences));
  nx <- length(x); ny <- length(y); 
  x.d <- demog(x); if(!is.null(y)) y.d <- demog(y);
  x.d <- with(x.d,data.frame(Days=time,`At Risk (Nx)`=nx-cumsum(deaths),`Deaths (Dx)`=deaths,`Censoring (Cx)`=NA, `Fraction Surviving`=lx, `Fraction Dying`=1-lx,
                             `Empirical Hazard`=ux,`Smoothed Deaths`=smoother(deaths,smooth),`Smoothed Hazard`=smoother(ux,smooth),
                             `Predicted Survivorship`=srvshp(x.d$time,
                               a=x.m[x.m$model==suggested.models[1],'a'],
                               b=x.m[x.m$model==suggested.models[1],'b'],
                               c=x.m[x.m$model==suggested.models[1],'c'],
                               s=x.m[x.m$model==suggested.models[1],'s'],model=suggested.models[1]),
                             `Predicted Hazard`=srvhaz(x.d$time,
                               a=x.m[x.m$model==suggested.models[1],'a'],
                               b=x.m[x.m$model==suggested.models[1],'b'],
                               c=x.m[x.m$model==suggested.models[1],'c'],
                               s=x.m[x.m$model==suggested.models[1],'s']),check.names=F));
  x.d$`Predicted Deaths` <- x.d$`Predicted Survivorship`*x.d$`Predicted Hazard`*nx;

  if(!is.null(y)) {
    y.d <- with(y.d,data.frame(Days=time,`At Risk (Nx)`=ny-cumsum(deaths),`Deaths (Dx)`=deaths,`Censoring (Cx)`=NA, `Fraction Surviving`=lx, `Fraction Dying`=1-lx,
                             `Empirical Hazard`=ux,`Smoothed Deaths`=smoother(deaths,smooth),`Smoothed Hazard`=smoother(ux,smooth),
                             `Predicted Survivorship`=srvshp(y.d$time,
                               a=y.m[y.m$model==suggested.models[1],'a'],
                               b=y.m[y.m$model==suggested.models[1],'b'],
                               c=y.m[y.m$model==suggested.models[1],'c'],
                               s=y.m[y.m$model==suggested.models[1],'s'],model=suggested.models[1]),
                             `Predicted Hazard`=srvhaz(y.d$time,
                               a=y.m[y.m$model==suggested.models[1],'a'],
                               b=y.m[y.m$model==suggested.models[1],'b'],
                               c=y.m[y.m$model==suggested.models[1],'c'],
                               s=y.m[y.m$model==suggested.models[1],'s']),check.names=F));
    y.d$`Predicted Deaths` <- y.d$`Predicted Survivorship`*y.d$`Predicted Hazard`*ny;
    o <- list(x.m=x.m,y.m=y.m,xy.sm=xy.sm,par.differences=par.differences,x=x,y=y,cx=cx,cy=cy,
              x.d=x.d,y.d=y.d,suggested.models=suggested.models,nx=nx,ny=ny);
  } else {
    o <- list(x.m=x.m,x=x,cx=cx,x.d=x.d,nx=nx);
  }
  class(o) <- 'survomatic';
  if(!silent) { print(o); return(o);} else invisible(o);
}
