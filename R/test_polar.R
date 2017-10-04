#' ---
#' title: "Test Functions for PowerTrip"
#' author: "Alex F. Bokov"
#' date: "10/03/2017"
#' ---

#' Generates a random binomial value whose probability of being
#' 1 is a bivariate normal distribution scaled such that it is unity
#' at 0,0 and diminishes asymptotically in all directions from there
#' The purpose of this is to generate an underlying smooth, 
#' "radially monotonic" distribution for which only binomial outcomes
#' be observed, for testing powertrip functionality in isolation from
#' any specific statistical tests, etc.
#' @param xx numeric, the closer to 0, the more likely to return 1
#' @param yy numeric, the closer to 0, the more likely to return 1
gen_binom<-function(xx,yy){
  rbinom(1,1,mvtnorm::dmvnorm(c(xx,yy))/mvtnorm::dmvnorm(c(0,0)))
}

## convert cartesian to polar coordinates
## (why isn't this part of base R?!)
## returns a matrix of same dimension
crt2pol <- function(xx,center=kmeans(xx,1)$centers,sort=c('angle','radius','none')){
  ## xx = cartesian coordinates matrix [x,y]
  ## center = the center relative to which the polar coordinates will be returned
  ## sort = how/whether to sort the coordinates
  ##        (so the output matrix forms a closed shape when used by lines(...) )
  cxx<-scale(xx,center=center,scale=F);
  ## note the multiplication by sign--
  ## if not done, results will only be correct for x > center[1]
  oo<-cbind(r=sqrt(rowSums(cxx^2))*ifelse(sign(cxx[,1]),sign(cxx[,1]),1),theta=atan(cxx[,2]/cxx[,1]));
  ## retain the center for later use
  attr(oo,'center')<-center;
  ## preliminary testing suggests the below will work even in the degenerate case where
  ## all oo[,1] has the same sign, so don't worry about that for now
  oo[]<-switch(match.arg(sort),
               angle={
                 ooa<-oo[oo[,1]>=0,];oob<-oo[oo[,1]<0,];
                 ooa<-ooa[order(ooa[,2],decreasing=T),];
                 oob<-oob[order(oob[,2],decreasing=T),];
                 rbind(oob,ooa);},
               radius=oo[order(oo[,1]),],
               none=oo);
  oo;
}

## convert polar coordinates back to cartesian
## returns a matrix of same dimension
pol2crt <- function(pp,center=if('center'%in%names(attributes(pp))) as.numeric(attr(pp,'center')) else rep(0,len=ncol(pp))){
  ## pp = polar coordinate matrix [r,theta]
  ## center = point on which the cartesian result should be centered
  oo<-scale(cbind(x=pp[,1]*cos(pp[,2]),y=pp[,1]*sin(pp[,2])),center= -center,scale=F);attr(oo,'center')<-center;oo;
}

#' ### Here we try it
#' 
#' Sample from a polar space immediately converting to cartesian
sampled<-pol2crt(cbind(r=runif(1000,-4,4),theta=runif(100,-pi/2,pi/2)));
#' generate binomial outcomes (we are going to look for the 0.8 contour)
sampled<-cbind(sampled,res=apply(sampled,1,function(zz) gen_binom(zz[1],zz[2])));
#' Fit a model
prmod <- glm(res~x*y,data.frame(sampled),family='binomial');
sampled0 <- pol2crt(cbind(r=runif(1000,-4,4),theta=runif(100,-pi/2,pi/2)));
sampled0 <- sampled0[with(predict(update(prmod,data=data.frame(sampled)),data.frame(sampled0),type='response',se.fit=T),which(abs(.GlobalEnv$.dbg<-(fit-0.2))<0.1&se.fit>quantile(se.fit,.5))),];
sampled <- rbind(sampled,cbind(sampled0,res=apply(sampled0,1,function(zz) gen_binom(zz[1],zz[2]))));
plot(sampled,pch='.');
