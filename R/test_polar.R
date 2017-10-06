#' ---
#' title: "Test Functions for PowerTrip"
#' author: "Alex F. Bokov"
#' date: "10/03/2017"
#' ---

#' Rescale cartesian parameters to span the space 
#' between maxs and mins (both of them vectors equal to
#' number of dimensions over which we are simulating)
#' xx is a matrix of cartesian coordinates with as many
#' columns as there are dimensions
#' Apparently we should just leave the ctr parameter alone
rescale <- function(xx,maxs=0,mins=0,ctr=(maxs+mins)/2,smaxs=1,smins=-1,sctr=0){
  sranges <- smaxs - smins;
  ranges <- maxs - mins;
  # first we center and scale to unit circle (cartesian)
  # then by multiplying by ranges and subtracting ctr we get
  # it to the new scale and center
  oo <- t(apply(xx,1,function(zz) ranges*(zz - sctr)/sranges+ctr));
}

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
  rbinom(1,1,mvtnorm::dmvnorm(c(xx,yy))/mvtnorm::dmvnorm(c(0,0)));
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
pol2crt.old <- function(pp
                    ,center=if('center'%in%names(attributes(pp))) as.numeric(attr(pp,'center')) else rep(0,len=ncol(pp))
                    #,maxes=Inf,mins=-Inf,denom=2
                    ){
  ## pp = polar coordinate matrix [r,theta]
  ## center = point on which the cartesian result should be centered
  oo<-scale(cbind(x=pp[,1]*cos(pp[,2]),y=pp[,1]*sin(pp[,2])),center= -center,scale=F);attr(oo,'center')<-center;
  # if(!isTRUE(abs(maxes)==Inf)&!isTRUE(abs(mins)==Inf)){
  #   # find the out-of-range values
  #   oor<-apply(oo,1,function(xx) any(xx>maxes)||any(xx<mins));
  #   oo[oor,]<-oo[sample(which(!oor),sum(oor),rep=T),];
  #   cat('\nReplacing ',sum(oor),' out-of-range values via resample')
  #   rngs <- (maxes - mins);
  #   cat('\nScaling...\n');
  #   oo <- t(apply(oo,1,function(xx) (xx-mins)/rngs));
  # }
  oo;
}

pol2crt <- function(coords,center=0,names){
  # Try to save some effort by accepting both matrices and vectors
  if(isTRUE(nrow(coords)>0)) {
    thisfun <- match.fun(as.character(sys.call()[1]));
    return(t(apply(coords,1,thisfun
              ,center=center,names=names)));
    } else {
      # we assume that the first coordinate is the radial one
      if(missing(names)) names<-make.names(seq_along(coords));
      setNames(coords[1] * 
                 # the rest are assumed to be thetas
                 c(cos(coords[-1]),1) *
                 cumprod(c(1,sin(coords[-1]))) + 
                 center,names);
    }
  }

sample_polar <- function(nn=1000,rmax=3,thetawrap=0.5
                         ,ntheta=1,center=c(0,0),maxes,mins){
  oo <- cbind(r=runif(nn,0,rmax),theta=runif(nn,0,(2+thetawrap)*pi));
  oocr <- pol2crt(oo,center=center);
  if(!missing(maxes)&!missing(mins)){
    oor <- apply(oocr,1,function(xx) any(xx>maxes)||any(xx<mins));
    oo <- oo[!oor,]; oocr <- oocr[!oor,];
  }
  return(list(spol=oo,scrt=oocr));
}

#' ### Here we try it
#' 
#maxes <- c(1.5,1.8); mins<-c(-1.65,-2.3);
maxes <- c(4,4); mins <- c(-4,-4);
#ctr <- c(0.3,-0.2); # okay, so the center parameter works...
ctr <- c(0,0);
#' Sample from a polar space immediately converting to cartesian
#' TODO: create a sample_polar function
spol<-sample_polar(rmax=1,maxes=maxes,mins=mins,center=ctr)$spol;
#spol<-cbind(r=runif(1000,0,3),theta=runif(1000,0,2.5*pi));
#' generate binomial outcomes (we are going to look for the 0.8 contour)
spol<-cbind(spol,res=apply(pol2crt(spol,center = ctr),1,function(zz) gen_binom(zz[1],zz[2])));
#' Fit a model
prmod <- glm(res~r*sin(theta)*cos(theta),data.frame(spol),family='binomial');
logenv <- new.env();
#' The following parts get repeated many times
for(ii in 100000){
  newsmp <- sample_polar(rmax=1,maxes=maxes,mins=mins,center=ctr);
  bestfit <- with(predict(
    update(prmod,data=data.frame(spol))
    ,data.frame(newsmp$spol),type='response',se.fit=T)
    ,{prs <- abs(fit-0.2)^2/se.fit;
    #prs <- abs(fit-0.2);
    which(prs<quantile(prs,.15));});
  newsmp <- sapply(newsmp,function(xx) xx[bestfit,],simplify=F);
  rownames(newsmp$scrt) <- apply(newsmp$scrt,1
                                 ,function(xx) paste0(c('_',xx),collapse='_'));
  data <- sapply(rownames(newsmp$scrt)
                 ,function(ii) ptsim_binom(newsmp$scrt[ii,])
                 ,simplify=F);
  outcome<-c();
  for(ii in names(data)) {
    logenv[[ii]]<-ptpnl_passthru(data[[ii]],coords=newsmp$scrt[ii,]);
  }
  #spol <- rbind(spol,cbind(newsmp$spol,res=apply(newsmp$scrt,1,function(zz) gen_binom(zz[1],zz[2]))));
  spol <- rbind(spol,cbind(newsmp$spol
                           ,res=sapply(names(data),function(xx) logenv[[xx]]$outcome[1])))
}
#spol0 <- cbind(r=runif(1000,0,3),theta=runif(1000,0,2.5*pi));
#scrt0 <- pol2crt(spol0,center=ctr);
#oor <- apply(scrt0,1,function(xx) any(xx>maxes)||any(xx<mins));
#scrt0 <- scrt0[!oor,]; spol0 <- spol0[!oor,];
#spol <- spol[bestfit,]; scrt0 <- scrt0[bestfit,];
#' End repeat/update part
#' 
#' Below are the visualizations that can be done on any iteration
plot(spol[,-3],pch='.',col='#00000050'); #,xlim=c(1.5,2.5));
points(newsmp$spol,pch='+',col='red');
#' How good is the fitted contour?
bar <-abs(predict(update(prmod),type='response')-0.2);
foo<-pol2crt(spol[which(bar<0.01),-3],center=ctr);plot(foo,pch='.',col='#00000099');dim(foo);
#' Run the following once only after the first few iteration, for reference
# bar.bak <- bar; foo.bak <- foo;
points(foo.bak,col='red',pch='+');
#' How does this look on polar coordinates?
plot(spol[which(bar.bak<quantile(bar.bak,.1)),-3],pch='+',col='red',xlim=range(spol[,1]),ylim=range(spol[,2]));
points(spol[which(bar<quantile(bar,.1)),-3]);
#' Note: if the distribution is not logistic, so far it's in a way
#' that does not cause it to be zero-inflated, overdispersed, or 
#' heteroscedastic
#' 
#' Let's see if this holds up when we are using real p-values.
#' But first we need to get bounds and centering to work.
#' 
