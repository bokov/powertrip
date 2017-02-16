cpars <- function(x,y,xpars,ypars,model,LL,const=NULL,cx=rep(1,length(x)),cy=rep(1,length(y))){
  # x and y are vectors of survival times; cx and cy are corresponding censoring verctors
  # xpars and ypars are vectors of length for with starting values for the respective a, b, c, and s
  # model is a character indicating which model
  # LL is the log-likelihood of the unconstrained model
  # const is an k x 4 matrix with a constraint of interest in each row
  # strip out the constraints that are not possible for this model
  if(is.null(const)) const <- matrix(c(1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,1,1,1,1,0,0,0,0,0),nrow=6,byrow=T);
  const[,is.na(xpars)|is.na(ypars)]<-1; const <- unique(const);
  pcols <- c('a1','b1','c1','s1','a2','b2','c2','s2');
  # initialize a data frame
  o <- data.frame(const,df=4-rowSums(const),LL=NA,a1=NA,b1=NA,c1=NA,s1=NA,a2=NA,b2=NA,c2=NA,s2=NA,LR=NA,chi2=NA,p=NA);
  names(o)[1:4] <- c('a','b','c','s');
  # the unconstrained parameters (xpars and ypars) and the unconstrained likelihood (LL) go in the first row
  o[1,c('LL',pcols)] <- c(LL,xpars,ypars);
  # do findpars for each constraint, put resulting maxes and estimates into data.frame
  # for each LL, calculate LR, chi2, and p against the first row
  for(i in 2:nrow(o)) {
    ## Holy shit! I've been ignoring censoring in the hypothesis tests all this time!!! 130617
    oi<-opsurv(x,y,model,cx=cx,cy=cy,cons=as.logical(o[i,1:4]),par=modpars(as.numeric(o[1,pcols]),model,cno=as.logical(o[i,1:4]),pf=mean,trim=T));
    ##oi <- opsurv(x,y,model,modpars(as.numeric(o[1,pcols]),modeli=model,cno=as.logical(o[i,1:4])),cons=o[i,1:4],lb=c(1e-14,0,0,0));
    o[i,'LL'] <- oi$max; o[i,pcols] <- modpars(oi$estimate,modeli=model,cni=o[i,1:4]);
  }
  o$LR <- o$LL[1]-o$LL; o$chi2 <- 2*o$LR; o$p <- pchisq(o$chi2,o$df,lower.tail=F);
  rownames(o) <- c("",apply(o[-1,1:4],1,function(x) paste('Is',paste(c('a','b','c','s')[!as.logical(x)],collapse=','),'different?')));
  return(o);
}
