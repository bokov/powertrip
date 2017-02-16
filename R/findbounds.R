findbounds<-function(n=n,cx=cx,xpars,mypars,xmodel=xmodel,pow=pow,sig=sig,
                     mbounds,tests,ps,sim,fun=fun,flipstop=3,mindiff=1e-10){
  # mbounds is a matrix of afac, bfac, abounds, bbounds (in columns)
  # on each round, modify mypars[i] if the ith row in mbounds is not NA
  ###### the below needs to be modified to take bounds from a matrix
  mypars.good<-mypars; mbounds.good<-mbounds; sim.good<-sim; ps.good<-ps;
  mbounds.prev<-mbounds; mypars.prev<-mypars;
  use<-apply(mbounds,2,function(x) all(!is.na(x)));
  mbounds.prev[2,]<-mbounds.prev[2,]*0+Inf;
  mbounds.prev[3,]<-mbounds.prev[3,]*0;
  print(abs(diff(mbounds[-1,]))); print(abs(diff(mbounds.prev[-1,])));
  flip<-0; nflips<-0;
  xsim<-matrix(simsurv(n*tests[1],xmodel,xpars),nrow=n);
  #browser();
  while(nflips<flipstop&max(c(abs(diff(mbounds[-1,])),abs(diff(mbounds[-1,]))),na.rm=T)>mindiff){
    # create an 'n' by 'tests' matrix of simulated data using mypars
    sim<-matrix(simsurv(n*tests[1],xmodel,mypars),nrow=n);
    ps<-numeric(tests[1]);
    # apply the fun to each column, returning the pval
    for(i in 1:tests[1]) ps[i]<-fun(xsim[,i],sim[,i]);
    #ps<-apply(sim,2,fun); 
    mbounds.prev<-mbounds; mypars.prev<-mypars;
    cat('Min Diff:',mbounds[2,use],'  Max Diff:',mbounds[3,use],'\n');
    cat('Power: ',sum(ps<sig)/tests[1],'  Flips:',nflips,'\n');
    # if the fraction of the pvals below 'sig' is greater than 'pow'
    if(sum(ps<sig)/tests[1] >= pow){
      mypars.good<-mypars; mbounds.good<-mbounds; sim.good<-sim; ps.good<-ps;
      mbounds[3,use]<-mypars[1:2][use];
      #mbounds[2,use]<-mypars[1:2][use];
      if(flip>0) {nflips<-nflips+1;}; flip<- -1;
    } else {
      mbounds[2,use]<-mypars[1:2][use];
      #mbounds[3,use]<-mypars[1:2][use];
      if(flip<0) {nflips<-nflips+1;}; flip<- 1;
    }
    mypars[1:2][use]<-apply(mbounds[-1,],2,mean)[use];
  }
  if(length(tests)>1){
    return(findbounds(n,cx,xpars,mypars.good,xmodel,pow,sig,mbounds.good,
                      tests[-1],ps,sim,fun,flipstop,mindiff));
  } else return(list(mypars=mypars.good,mbounds=mbounds.good,ps=ps,sim=sim.good));
}
