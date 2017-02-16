sp<-function(x,a,b,burnin,pw,sg,nn,parx=NULL,stf,stargs,
             maxit=Inf,stopfile='stop.txt',outfile='sp.rdata'){
  # x: collection of survival times
  # a: vector of a value boundaries
  # b: vector of b value boundaries
  # burnin: how many cycles to search for the a and b anchor points
  # pw: desired resolving power
  # nn: sample size
  # parx: fitted pars; will be found from x if it's given and parx is not
  # stf: statistical function that returns a vector starting with the p-value
  # stargs: list of arguments to stf
  # maxit: maximum number iterations
  # stopfile: if a file by this name is found to exist, this function exits
  # outfile: file to which this function will periodicall save its data
  # create the global container
  g$la<-length(a); g$lb<-length(b);
  g<-new.env(); g$x<-x; g$pw<-pw; g$sg<-sg; g$nn<-nn; g$parx<-parx;
  g$stf<-stf; g$stargs<-stargs;
  # populate the gridlist
  # the grid will have a as the horizontal axis and b as the vertical
  g$g<-list();
  for(i in 2:la) for(j in 2:lb){
    g$g[[as.character(c(a[i],b[j]))]]<-
      list(t=c(),arn=a[c(i,i-1)],brn=b[c(j,j-1)]);
  }

  # define the cpow(ab,g) function
  # sum(g$g[[ab]]$t[1,])/nrow(g$g[[ab]]$t
  
  # define the of(ab,g) function
  # x<-sample(g$x,nn,rep=nn>length(g$x));
  # pary<-g$parx;
  # pary[1]<-runif(1,g$g[[ab]]arn[1],g$g[[ab]]arn[2]);
  # pary[2]<-runif(1,g$g[[ab]]brn[1],g$g[[ab]]brn[2]);
  # y<-simsurv(nn,'lm',p=pary);
  # g$g[[ab]]$t<-rbind(g$g[[ab]]$t,do.call(g$stf,g$stargs));

  # define the pathfinder(ab,g) function
  # if any(ab==c(la,lb)) return;
  # of(ab+0:1); of(ab+1:0); of(ab+1);
  # next<-ab; next[2]<-next[2]+1;
  # 
  # at the maximum b value, compute of() at each a and append to g$g
  # set j to 1; 
  # for(i in 1:burnin){
  # moving from the left (i.e. all 1s), compute of() and append to g$g again
  # if(j+1 >= la) {j<-1; next;}
  # of(c(bmax,j));
  # if(cpow(c(bmax,j)) > g$pw & cpow(c(bmax,j+1)) <= g$pw) j<-j+1 else j<-1;
  # repeat the above at maximum a value for b
  # i<-0;
  # while(i < maxit && !file.exists(stopfile)){
  # find best a
  # pathfinder(c(a,lb),g);
  # find best b
  # pathfinder(c(la,b),g);
  # write to rfile
  # }
}
