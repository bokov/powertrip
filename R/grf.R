# This is the old, interpreted version, which will be removed once the new, compiled version is verified
`grf2` <-
function(par,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,model='g'){
 par1<-par2<-rep.int(NA,4); par1[keep]<-par2[keep]<-par[1:np]; cons<-cons[keep];
 par2[keep][cons>0]<-par[(np+1):length(par)];
 env1<-list(a=par1[1],b=par1[2],c=par1[3],s=par1[4],x=x); 
 env2<-list(a=par2[1],b=par2[2],c=par2[3],s=par2[4],x=y);
 out1<-out2<-c();
 dex<-getOption('dex');
 for(i in dex[[model]]){
 	out1<-c(out1,sum(eval(i,envir=env1)));out2<-c(out2,sum(eval(i,envir=env2)));
 }
 out1<-c(out1,out1); out2<-c(out2,out2);
 id1<-c(rep.int(1,np),1-cons); id2<-c(1-cons,rep.int(1,np));
 out<-(id1*out1)+(id2*out2); 
 if(rm.temp$debug){cat("gradient:");print(out);}
 out[c(1:np,(np+1:np)[as.logical(cons)])];
}

`grf` <-
function(par,cons=rep.int(1,4),x,y=NULL,keep=1:2,np=2,model="gg",nx,ny,cx,cy,tlog=F){
 par1<-par2<-rep.int(0,4); par1[keep]<-par2[keep]<-par[1:np]; cons<-cons[keep];
 par2[keep][cons>0]<-par[(np+1):length(par)];
 if(tlog){par1<-exp(par1);par2<-exp(par2);}
 out1<-.C(model,a=as.double(par1[1]),b=as.double(par1[2]),c=as.double(par1[3]),s=as.double(par1[4]),x=as.integer(x),size=as.integer(nx),censor=as.integer(cx),ans=double(4),PACKAGE="Survomatic")$ans[1:np];
 if(ny>0){out2<-.C(model,a=as.double(par2[1]),b=as.double(par2[2]),c=as.double(par2[3]),s=as.double(par2[4]),x=as.integer(y),size=as.integer(ny),censor=as.integer(cy),ans=double(4),PACKAGE="Survomatic")$ans[1:np];} else {out2<-rep(0,np);}
 out1<-c(out1,out1); out2<-c(out2,out2);
 id1<-c(rep.int(1,np),1-cons); id2<-c(1-cons,rep.int(1,np));
 out<-(id1*out1)+(id2*out2); 
 if(rm.temp$debug){cat("gradient:");print(out);}
 out[c(1:np,(np+1:np)[as.logical(cons)])];
}

