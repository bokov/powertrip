`simsurv`<-
function(n,type='g',p=c(2.433083e-05,0.005,3e-11,0.0015)){
  # beginning of new and faster implementation of lm and l
  # if new problems, they may start here
  if(type=='lm'){if(p[4]==0) type<-'gm' else if(p[3]==0) type<-'l' else{
    return(ceiling(rlogmake(n,a=p[1],b=p[2],c=p[3],s=p[4])))
    }
  }
  if(type=='l'){if(p[4]==0) type<-'g' else return(ceiling(rlgst(n=n,a=p[1],b=p[2],s=p[4])));}
  # end new stuff
  if(type=='gm'){if(p[3]==0) type<-'g' else return(ceiling(rmakeham(n=n,shape=p[c(1,3)],scale=1/p[2])));}
  if(type=='g'){if(p[2]==0) type<-'e' else return(ceiling(rgompertz(n=n,scale=1/p[2],shape=p[1])));}
  if(type=='w'){if(p[2]==0) type<-'e' else return(ceiling(rweibull(n=n,scale=1/p[1],shape=p[2])));}
  if(type=='e'){return(rexp(ceiling(n=n,rate=p[1])));}
  # next line added temporarily for testing out the new lm and l above
  if(type=='l.old') type<-'l'; if(type=='lm.old') type<-'lm';
  # the below is slower than the original, compiled method, but if a compiled version 
  # of the NON summed likelihood function is written, this way might be better after all
  if(type=='l2'){
    prob<-srvshp(1:1e5,a=p[1],b=p[2],s=p[4],model='l')*srvhaz(1:1e5,a=p[1],b=p[2],s=p[4]);
    prob[is.na(prob)]<-0;
    return(sample(1:1e5,n,rep=T,prob=prob));
  }
  model<-switch(type,w2=1,l=,lm=2,g2=,gm2=3,l3=4,e=5);
  if((model==2|model==4)&p[4]<6e-10){model=3;}
  if((model==3|model==1)&p[2]<2e-19){model=5;}
  if(model==5){out<-as.integer(rexp(n,as.double(p[1])));} else {
    out<-
      .C('simsurv',a=as.double(p[1]),b=as.double(p[2]),c=as.double(p[3]),s=as.double(p[4]),
         size=as.integer(n),model=as.integer(model),ans=double(n),PACKAGE="Survomatic")$ans
  }
  attr(out,'type')<-type;attr(out,'pars')<-p;
  return(out);
}
