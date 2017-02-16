`empdist` <-
function(x,y,label,cx=NULL,cy=NULL,nil=0,bnil=0,wbnil=1,pf=mean,rounds=5000,type='p',models=NULL){
	d<-int2tab(x,y); lx<-length(x); ly<-length(y); lxy<-lx+ly;
	if(is.null(models)){
		models<-c('w','wa','wb','g','ga','gb','gm','gma','gmb','gmc','l','la','lb','ls',
		          'lm','lma','lmb','lmc','lms');
	}
# 	labels<-list(main=paste(label,'.main.txt',sep=''),
# 		     unc=paste(label,'.unc.txt',sep=''),
# 		     gml=paste(label,'.gml.txt',sep=''));
	if(type=='p'){
		m<-matrix(nrow=lxy,ncol=rounds);
		for(i in 1:rounds){m[,i]<-sample(d[,2],lxy);}
		m<-cbind(d[,2],m);
		write.table(m,paste(label,'_perms.txt',sep=''),col.names=F,row.names=F,sep='\t');
		cat('\nCommencing permutations.\n');
		for(i in 1:(rounds+1)){
			findpars(d[m[,i]==0,1],d[m[,i]==1,1],nil=nil,bnil=bnil,wbnil=wbnil,
				 models=models,pf=pf,label=label,id=i-1,cx=cx,cy=cy);
		}
	}
	if(type=='r'){
		m<-matrix(nrow=lxy,ncol=rounds);
		for(i in 1:rounds){m[,i]<-c(sample(x,lx,rep=T),sample(y,ly,rep=T));}
		m<-cbind(d[,1],m);
		write.table(m,paste(label,'_resam.txt',sep=''),col.names=F,row.names=F,sep='\t');
		cat('\nCommencing resampling.\n');
		for(i in 1:(rounds+1)){
			findpars(m[d[,2]==0,i],m[d[,2]==1,i],nil=nil,bnil=bnil,wbnil=wbnil,
				 models=models,pf=pf,label=label,id=i-1,cx=cx,cy=cy);
		}
	}
# 	if(type=='b'){
# 		rsx<-matrix(nrow=lx,ncol=rounds);
# 		rsy<-matrix(nrow=ly,ncol=rounds);
# 		for(i in 1:rounds){
# 			rsx[,i]<-sample(x,lx,replace=T);
# 			rsy[,i]<-sample(y,ly,replace=T);
# 		}
# 		rsx<-cbind(x,rsx);rsy<-cbind(y,rsy);
# 		write.table(rsx,paste(label,'_xboot.txt',sep=''),col.names=F,row.names=F,sep='\t');
# 		write.table(rsy,paste(label,'_yboot.txt',sep=''),col.names=F,row.names=F,sep='\t');
# 		cat('\nCommencing bootstrap.\n');
# 		for(i in 1:(rounds+1)){
# 			findpars(rsx[,i],rsy[,i],nil=nil,bnil=bnil,models=models,pf=pf,
# 				 label=label,id=i);
# 		}
# 	}
}
