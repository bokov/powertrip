`trypars` <-
function(d,x,mod=2,line=16,par=18:21){
	mods<-paste('ex',as.character(d[,mod]),sep='');
	out<-c();
	for(i in 1:dim(d)[1]){
		envir<-list(a=d[i,par[1]],b=d[i,par[2]],c=d[i,par[3]],s=d[i,par[4]],x=x[[d[i,line]]]);
		out[i]<-sum(eval(get(mods[i]),envir=envir));
	}
	out;
}