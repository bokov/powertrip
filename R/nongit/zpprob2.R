`zpprob2` <-
function(x1,x2,n1,n2,th){
	if((ln<-length(x1))!=length(x2)&length(x2)>1){
		stop('The length of the two input vectors must be the same')
	}
	out<-.C('zpprob',x1=as.double(x1),x2=as.double(x2),n1=as.double(n1),n2=as.double(n2),ln=as.integer(ln),th=as.double(th),ans=double(ln),PACKAGE='Survomatic');
	return(out$ans);
}