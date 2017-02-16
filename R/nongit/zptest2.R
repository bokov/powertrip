`zptest2` <-
function(x1,x2,n1,n2,n,zponly=T){
	if((ln<-length(x1))!=length(x2) ){stop('The length of the two input vectors must be the same')}
	out<-.C('zptest',x1=as.double(x1),x2=as.double(x2),n1=as.double(n1),n2=as.double(n2),n=as.double(n),ln=as.integer(ln),that1=double(ln),that2=double(ln),ttil=double(ln),zp=double(ln),PACKAGE='Survomatic');
	#if(length(out)==0|is.null(out$ans)) browser();
	#cat('.'); browser();
	if(zponly){return(out$zp)} else return(data.frame(x1,out$that1,x2,out$that2,out$ttil,out$zp));
}

# The compiled version needs to operate on vectors, not atomic values; so it should do everything as vectors