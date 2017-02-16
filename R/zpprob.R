`zpprob` <-
function(x1,x2,n1,n2,th){
	#choose(n1,x1)*choose(n2,x2)*th^(x1+x2)*(1-th)^(n1+n2-x1-x2);
  exp(lchoose(n1,x1)+lchoose(n2,x2)+(x1+x2)*log(th)+(n1+n2-x1-x2)*log(1-th));
}

# This one returns correct answers now, but is SLOWER than the scripted version. ?!?!
`zpprob.compiled` <- function(x1,x2,n1,n2,th){
  ln <- length(x1); if(ln!=length(x2)) warning('x1 and x2 length mismatch in zpprob, will cause an error or an incorrect answer');
  .C('zpprob',x1=as.double(x1),x2=as.double(x2),n1=as.double(n1),n2=as.double(n2),ln=as.integer(ln),th=as.double(th),ans=double(ln),PACKAGE='Survomatic')$ans;
}
