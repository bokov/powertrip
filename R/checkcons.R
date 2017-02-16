`checkcons` <-
function(cons,np,keep){
	mykeep<-keep[1:np];
	cons<-as.logical(cons);
	if(length(cons)==1){return(rep.int(T,np));}
	if(length(cons)==4){return(cons[mykeep]);}
	if(length(cons)==np){return(cons);}
	stop("If a constraint (cons) is specified, is should be
  a vector of logical values equal in length either
  to 4 or to the number of parameters in your model.");
}

