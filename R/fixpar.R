`fixpar` <-
function(par,np,keep){
	lp<-length(par);cons<-rm.temp$cons;
	if(is.null(rm.temp$y)){
		if(lp==4|lp==8){return(par[1:4][keep][1:np]);}
		if(lp>np){
			warning("You gave more parameters than necessary
for a single-model fit. Only the first ",np," parameters were used.");
			return(par[1:np]);
		}
		if(lp==np){return(par);}
		if(lp==1){return(rep(par,np));}
	}
	if(lp==1){return(rep(par,2*np)[c(rep.int(T,np),cons)]);}
	if(lp==np+sum(cons)){return(par);}
	if(lp==np){return(c(par,par)[c(rep.int(T,np),cons)]);}
	if(lp==2*np){return(par[c(rep.int(T,np),cons)]);}
	if(lp==4){return(fixpar(c(par,par)[keep],np,keep));}
	if(lp==8){return(fixpar(c(par)[keep],np,keep));}
	#browser();
# 	stop("If starting parameters (par) are specified,
# they should be a numerical vector equal in length either
# to 4, or 8, or the number of parameters in your model, or
# twice the number of parameters in your model, or the
# total number of *unique* parameters in your model (i.e.
# two starting values for each unconstrained parameter and
# one starting value for each constrained one).");
}

