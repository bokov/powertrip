`tidyfits` <-
function(d,rm.rejected.comparisons=T,rm.rejected.models=F){
	d<-d[-grep('e|gmlm|llm',d$model),];
	if(rm.rejected.comparisons){
		if(any(d[grep('^gm$|^l$',d$model),'sig?']) ){
			if(any(filter<-grep('ga|gb',d$model))){d<-d[-filter,];}
			if(rm.rejected.models){d<-d[-grep('^g$',d$model),];}
		}
		if(any(d[grep('^lm$',d$model),'sig?']) ){
			if(any(filter<-grep('gma|gmb|gmc|la|lb|ls',d$model))){d<-d[-filter,];}
			if(rm.rejected.models){d<-d[-grep('^gm$|^l$',d$model),];}
		} else {
			if(!any(d[grep('^gm$',d$model),'sig?']) ){
				if(any(filter<-grep('gma|gmb|gmc',d$model))){d<-d[-filter,];}
				if(rm.rejected.models){d<-d[-grep('^gm$',d$model),];}
			} 
			if(!any(d[grep('^l$',d$model),'sig?']) ){
				if(any(filter<-grep('la|lb|lc',d$model))){d<-d[-filter,];}
				if(rm.rejected.models){d<-d[-grep('^l$',d$model),];}
			} 
		}
	}
	return(d);
}