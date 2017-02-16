`maxmodel`<-function(x,y){
    out<-c();
    if('w' %in% x & 'w' %in% y){out<-'w';} else {if(xor('w' %in% x, 'w' %in% y)){
	warning('One of the model lists doesn\'t contain the Weibull model while the other does.');
    }}
    xy<-grep('w',unique(c(x,y)),value=T,invert=T); modtemp<-options('dmodels')[[1]][,c('complex','dfc')]; 
    if(('l' %in% xy & 'gm' %in% xy)|'lm' %in% xy){return(c(out,'lm'));}
    modtemp<-unique(subset(modtemp,complex %in% xy));
    return(c(out,modtemp[which.max(modtemp$dfc),'complex']));
}
