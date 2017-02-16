modelinfo<-function(x,what=c('comp', 'rcomp', 'dfcomp','dfrcomp', 'tblcomp', 'tblrcomp', 'tblcomp2', 'deps', 'rdeps', 'tbldeps', 'tblrdeps', 'vars', 'keep', 'fullname'),extended=T){
# this function returns various information about supported models, drawn from tables this package places in options
# extended determines what comparisons should be made between models:
# the default (T) will compare each model to the next simpler one except 'lm', which will get compared to 'g', 'gm', and 'l'
# setting extended to F will only compare each model to the next simpler one (so 'lm' will get compared only to 'gm' and 'l')
what<-match.arg(what);
dmodels<-options('dmodels')[[1]];
if(extended){xt=dmodels$compare;} else {xt=dmodels$nochain;}
return(switch(what,
    comp=subset(dmodels,complex %in% x & xt)$simple,
    rcomp=subset(dmodels,simple %in% x & xt)$complex,
    dfcomp=subset(dmodels,complex %in% x & xt)$df,
    dfrcomp=subset(dmodels,simple %in% x & xt)$df,
    tblcomp=subset(dmodels,complex %in% x & xt)[,-(6:7)],
    tblrcomp=subset(dmodels,simple %in% x & xt)[,-(6:7)],
    # horrible kludge below... when whatever else uses options('dmodels')
    # is updated or removed, dmodels itself should include the exponential
    # model
    tblcomp2={o <- rbind(NA,subset(dmodels,complex %in% x & xt)[,-c(3,4,6,7)]);o[1,1:2] <- 'e'; o[1,3] <- 0; o;},
    deps=unique(subset(dmodels,complex %in% x )$simple),
    rdeps=unique(subset(dmodels,simple %in% x )$complex),
    tbldeps=subset(dmodels,complex %in% x )[,-(6:7)],
    tblrdeps=subset(dmodels,simple %in% x )[,-(6:7)],
#     comp=unlist(apply(as.matrix(options('modeldeps')[[1]][,x,1]),2,function(y) names(y)[y])),
#     deps=rownames(options('modeldeps')[[1]])[apply(as.matrix(options('modeldeps')[[1]][,x,2]),1,function(y) any(y))],
#     rdeps=rownames(options('modeldeps')[[1]])[as.matrix(options('modeldeps')[[1]][x,,2])],
    vars=rownames(options('modelpars')[[1]])[(options('modelpars')[[1]][,c(x[1])])],
    keep=(1:8)[(options('modelpars')[[1]][,c(x[1])])],
    fullname={
	x<-gsub('^e','Exponential',x);x<-gsub('^w','Weibull',x);x<-gsub('^gm','Gompertz-Makeham',x);
	x<-gsub('^lm','Logistic-Makeham',x);x<-gsub('^g','Gompertz',x);x<-gsub('^l','Logistic',x); x;
    }
));
}
