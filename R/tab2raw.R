tab2raw<-function(d,output=c('times','censor'),timecol=1,eventcol=2,groupcol=3,censcol=4,choosegroup=0){
	if(d[nrow(d),1]== -1) {d<-d[-nrow(d),];}
        totalevents <- d[,eventcol];
        if(!is.null(censcol)) totalevents <- totalevents + d[,censcol];
	if(!is.null(groupcol)&&length(uu <- unique(d[,groupcol]))>1){
          if(!is.null(choosegroup)){
            if(!choosegroup %in% uu) {
              ocg <- choosegroup;
              choosegroup <- uu[1];
              warning(paste(ocg,"not found in column",groupcol,"returning group",choosegroup,"instead"));
            }
            totalevents<-totalevents[d[,groupcol]==choosegroup];
            d<-d[d[,groupcol]==choosegroup,];
          } else if(is.null(choosegroup)){
            groups <- rep(d[,groupcol],totalevents);
          }
        }
        #browser();
        events <- rep(d[,timecol],totalevents);
        if(!is.null(censcol)) cens <- do.call(c,apply(d,1,function(ii) rep(1:0,c(ii[eventcol],ii[censcol]))));
        out <- data.frame(time=events);
        #browser();
        if(exists('cens')) out$censor <-  cens;
        if(exists('groups')) out$groups <- groups;
        if('times'%in%output & 'censor' %in% output ) return(out);
	if('times'%in%output) return(out$time);
        if('censor'%in%output) return(out$censor);
}
