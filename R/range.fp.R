range.fp<-function(...,smooth=10,what=c('time','hazard','density')){
  # Give a range for the parameter of interest for a group of fp objects.
  # Need to at some point also have the option of including fitted values
  # in the range calculation.
  what<-match.arg(what);
  objs<-list(...);

  column<-switch(what,'time' = 'time','hazard' = 'ux',
                 'density' = 'deaths');

  if(smooth>1){
    smoother<-function(x,s=smooth){
      c(rep(NA,s-1),apply(embed(x,s),1,mean,na.rm=T));
    }
  } else smoother<-function(x,s=smooth) x;
  
  d<-c();
  for(i in objs){
    d<-range(d,smoother(demog(i$x)[,column],smooth),na.rm=T);
    d<-range(d,smoother(demog(i$y)[,column],smooth),na.rm=T); 
  }
  return(d);
}
