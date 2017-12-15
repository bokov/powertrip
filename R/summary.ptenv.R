read_ptenv<-function(logenv,...){
  oo<-list(logenv=if(is.character(logenv)) logenv<-load.ptenv(logenv) else logenv);
  nfnames <- paste0('lims.',logenv$names$notfailed);
  radnames <- logenv$names$radnames;
  oo$polradcols <- c(radnames,'maxrad','lims.min','lims.max');
  oo$phinames <- logenv$names$phinames;
  oo$crtcolumns <- sapply(oo$polradcols,paste0,c('.X2','.X3','.X1'));
  colnames(oo$crtcolumns) <- paste0('c.',colnames(oo$crtcolumns));
  oo$maindata <- pt2df(logenv); #head(oo$maindata);
  oo$points <-do.call(cbind,sapply(oo$polradcols
                                   ,function(xx) dfcrt(oo$maindata,xx),simplify=F));
  oo$subsets <- data.frame(sapply(setNames(sprintf('%s==1',nfnames),paste0('s.',radnames)),function(xx) eval(parse(text=xx)[[1]],envir=oo$maindata)));
  oo$subsets$notok<-with(oo$maindata,lims.status!=1);
  oo$subsets$ok<-with(oo$maindata,lims.status==1);
  #with(dat4,addmargins(table(lims.notfailed.cx,lims.notfailed.gm)))/nrow(dat4);
  #oo$points$extents <- dfcrt(transform(oo$maindata
  #                                     ,rad=pmax(ifelse())))
  #pt4all<-dfcrt(transform(dat4,rad=pmax(ifelse(lims.notfailed.cx==1,rad.cx,0),ifelse(lims.notfailed.gm==1,rad.gm,0))),'rad',subset=lims.status==1);
  #oo$points$minpts <- dfcrt(oo$maindata,'lims.min',subset=lims.status==1&phicycle>1);
  #oo$points$maxpts <- dfcrt(oo$maindata,'lims.max',subset=lims.status==1&phicycle>1);
  #ptcx4.0 <- dfcrt(dat4,'rad.cx',subset=lims.notfailed.cx==1);
  #ptgm4.0 <- dfcrt(dat4,'rad.gm',subset=lims.notfailed.gm==1);
  #ptgoodbox4 <- dfcrt(dat4,'maxrad',subset=lims.status==1);
  #ptbadbox4 <- dfcrt(dat4,'maxrad',subset=lims.status!=1);
  oo$by_cycle <- summarize(group_by(dat4,phicycle)
                           ,ttime=sum(time),hits=sum(lims.status==1)
                           ,hitrate=mean(lims.status==1),nn=length(time)
                           ,timeperhit=ttime/hits);
  attr(oo,'call') <- match.call();
  class(oo) <- c('summary.ptenv',class(oo));
  oo;
  }


print.summary.ptenv<-summary.summary.ptenv<-function(x,...){
  show(list(column_groups=c(colnames(x$maindata),colnames(x$crtcolumns))
            ,named_subsets=colnames(x$subsets)
            ,last_phicycle=max(x$maindata$phicycle)
            ,points=nrow(x$maindata)));
  #show(x$by_cycle);
  #head(x);
}

head.summary.ptenv<-function(x,...){
  head(x$maindata);
}

length.summary.ptenv <- function(x,...) nrow(x$maindata);

# x : summary.ptenv object
# subset : either a logical statement to evaluate in the context of x$maindata or
# a vector of subset names that exist in x
# select : a vector of one or more columns or column groups in x
# ... : put subset names here if desired when the subset argument is already used
# for a logical expression
# Also copying this into the [ method
`[.summary.ptenv`<-subset.summary.ptenv<-function(x,subset=T,select=T,...,minphicycle=1,maxphicycle=length(x)){
  dots <- list(...);
  subset<-substitute(subset);
  oo <- cbind(x$maindata,x$points);
  if(is.character(subset)||subset[[1]]==quote(c)||subset[[1]]==quote(list)) {
    dots<-c(dots,as.list(subset[-1])); tfsubset <-T} else {
    tfsubset <- eval(subset,envir=oo)};
  for(ii in dots) if(is.character(ii)&&ii %in% colnames(x$subsets)) tfsubset <- tfsubset & x$subsets[,ii];
  tfsubset <- tfsubset & x$maindata$phicycle > minphicycle & x$maindata$phicycle <= maxphicycle;
  carts <- intersect(select,colnames(x$crtcolumns));
  select <- setdiff(select,carts);
  for(ii in carts) select <- c(x$crtcolumns[,ii],select);
  oo[tfsubset,select,drop=F];
}

# ... : name/s of column/s for which to get ranges
range.data.frame <- function(df,...){
  if(length(dots<-list(...,quote(.)))==2) return(t(apply(df,2,range,na.rm=T)));
  rnames<-sapply(oo<-lapply(Filter(function(xx) 
    mode(xx) %in% c('character','numeric'),unlist(dots))
    ,function(zz) df[,zz,drop=F]),colnames);
  oo<-t(sapply(oo,range.default,na.rm=T)); rownames(oo)<-rnames;
  oo;
}

# lazy 3d plotting
# which: name of a radial variable (maxrad, lims.max, lims.min, or one of the model radii)
# subset : something that can be passed to the subset argument of subset.summary.ptenv()
# ... : passed to plot3d
# polar : whether to plot the polar or the cartesian version (default)
pt3d<-function(ptenv,which,subset,...,polar=F){
  which<-which[1];
  which1 <- if(which %in% ptenv$polradcols) which else {
    if((which0 <- paste0('rad.',which))%in%ptenv$polradcols)  which0 else {
      if((which0 <- paste0('lims.',which))%in%ptenv$polradcols) which0
    }};
  select <- if(polar) c(which1,ptenv$phinames) else paste0('c.',which1);
  if(missing(subset)){
    subset <- if((sub0<-paste0('s.',which1))%in%names(ptenv$subsets)) sub0 else 'ok';
  }
  rgl::plot3d(ptenv[subset,select],...);
  };
