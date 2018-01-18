# DONE: allow configuring wait-times
read_ptenv<-function(logenv,...){
  oo<-list(logenv=if(is.character(logenv)) logenv<-load.ptenv(logenv) else logenv);
  nfnames <- paste0('lims.',oo$logenv$names$notfailed);
  radnames <- oo$logenv$names$radnames;
  oo$polradcols <- c(radnames,'maxrad','lims.min','lims.max');
  oo$phinames <- oo$logenv$names$phinames;
  polcols <- sapply(oo$logenv$names$radnames,c,oo$phinames);
  colnames(polcols)<-paste0('p.',colnames(polcols));
  crtsuffixes<-paste0('.X',seq_len(1+length(oo$phinames)));
  # rearrange to put 'most interesting' dimension on Z axis?!
  crtsuffixes <- c(crtsuffixes[-1],crtsuffixes[1]);
  crtcols <- sapply(oo$polradcols,paste0,crtsuffixes);
  colnames(crtcols) <- paste0('c.',colnames(crtcols));
  oo$crtcolumns <-cbind(crtcols,polcols);
  oo$maindata <- pt2df(oo$logenv); #head(oo$maindata);
  oo$points <-do.call(cbind,sapply(oo$polradcols
                                   ,function(xx) dfcrt(oo$maindata,xx),simplify=F));
  oo$subsets <- data.frame(
    sapply(setNames(sprintf('%s==1',nfnames)
                    ,paste0('s.',radnames))
           ,function(xx) eval(parse(text=xx)[[1]],envir=oo$maindata)));
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
  oo$by_cycle <- summary(group_by(oo$maindata,phicycle)
                           ,ttime=sum(time),hits=sum(lims.status==1)
                           ,hitrate=mean(lims.status==1),nn=length(time)
                           ,timeperhit=ttime/hits);
  attr(oo,'call') <- match.call();
  class(oo) <- c('summary.ptenv',class(oo));
  print(oo); invisible(oo);
  }


print.summary.ptenv<-summary.summary.ptenv<-function(x,...){
  show(list(column_groups=c(colnames(x$maindata),colnames(x$crtcolumns))
            ,named_subsets=colnames(x$subsets)
            ,last_phicycle=max(x$maindata$phicycle)
            ,points=nrow(x$maindata)
            ,rowsets=cbind(sapply(x$logenv$subsets,length))));
  #show(x$by_cycle);
  #head(x);
}

head.summary.ptenv<-function(x,...){
  head(x$maindata);
}

tail.summary.ptenv<-function(x,...){
  tail(x$maindata);
}

length.summary.ptenv <- function(x,...) nrow(x$maindata);

# x : summary.ptenv object
# subset : either a logical statement to evaluate in the context of x$maindata or
# a vector of subset names that exist in x
# select : a vector of one or more columns or column groups in x
# ... : put subset names here if desired when the subset argument is already used
# for a logical expression
# Also copying this into the [ method
subset.summary.ptenv<-function(x,subset=T,select=T,...,minphicycle=1,maxphicycle=length(x)){
  dots <- list(...);
  subset<-substitute(subset);
  oo <- cbind(x$maindata,x$points);
  lsubset <- as.list(subset);
  tfsubset <- T;
  if(is.character(subset)){
    dots<-c(dots,subset); 
    } else if(lsubset[[1]]==quote(c)||lsubset[[1]]==quote(list)) {
      dots<-c(dots,lsubset[-1]); 
      } else if(is.language(subset)||is.expression(subset)) tfsubset <- eval(eval(subset,envir=oo));
  for(ii in dots) {
    if(is.character(ii)&&ii %in% colnames(x$subsets)) tfsubset <- tfsubset & x$subsets[,ii];
    if(is.language(ii)) tfsubset <- tfsubset & eval(eval(ii,envir = oo));
  }
  tfsubset <- tfsubset & x$maindata$phicycle > minphicycle & x$maindata$phicycle <= maxphicycle;
  carts <- intersect(select,colnames(x$crtcolumns));
  select <- setdiff(select,carts);
  for(ii in carts) select <- c(x$crtcolumns[,ii],select);
  # plan to later learn more about data.table and see if the stuff that uses the 
  # output from subset.summary.ptenv would benefit from its features but sticking
  # with the clumsy approach for now and coercing back to data.frame() on output
  data.frame(oo)[tfsubset,select,drop=F];
}

`[.summary.ptenv`<-subset.summary.ptenv;

# ... : name/s of column/s for which to get ranges
range.data.frame <- function(df,...){
  if(length(dots<-list(...,quote(.)))==2) return(t(apply(df,2,range,na.rm=T)));
  rnames<-sapply(oo<-lapply(Filter(function(xx) 
    mode(xx) %in% c('character','numeric'),unlist(dots))
    ,function(zz) df[,zz,drop=F]),colnames);
  oo<-t(sapply(oo,range.default,na.rm=T)); rownames(oo)<-rnames;
  oo;
}

# model-specific point-clouds
#' Title
#'
#' @param ptenv  A summary.ptenv object
#' @param paramnames
#'
#' @return  A `list` with a `data.frame` for each member of `ptenv$logenv$names$pnfit`
#' (i.e. the verdict-returning members of the pnlst panel in `powertrip`)
#' @export
#'
#' @examples
getresults<-function(ptenv,paramnames,...){
  pnfit <- ptenv$logenv$names$pnfit;
  oo <- list();
  # ugh, argument-magic comes back to bite me
  eval(parse(text=sprintf("oo$%1$s<-d0[c('ok','s.rad.%1$s'),c('p.rad.%1$s','c.rad.%1$s','ID')];
                          rownames(oo$%1$s)<-oo$%1$s$ID;",pnfit)));
  if(!missing(paramnames)) for(ii in names(oo)) colnames(oo[[ii]])[seq_along(paramnames)]<-paramnames;
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

#' moves specified subsets to a new ptenv
env_splitoff<-function(ptenv,sets){
  targets <- intersect(sets,allsubsets<-names(ptenv$subsets));
  if(length(targets)==0) stop('None of the requested subsets were found.');
  if(length(targets)<length(sets)) warning('Some of the requested subsets were not found.');
  # names of objects to simply copy as-is
  verbatimcopy <- setdiff(names(ptenv),c('subsets','coords','allpoints','fits'));
  # new ptenv to copy to
  out <- new.env();
  # copying
  for(ii in verbatimcopy) out[[ii]] <- ptenv[[ii]];
  labscoords <- unique(names(ptenv$coords));
  labsallpts <- unique(names(ptenv$allpoints));
  # copy over the requested subsets
  out$subsets <- lapply(ptenv$subsets[targets]
                        ,function(xx) intersect(xx,c(labscoords,labsallpts)));
  # copy over the data these subsets represent
  labstocopy <- unlist(out$subsets);
  out$coords <- ptenv$coords[intersect(labstocopy,labscoords)];
  out$allpoints <- ptenv$allpoints[intersect(labstocopy,labsallpts)];
  # delete the original copies of the entries that are only referenced by the 
  # copied-over labs
  labstokeep <- unlist(ptenv$subsets[setdiff(allsubsets,targets)]);
  labstodelete <- setdiff(labstocopy,labstokeep);
  ptenv$coords[intersect(labscoords,labstodelete)] <- NULL;
  ptenv$allpoints[intersect(labscoords,labstodelete)] <- NULL;
  ptenv$subsets[targets] <- NULL;
  return(out);
}
