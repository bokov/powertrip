# The purpose of this function is to obtain an empirical null distribution
# for any test statistic comparing two groups. By default it is the t.
empnull <- function(x,y,test=function(a,b) return(t.test(a,b)$statistic),
                    nboot=1000,keepperm=F){
  if(is.vector(x)&is.vector(y)){
    nx <- length(x); ny <- length(y); xy <- c(x,y);
  } else if((is.matrix(x)|is.data.frame(x))&(is.matrix(y)|is.data.frame(y))){
    nx <- nrow(x); ny <- nrow(y); xy <- rbind(x,y);
  } else stop('x and y must both be either vectors, matrices, or data frames');
  nxy <- nx+ny;
  mboot <- replicate(nboot,sample(1:nxy,nxy));
  if(is.vector(x)&is.vector(y)) dboot <- apply(mboot,2,function(a) test(xy[a[1:nx]],xy[a[nx+(1:ny)]]))
  else dboot <- apply(mboot,2,function(a) test(xy[a[1:nx],],xy[a[nx+(1:ny)],]));
  t <- test(x,y); p <- {if(all(abs(dboot)==abs(t))) 1 else mean(abs(dboot)>abs(t))};
#  list(dboot=dboot,test=t,p=min(mean(dboot>t),mean(dboot<t)));
  list(dboot={if(keepperm) dboot else NULL},test=t,p=mean(abs(dboot)>abs(t)));
}
