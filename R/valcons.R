## A function that returns a valid matrix of parameter contrasts
valcons <- function(cons,groups,model,customrows){
  ## cons = either a T/F or 1/0 vector or a matrix with groups columns and params rows
  ##        identical values in a given row mean those parameters are constrained for those groups
  ##        can be either numeric or character
  ## model = character scalar indicating model
  ## groups = character vector of group names
  ngrp<-length(groups);
  if(ngrp==1) return(matrix(
       cumsum(modparsng(rep(1,4),min='lm',mout=model)),ncol=1,
       dimnames=list(modparsng(c('a','b','c','s'),min='lm',mout=model),groups)));
  if(!is.matrix(cons)) {
    cons <- t(sapply(as.logical(cons),function(ii) if(ii) seq(1,ngrp) else rep(1,ngrp)));
  }
  nr <- nrow(cons);
  ## rnames<-c('a','b','c','s');
  cons <- cbind(t(apply(cons,1,function(ii) as.numeric(factor(ii)))));
  for(ii in 2:nr) cons[ii,]<-cons[ii,,drop=F]+max(cons[ii-1,,drop=F]);
  ## switch(model,
  ##        g={cons<-cons[1:2,]; cons<-rbind(cons,NA,NA);},
  ##        gm={if(nr<3) stop('Wrong number of parameters for GM'); cons<-cons[1:3,]; cons<-rbind(cons,NA);},
  ##        l={if(nr==4) {
  ##          cons[3,] <- NA; cons[4,]<-cons[4,]-min(cons[4,]-cons[2,])+1;
  ##          ## cons<-cons[c(1:2,4),]; cons[3,]<-cons[3,]-min(cons[3,]-cons[2,])+1;
  ##        } else if(nr(cons)<3) stop('Wrong number of parameters for L') else cons<-rbind(cons[1:2,],NA,cons[3,]);},
  ##           ## cons<-cons[1:3,]; rownames(cons)<-c(rnames,'s');}, 
  ##        lm={if(nr<4) stop('Wrong number of parameters for LM') else cons<-cons[1:4,];},
  ##        {if(missing(customrows))
  ##           stop("If an LM-family model not specified you need to manually specify rownames") else cons<-cons[1:(length(customrows)),];
  ##         rownames(cons)<-customrows;}
  ##        );
  if(missing(customrows)) rownames(cons)<-modparsng(c('a','b','c','s'),min='lm',mout=model);
  if(ncol(cons)!=ngrp) stop('Wrong number of groups for constraint matrix.') else colnames(cons)<-groups;
  cons;
}
