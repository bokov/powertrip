#' Trying again at samplephis from the  inside out
#' 
#' To work on these individually, load `data/testdata.rda` the following hopefully
#' self-explanatory objects currently exist:
#' `test_iisims`,`test_phis`,`test_pnlst`,`test_radii`,`test_tfresp`
#' 
resp_preds <- function(tfresp,radii,glmfit,saveglm=F,power=0.8){
  if(missing(glmfit)) {glmfit <- glm(tfresp~radii,family='binomial') } else {
    glmfit <- update(glmfit) };
  out <- MASS::dose.p(glmfit,p=power);
  out <- setNames(c(out[1],attr(out,'SE')
           ,with(predict(glmfit,newdata=data.frame(radii=out[1]),type='resp',se.fit=T)
                 ,c(fit[1],se.fit[1]))
           ,glmfit$converged)
           ,c('radest','radse','respest','respse','conv'));
  if(saveglm) attr(out,'glmfit') <- glmfit;
  out;
  # fits glmfit object to tfresp ~ radius or updates optionally passed glmfit
  # returns vector that contains reverse predictor for radius, se for reverse
  # predictor, and se of prediction of glmfit when evaluated at reverse predictor
}

#' given a vector of radii and data.frame or matrix of responses with the same 
#' number of columns as the length of radii, res_preds() should iterate over each
#' column

preds_lims <- function(preds,tolse=1,tol=0.01,limit=Inf,radci=2,...){
  # takes the result of iterating res_preds() over each column of tfresps returned
  # by the panel
  # determines which panels failed on latest round 
  notfailed <- preds['conv',]==1 & preds['radest',]>0 & preds['radest',] < limit;
  # determines which panels not yet finished::
  #     sepred*setol > tol  | min < 0 | max > limits
  notdone <- preds['respse',]*tolse > tol;
  # if all panels failed, catch failure and signal accordingly
  if(!any(notfailed)) return(c(min=NA,max=NA,status=-1)); #TODO: log failure
  # if all non-failed panels converged, return final estimates with failed ones flagged
  # if there are any non-failed non-converged panels
  # returns a single max and min for the next round of radii based on those
  if(status<-!any(select <- notfailed & notdone)) select <- notfailed;
  return(c(
     min=max(min(preds['radest',select] - preds['radse',select]*radci),0)
    ,max=min(max(preds['radest',select] + preds['radse',select]*radci),limit)
    ,status=status));
}

#' Not needed, one liner: runif(nn,lims[1],lims[2]);
#lims_radii <- function(lims,phis,opts,...){
  # using limits as provided by preds_lims() and a set of phis, uniformly sample
  # radii within those limits
  # returns a vector of radii
#}

radii_res <- function(radii,phis,opts,ptsim,pnlst,...){
  # using radii as provided by lims_radii() and phis, a vector of angles
  # simulates data at each value of radii and the entire vector of phis using
  # the ptsim() function and iterates over all the functions in pnlst()
  # returns a data structure which for each value of radii has the verdict returned
  # by each function of pnlst() (T/F)
}

testenv <- new.env();
load('data/testdata.rda',envir = testenv);

test_harness <- function(list_tfresp=testenv$test_tfresp,radii=testenv$test_radii
                         ,phi=min(testenv$test_phis)
                         ,maxs=c(2,4.5),mins=c(-3.1,-1.3)
                         ,nrads=20,...){
  # The function which will plug the above modules into each other and test them
  # jointly
  # First obtain the maximum allowed radius for the current phis
  # (derived from the cartesian limits maxs and mins)
  maxrad<-pollim(phi,maxs=maxs,mins=mins);
  # then fit models on the panel verdicts (T/F), tfresp
  preds<-sapply(data.frame(do.call(rbind,list_tfresp)),resp_preds,radii=radii);
  lims <- preds_lims(preds,limit=maxrad);
  # create radii within the latest version of lims
  oldradii <- radii;
  radii <- runif(nrads,lims[1],lims[2]);
  browser();
}
