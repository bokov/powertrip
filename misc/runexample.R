options(error=recover);
library(Survomatic);
logenv_file <- 'pt_result.rdata';
if(file.exists(logenv_file)){
   logenv <- load.ptenv(logenv_file,savewait = 0);
   #pnlst_gmcx<-logenv$state$powertrip$pnlst;
 } else {
   logenv<-new.env();
   #pnlst_gmcx <- list(cx=ptpnl_sr,gm=ptpnl_gm,sims=ptpnl_simsumm);
 }
#tol<-0.1;
lrefcoords<-log(refcoords<-c(3.02495622562167e-06,0.00766970877053115,1.97042457165941e-05));
# lrelmaxs<-((lmaxs<-log(maxs <- c(1e-04, 0.0195827842383701, 0.001145822)))-lrefcoords);
# phicycle 142: expanding maxs to... 
lrelmaxs<-((lmaxs<-log(maxs <- c(0.000448943309591518, 0.0195827842383701, 0.058737528647514)))-lrefcoords);
# lrelmins<-((lmins<-log(mins <- c(2.29371429246515e-08, 0.00445842550313631, 1e-16)))-lrefcoords);
# phicycle 142: expanding to mins to.
#lrelmins<-((lmins<-log(mins <- c(1.01475976473711e-09, 0.00103798220880217, 1e-16)))-lrefcoords);
# might need to tighten up lower bound on IMR so make_phis() will pay more attention to the
# extremes of the RoA range for the Weibull model (cx) (starting with IMRminup instance)
#lrelmins<-((lmins<-log(mins <- c(1.42888814448938e-06, 0.00103798220880217, 1e-16)))-lrefcoords);
# oops, that was wrong, should have been the RoA param
lrelmins<-((lmins<-log(mins <- c(1.01475976473711e-09, 0.00362291389246332, 1e-16)))-lrefcoords);
#' Exploring 'left' edge of surface
# leftinterestingmaxs<-c(-2.5,0.5,1);
# leftinterestingmins<-c(-8,0.1,-26);
# simulating over one specific zone
#' Exploring 'right' edge of surface
#rzmaxs <- c(5,-0.25,2);
#rzmins <- c(1,-0.75,-15);
#fzmaxs <- c(5,0.9373722,8);
#fzmins <- c(-8,-.75,.1);
#' ## Fresh start -- no training data
#fsmaxs <- c(4,0.9373722,8);
#fsmins <- c(-8,-.75,-20);
#' ## Cold start with off-center boundaries (lower half of front curvy part)
#lfzmaxs <- c(4,-0.005,6);
#lfzmins <- c(-2,-0.75,1);
nnrefs <- c(setNames(lrefcoords,c('imr','roa','eh')),nn=0);
nnmaxs <- c(imr=6.609438,roa=1.630519,eh=10,nn=6.39692965521615);
nnmins <- c(imr=-14.907755,roa=-2.359438, eh=-26, nn=2.30258509299405);
nnptsim <- ptsim_srvn;
#logenv<-new.env();
pnlst_fresh <- list(sr=ptpnl_sr,gm=ptpnl_gm,diff=ptpnl_diff,sims=ptpnl_simsumm);
tol<-0.05;
.out <- powertrip(logenv #,refcoords=lrefcoords
                  ,refcoords = nnrefs
                  #,maxs=leftinterestingmaxs,mins=leftinterestingmins
                  #,maxs=rzmaxs,mins=rzmins
                  #,maxs=fzmaxs,mins=fzmins
                  #,maxs=fsmaxs,mins=fsmins
                  #,maxs=lfzmaxs,mins=lfzmins
                  ,maxs=nnmaxs,mins=nnmins
                  ,npoints=100
                  #,pnlst=pnlst_gmcx
                  ,pnlst=pnlst_fresh
                  #,ptsim=ptsim_surv
                  ,ptsim=nnptsim
                  ,nrads=60
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Mleftzone')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Mlzn_aftersplitoff')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Mrzn_aftersplitoff')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Mfzn_aftersplitoff')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Mfresh')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%Moffctr')
                  ,instance=as.character(Sys.time(),'i%y%m%d%I%Mnn')
                  ,backtrans=exp,type='gm',tol=tol);

# .out <- powertrip(logenv,refcoords=lrefcoords,maxs=lrelmaxs,mins=lrelmins
#                   ,npoints=100,pnlst=pnlst_gmcx,ptsim=ptsim_surv,nrads=60
		  #,instance='i18010909lessrunifbias' 
		  #,instance=as.character(Sys.time(),'i%y%m%d%Ilessrunifbias')
		  #,instance=as.character(Sys.time(),'i%y%m%d%I%MIMRminup')
	     	  #, instance=as.character(Sys.time(),'i%y%m%d%I%MIMRRewnkeeps')
                  #,instance=as.character(Sys.time(),'i%y%m%d%I%MRoAnewnkeeps')
		  # ,instance=as.character(Sys.time(),'i%y%m%d%I%Mnkeep127')
		  #                   ,backtrans=exp,type='gm',tol=0.1);
#savehistory(file='runscrap.R')
