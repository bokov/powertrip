options(error=recover);
library(Survomatic);
logenv_file <- 'pt_result.rdata';
if(file.exists(logenv_file)){
  logenv <- load.ptenv(logenv_file,savewait = 0);
  pnlst_gmcx<-logenv$state$powertrip$pnlst;
} else {
  logenv<-new.env();
  pnlst_gmcx <- list(cx=ptpnl_sr,gm=ptpnl_gm,sims=ptpnl_simsumm);
}
lrefcoords<-log(refcoords<-c(3.02495622562167e-06,0.00766970877053115,1.97042457165941e-05));
# lrelmaxs<-((lmaxs<-log(maxs <- c(1e-04, 0.0195827842383701, 0.001145822)))-lrefcoords);
# phicycle 142: expanding maxs to... 
lrelmaxs<-((lmaxs<-log(maxs <- c(0.000448943309591518, 0.0195827842383701, 0.058737528647514)))-lrefcoords);
# lrelmins<-((lmins<-log(mins <- c(2.29371429246515e-08, 0.00445842550313631, 1e-16)))-lrefcoords);
# phicycle 142: expanding to mins to... 
lrelmins<-((lmins<-log(mins <- c(1.01475976473711e-09, 0.00103798220880217, 1e-16)))-lrefcoords);
.out <- powertrip(logenv,refcoords=lrefcoords,maxs=lrelmaxs,mins=lrelmins
                  ,npoints=100,pnlst=pnlst_gmcx,ptsim=ptsim_surv,nrads=60
		  #,instance='i18010909lessrunifbias' 
		  ,instance=as.character(Sys.time(),'i%y%m%d%Ilessrunifbias')
                  ,backtrans=exp,type='gm',tol=0.1);
#savehistory(file='runscrap.R')
