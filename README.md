Please note: currently active development is being done on the `test_polar` branch. Occasionally working snapshots get merged into `integration`. The `master` branch will be updated only after the current version is out of alpha and I have time to move the Survomatic code out into a separate repository and package.

For now, the usage is:

* Clone this repo.
* Check out the `test_polar` branch.
* Install the prerequisite packages listed in the DESCRIPTION file
* In the directory containing this repository (by default, `powertrip`) run this: `R CMD INSTALL R_COMPILE_PKGS=1 --byte-compile powertrip`. It will create a package called `Survomatic` (not called powertrip, I realize this is confusing, it will be fixed after the higher priority issues are dealt with).
* Run R and execute the following commands:

```R
library(compiler); setCompilerOptions(optimize=3); enableJIT(3); # optional, for speed
library(Survomatic);
# Below allows for continuation of previous saved sessions
if(file.exists('pt_result.rdata')) logenv <- load.ptenv('pt_result.rdata') else logenv <- new.env();
# these are the coordinates to set to simulate the distribution of human lifespans
lrefcoords<-log(refcoords<-c(3.02495622562167e-06,0.00766970877053115,1.97042457165941e-05));
lrelmaxs<-((lmaxs<-log(maxs <- c(1e-04, 0.0195827842383701, 0.001145822)))-lrefcoords);
lrelmins<-((lmins<-log(mins <- c(2.29371429246515e-08, 0.00445842550313631, 1e-16)))-lrefcoords);
# create a panel of test wrappers, in this case Gompertz-Makeham (ptpnl_gm) and Weibull (ptpnl_sr)
# plus a summary panel
pnlst_gmwx <- list(gm=ptpnl_gm,wx=ptpnl_sr,sims=ptpnl_simsumm);
powertrip(logenv,refcoords=lrefcoords,maxs=lrelmaxs,mins=lrelmins,npoints=100
,pnlst=pnlst_gmwx,ptsim=ptsim_surv,nrads=60,backtrans=exp,type='gm');
```
