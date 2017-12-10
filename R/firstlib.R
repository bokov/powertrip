.onLoad <- function(lib,pkg)
# you get busy with a couple of other things for 7 or 8 years and the whole 
# language changes on you...
#.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	packageStartupMessage(paste(c(packageDescription("Survomatic")[c('Package','Version')],' loaded\n'),collapse=' '));
  # Holy crap! Did you know that you can put browser() in here and a) R CMD INSTALL
  # will happily let you build the package and b) when you load it you immediately
  # drop into the debugger when you load the package? That's so useful and 
  # powerful I can't believe nobody prevented it from happening! Try doing that
  # in SAS or SPSS!!!!
  #browser();
	#data(ctrl,dex,modelinfo,modelpars,envir=environment());
  data(ctrl, dex, modelinfo, modelpars,package=pkg,lib.loc=lib,envir=environment());
  # Okay, I know the below is a hacky way of doing this, I was young and didn't
  # know any better back then. I promise I'll fix this once I'm done getting 
  # enough of it working to finish powertrip
	options(ctrl=ctrl);
	options(dex=dex);
	options(dmodels=dmodels);
	options(modelpars=modelpars);
}
