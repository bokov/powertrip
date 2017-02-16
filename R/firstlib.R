.First.lib <- function(lib,pkg)
{
	library.dynam("Survomatic",pkg,lib.loc=lib);
	packageStartupMessage(paste(c(packageDescription("Survomatic")[c('Package','Version')],' loaded\n'),collapse=' '));
	data(ctrl,dex,modelinfo,modelpars,envir=environment());
	options(ctrl=ctrl);
	options(dex=dex);
	options(dmodels=dmodels);
	options(modelpars=modelpars);
}
