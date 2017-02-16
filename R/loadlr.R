`loadlr`<-
function(label,suffixes=c('g','gm','l','lm'),sep='\t',header=1,namesep='_',fpsuffix='.fp'){
	labels<-paste(label,suffixes,sep=namesep);
	for(i in labels){assign(i,read.table(i,header=header,sep=sep),envir=.GlobalEnv);}
	input<-load(paste(label,'.rdata',sep=''));
	assign(paste(label,fpsuffix,sep=''),get(input[1]),envir=.GlobalEnv);
}