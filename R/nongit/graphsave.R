`graphsave` <-
function(path,xynames,suffix,prompt){
	if(prompt<2){return();}
	fname<-paste(path,paste(xynames,collapse=""),suffix,sep="");
	message<-paste("Do you wish to save this graph as",fname,"?\n");
	choices<-c("Save under suggested name.",
		   "Give a different name (point and click).",
		   "Type in a name from the console.",
		   "Don\'t save.")
	cat(message);
	input<-menu(choices,graphics=T,title=message);
	switch(input,
		dev.copy2eps(file=fname),
		{fname<-tclvalue(tkgetSaveFile(initialdir=path,initialfile=paste(xynames,collapse=""),
			defaultextension=suffix));dev.copy2eps(file=fname);},
		{res<-0;class(res)<-'try-error';
		 while(class(res)=='try-error'){
			fname<-readline("Please type the path and name to which you want to save the file: ");
			print(fname); res<-try(dev.copy2eps(file=fname));
			if(class(res)=='try-error'){cat('Oops. Bad file name or path, try again please.\n')}}
		},
		cat("\nContinuing without saving.")
	)

}

