`survgui`<-
function(xtab=NULL,ytab=NULL,dat=NULL,xynames=c(deparse(substitute(xtab)),deparse(substitute(ytab))),surv_percentalive=TRUE,surv_mark=NA,surv_pch=23:24,surv_col=c('black','darkred'),logrank_perms=2000,mperms=20,qsigf=function(s){return(s)},pcx=6,qregylim=c(-3,3)){
xynames<-as.character(xynames);
plotmortlabels<-c('Survival, model fit', 'Survival, both' ,'Hazard, actual data' ,'Hazard, model fit' ,'Hazard, both' ,'Density, actual data' ,'Density, model fit' ,'Density, both');
survlabel<-'Survival, actual data';
plotqrlabel<-'Quantile regression';
dosummlabel<-'Run or update summary report'<-prsummlabel<-savesummlabel<-'Summary';
doscorelabel<-'Score test (time consuming)';prscorelabel<-savescorelabel<-'Score test';
doqrlabel<-prqrlabel<-saveqrlabel<-'Quantile regression';
domortlabel<-'Fit mortality models and compare parameters between the two groups (time consuming)'; prmortlabel<-savemortlabel<-'Mortality models';
doalllabel<-'Run all tests (time consuming)';pralllabel<-'All results';savealllabel<-'All results (as .Rdata file)';
saveimglabel<-'Current image (as PDF)';
# make sure the environment can support this GUI
require(tcltk);
tclRequire('Tktable');
tcl("set","surv_percentalive",as.numeric(surv_percentalive));
# surv_mark is the mark parameter of plot.survfit; if it's NA then censoring will not be printed
# otherwise the corresponding pch character will be used
tcl("set","surv_mark",ifelse(is.na(surv_mark),'',surv_mark));
tcl("set","surv_pch1",ifelse(is.na(surv_pch[1]),'',surv_pch[1]));
tcl("set","surv_pch2",ifelse(is.na(surv_pch[2]),'',surv_pch[2]));
tcl("set","qregylim1",qregylim[1]);
tcl("set","qregylim2",qregylim[2]);
tcl("set","surv_col1",surv_col[1]);tcl("set","surv_col2",surv_col[2]);
tcl("set","logrank_perms",logrank_perms);
tcl("set","pcx",pcx);

# declare a bunch of tcl variables for later use
a<-tktoplevel(); n<-20;
if(is.vector(xtab)){xtab<-data.frame(time=xtab,censor=1);}
if(is.vector(ytab)){ytab<-data.frame(time=ytab,censor=1);}
frame1<-tkframe(a,pady=5,padx=5);
frame2<-tkframe(a);frame3<-tkframe(a);frame4<-tkframe(frame2);
ar10<-tkfont.create(family='arial',size=10); ar8<-tkfont.create(family='arial',size=8);
outputfont<-tkfont.create(family='courier',size=8);
qincr<-c(.01,.05,.1,.2); tcl('set','qrse','boot'); penv<-environment(); cursor<-1;
digits=4; data(references,envir=penv); myrefs<-c();
oldwidth<-options('width'); options(width=280);

# core functions
cleanUpObjectsFunc<-function(){
    for(i in c('summ.out','lrnk.out','scor.out','qreg.out','qreg.sum','mort.out')){
	if(exists(i,envir=penv)){rm(list=i,envir=penv);}
    }
#     for(i in c('plSurvBut', 'svSummBut', 'svScorBut', 'svQregBut', 'plQregBut', 'svMortBut', 'svAllBut', 'prSummBut', 'prScorBut', 'prQregBut', 'prMortBut', 'prAllBut')){
# 	tkconfigure(get(i),state='disabled');
#     }
    for(i in c(survlabel,plotmortlabels,plotqrlabel)){tkentryconfigure(plotMenu,i,state='disabled');} 
    tkentryconfigure(doMenu,dosummlabel,state='disabled');
}

loadWMFunc<-function(filename=NULL){
    if(is.null(filename)){
	filename<-tclvalue(tkgetOpenFile(initialdir=tclvalue('path'), filetypes="{{WinModest} {.dat .txt}}"));
    }
    if(filename!=""){
	d<-read.delim(filename,sep='\t',header=F); 
	if(d[nrow(d),1]== -1) {d<-d[-nrow(d),];}
	if(ncol(d)==4 & nrow(d) > 4 & all(sort(unique(d[,3]))==1:2)){
	    cleanUpObjectsFunc();
	    assign('xtab',tab2raw(d,output=c('times','censor'),choosegroup=1),env=penv);
	    assign('ytab',tab2raw(d,output=c('times','censor'),choosegroup=2),env=penv);
	    assign('n',max(dim(xtab)[1],dim(ytab)[1],n),env=penv);
	    #tkconfigure(plSurvBut,state='active');
	    tkentryconfigure(plotMenu,survlabel,state='normal');
	} else warning('File is not in WinModest format. Please try again.');
    }
}
loadFunc<-function(filename=NULL){ 
    tkconfigure(a,cursor='watch');
    if(is.null(filename)){
	filename<-tclvalue(tkgetOpenFile(initialdir=tclvalue('path'),filetypes="{{R Data} {.rdata .Rdata}}"));
    }
    if(filename!=""){
	cleanUpObjectsFunc();
	load(filename); objcount<-0;
	if(dev.cur()>1){dev.off(dev.cur());}
	# this used to test for existance but since xtab and ytab
	# have been added as parameters defaulting to null, this
	# now tests for nullness
	tcl("set","surv_percentalive",as.numeric(surv_percentalive));
	tcl("set","surv_mark",ifelse(is.na(surv_mark),'',surv_mark));
	tcl("set","surv_pch1",ifelse(is.na(surv_pch[1]),'',surv_pch[1]));
	tcl("set","surv_pch2",ifelse(is.na(surv_pch[2]),'',surv_pch[2]));
	tcl("set","qregylim1",qregylim[1]);
	tcl("set","qregylim2",qregylim[2]);
	tcl("set","surv_col1",surv_col[1]);tcl("set","surv_col2",surv_col[2]);
	tcl("set","logrank_perms",logrank_perms);
	if(exists('xynames')){
	    tcl("set","name1",xynames[1]);tcl("set","name2",xynames[2]);
	    objcount<-objcount+1;
	}
	if(exists('summ.out')){
	    #tkconfigure(svSummBut,state='active',foreground='black', activeforeground='black');
	    #tkconfigure(svAllBut,state='active',foreground='black', activeforeground='black');
	    #tkconfigure(prSummBut,state='active',foreground='black', activeforeground='black');
	    #tkconfigure(prAllBut,state='active',foreground='black', activeforeground='black');
	    tkentryconfigure(prMenu,prsummlabel,state='active');
	    tkentryconfigure(prMenu,pralllabel,state='active');
	    tkentryconfigure(fileSaveMenu,savesummlabel,state='active');
	    tkentryconfigure(fileSaveMenu,savealllabel,state='active');
	    assign('summ.out',summ.out,env=penv)
	    objcount<-objcount+1;
	}
	if(exists('lrnk.out')){assign('lrnk.out',lrnk.out,env=penv);objcount<-objcount+1;}
	if(exists('scor.out')){
# 	    tkconfigure(svScorBut,state='active',foreground='black', activeforeground='black');
# # 	    tkconfigure(svAllBut,state='active',foreground='black', activeforeground='black');
# 	    tkconfigure(prScorBut,state='active',foreground='black', activeforeground='black');
# # 	    tkconfigure(prAllBut,state='active',foreground='black', activeforeground='black');
	    tkentryconfigure(prMenu,prscorelabel,state='active');
	    tkentryconfigure(prMenu,pralllabel,state='active');
	    tkentryconfigure(fileSaveMenu,savescorelabel,state='active');
	    tkentryconfigure(fileSaveMenu,savealllabel,state='active');
	    assign('scor.out',scor.out,env=penv);
	    objcount<-objcount+1;
	}
	if(exists('qreg.out')){
# 	    tkconfigure(svQregBut,state='active',foreground='black', activeforeground='black');
# 	    tkconfigure(svAllBut,state='active',foreground='black', activeforeground='black');
# 	    tkconfigure(prQregBut,state='active',foreground='black', activeforeground='black');
# 	    tkconfigure(prAllBut,state='active',foreground='black', activeforeground='black');
	    tkentryconfigure(prMenu,prqrlabel,state='active');
	    tkentryconfigure(prMenu,pralllabel,state='active');
	    tkentryconfigure(fileSaveMenu,saveqrlabel,state='active');
	    tkentryconfigure(fileSaveMenu,savealllabel,state='active');
	    assign('qreg.out',qreg.out,env=penv);
	    objcount<-objcount+1;
	}
	if(exists('qreg.sum')){
	    #tkconfigure(plQregBut,state='active',foreground='black', activeforeground='black');
	    tkentryconfigure(plotMenu,plotqrlabel,state='active');
	    assign('qreg.sum',qreg.sum,env=penv)
	    objcount<-objcount+1;
	}
	if(exists('mort.out')){
	    tkentryconfigure(prMenu,prmortlabel,state='active');
	    tkentryconfigure(prMenu,pralllabel,state='active');
	    tkentryconfigure(fileSaveMenu,savemortlabel,state='active');
	    tkentryconfigure(fileSaveMenu,savealllabel,state='active');
	    for(i in plotmortlabels) tkentryconfigure(plotMenu,i,state='normal');
	    assign('mort.out',mort.out,env=penv);
	    objcount<-objcount+1;
	}
    }
    tkfocus(a);
    tkconfigure(a,cursor='');
}

quitFunc<-function(){options(width=oldwidth[[1]]);tkdestroy(a);}

doGeneralFunc<-function(){
	# the [,1:2] subscript might be extraneous; remove it when nothing has broken in a while
# 	assign('xtab',na.exclude(tcl2nArray(vard1,1:n,0:1))[,1:2],env=penv);
# 	assign('ytab',na.exclude(tcl2nArray(vard2,1:n,0:1))[,1:2],env=penv);
	assign('x',xtab[xtab[,2]==1,1],env=penv);
	assign('y',ytab[ytab[,2]==1,1],env=penv);
	assign('xynames',c(tclvalue('name1'),tclvalue('name2')),env=penv);
	assign('tmain',paste(xynames,collapse=' vs. '),env=penv);
	assign('xy',c(x,y),env=penv);
	assign('lx',length(x),env=penv);assign('ly',length(y),env=penv);
	assign('group',factor(rep(xynames,c(lx,ly)),levels=xynames),env=penv);
	assign('surv_percentalive',as.numeric(tclvalue("surv_percentalive")),env=penv);
	assign('surv_mark',as.numeric(tclvalue("surv_mark")),env=penv);
	assign('surv_pch',as.numeric(c(tclvalue('surv_pch1'),tclvalue('surv_pch2'))),env=penv);
	assign('qregylim',as.numeric(c(tclvalue('qregylim1'),tclvalue('qregylim2'))),env=penv);		
	assign('surv_col',c(tclvalue('surv_col1'),tclvalue('surv_col2')),env=penv);
	assign('logrank_perms',as.numeric(tclvalue("logrank_perms")),env=penv);
}

doSummFunc<-function(){
    # xtab clean
    tkconfigure(a,cursor='watch');
    doGeneralFunc(); 
    sdx<-sd(xtab[,1]); sdy<-sd(ytab[,1]); mx<-mean(xtab[,1]); my<-mean(ytab[,1]); 
    qx<-round(quantile(xtab[,1],c(.5,.9)),0); qy<-round(quantile(ytab[,1],c(.5,.9)),0);
    qcix<-qci(xtab[,1],c(.5,.9)); qciy<-qci(ytab[,1],c(.5,.9));
    #xyt<-t.test(log(x[x>0]),log(y[y>0]));
    xytab<-rbind(cbind(xtab,group=1),cbind(ytab,group=2));
    xylr<-surv2.logrank(Surv(xytab[,1],event=xytab[,2]),xytab[,3],nsim=logrank_perms);
    summ.out<-data.frame(NA,lx,NA,NA,ly,NA,NA);
    colnames(summ.out)<-c(paste(xynames[1],'lci'),xynames[1],paste(xynames[1],'uci'),
			    paste(xynames[2],'lci'),xynames[2],paste(xynames[2],'uci'),'p');
    rownames(summ.out)<-'n';
    summ.out<-rbind(summ.out,mean=c(mx-sdx,mx,mx+sdx,my-sdy,my,my+sdy,NA));
    summ.out<-rbind(summ.out,median=c(qcix[1,1],qx[1],qcix[1,2],qciy[1,1],qy[1],qciy[1,2],
		    NA)); #if qreg.out exists, get p from there
    summ.out<-rbind(summ.out,`90th percentile`=c(qcix[2,1],qx[2],qcix[2,2],
		    qciy[2,1],qy[2],qciy[2,2],
		    NA)); #if qreg.out exists, get p from there
    if(exists('scor.out')){
	    summ.out['median','p']<-scor.out['0.5','p'];
	    summ.out['90th percentile','p']<-scor.out['0.9','p'];
    }
    assign('summ.out',summ.out,env=penv); assign('lrnk.out',xylr,env=penv);
    prSummFunc();
#     tkconfigure(svSummBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(prSummBut,state='normal',foreground='black',activeforeground='black');
#     tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
    tkentryconfigure(prMenu,prsummlabel,state='normal');
    tkentryconfigure(prMenu,pralllabel,state='normal');
    tkentryconfigure(fileSaveMenu,savesummlabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(fileSaveMenu,savealllabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(doMenu,dosummlabel,foreground='black',activeforeground='black');
    #tkconfigure(doSummBut,foreground='black',activeforeground='black');
    tkconfigure(a,cursor='');
}
doScorFunc<-function(){
    # xtab clean
    tkconfigure(a,cursor='watch');
    doGeneralFunc();
    assign('scor.out',ezz(xtab[,1],ytab[,1],xynames,xtab[,2],ytab[,2],quant=seq(0,.95,.05)),env=penv);
    prScorFunc();
#     tkconfigure(svScorBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(prScorBut,state='normal',foreground='black',activeforeground='black');
#     tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
    tkentryconfigure(prMenu,prscorelabel,state='normal');
    tkentryconfigure(prMenu,pralllabel,state='normal');
    tkentryconfigure(fileSaveMenu,savescorelabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(fileSaveMenu,savealllabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(doMenu,dosummlabel,state='normal',foreground='blue',activeforeground='blue');
#    tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');
# 		tkconfigure(plQSurvBut,state='normal',foreground='black',activeforeground='black');
    tkconfigure(a,cursor='');
}

doQregFunc<-function(){
    # xtab clean
    tkconfigure(a,cursor='watch');
    doGeneralFunc(); qreg.sum<-0;class(qreg.sum)<-"try-error"; i<-1;
    xytab<-rbind(cbind(xtab,group=1),cbind(ytab,group=2));
    iqreg<-crq(Surv(log(xytab[,1]),xytab[,2])~xytab[,3],method="Por");
    while(class(qreg.sum)=="try-error"){
	if(i>length(qincr)){
	    out<-"Unable to perform quantile regression.";
	    next;
	}
	# qincr is a series of bandwidths, declared at the beginning
	# here, it is used to construct a sequence of tau values for rq
	# that don't crash it (and keep trying until it exhausts the
	# reasonable bandwidths)
	#iqreg<-rq(log(xy)~group,tau=seq(qincr[i],1-qincr[i],qincr[i]));i<-i+1;
	#qreg.sum<-try(summary(iqreg,se=tclvalue('qrse'))); cat('.');
	qreg.sum<-try(summary(iqreg,seq(qincr[i],1-qincr[i],qincr[i])));i<-i+1;
    }
    if(class(qreg.sum)!="character"){
	#qreg.sig=F;
	qreg.out<-data.frame(t(sapply(qreg.sum,function(x){
	    return(c(quantile=x$tau,x$coefficients[2,]))
	})));
	qreg.out$adjp<-p.adjust(qreg.out[,7],'holm');
	#tkconfigure(plQregBut,state='normal');
	tkentryconfigure(plotMenu,plotqrlabel,state='normal');
# 	tkconfigure(svQregBut,state='normal',foreground='blue',activeforeground='blue');
# 	tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
# 	tkconfigure(prQregBut,state='normal',foreground='black',activeforeground='black');
# 	tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
	tkentryconfigure(prMenu,prqrlabel,state='normal');
	tkentryconfigure(prMenu,pralllabel,state='normal');	tkentryconfigure(fileSaveMenu,saveqrlabel,state='normal',foreground='blue',activeforeground='blue');
	tkentryconfigure(fileSaveMenu,savealllabel,state='normal',foreground='blue',activeforeground='blue');
	
	tkentryconfigure(doMenu,dosummlabel,state='normal',foreground='blue',activeforeground='blue');
	assign('qreg.sum',qreg.sum,env=penv);
    }else{qreg.out<-qreg.sum;}
    assign('qreg.out',qreg.out,env=penv);
    assign('qreg.sig',na.exclude(qreg.out[qreg.out[,8]<.05,]),env=penv);
    #tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');
    prQregFunc();
    tkconfigure(a,cursor='');
}
doMortFunc<-function(){
    tkconfigure(a,cursor='watch');
    doGeneralFunc();
    assign('mort.out',fp(xtab[,1],ytab[,1],cx=xtab[,2],cy=ytab[,2],groupnames=xynames,perms=mperms), env=penv);
    assign('mort.text.out', gsub('-{2,}', '\n', gsub(' {2,}', '\t', capture.output(mtemp<-print(mort.out)))), env=penv);
    assign('mort.summary.out',mtemp,env=penv);
    prMortFunc();
#     tkconfigure(svMortBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(svAllBut,state='normal',foreground='blue',activeforeground='blue');
#     tkconfigure(prMortBut,state='normal',foreground='black',activeforeground='black');
#     tkconfigure(prAllBut,state='normal',foreground='black',activeforeground='black');
#     tkconfigure(doSummBut,text='Update Report',foreground='blue',activeforeground='blue');

    tkentryconfigure(prMenu,prmortlabel,state='normal');
    tkentryconfigure(prMenu,pralllabel,state='normal');
    tkentryconfigure(fileSaveMenu,savemortlabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(fileSaveMenu,savealllabel,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(doMenu,dosummlabel,state='normal',foreground='blue',activeforeground='blue');
    for(i in plotmortlabels) tkentryconfigure(plotMenu,i,state='normal');
    tkconfigure(a,cursor='');
}

doAllFunc<-function(){ doSummFunc(); doScorFunc(); doQregFunc(); doMortFunc();}

doFixXFunc<-function(){fixFunc(xtab,'xtab');}
doFixYFunc<-function(){fixFunc(ytab,'ytab');}

fixFunc<-function(d,name){
    dtemp<-fix(d);
    #browser();
    if(!all(dtemp==get(name,envir=penv))){
	assign(name,dtemp,envir=penv);
	cleanUpObjectsFunc();
    }
}

plQregFunc<-function(){
    tkconfigure(a,cursor='watch');
    doGeneralFunc();
    plot(qreg.sum,nrow=1,ncol=1,parm=2,ols=F,ylim=qregylim,
	    main=paste(tmain,"Quantile Regression"),sigf=qsigf);
    #qreg.sig<-na.exclude(qreg.out[qreg.out[,5]<.05,]);
    #points(qreg.sig[,1],qreg.sig[,2],col='red',lwd=3);
#	tkconfigure(svPlotBut,state='normal',foreground='blue',activeforeground='blue');
    tkentryconfigure(plotMenu, plotqrlabel, foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
    tkconfigure(a,cursor='');
}

plSurvFunc<-function(){
#    tkconfigure(a,cursor='watch');
    doGeneralFunc();
    if(!surv_percentalive){
	plotsurv(xtab[,1],ytab[,1],xtab[,2],ytab[,2],legend=c(xynames), fun='event',lloc='bottomright',col=surv_col,bg=surv_col,pch=surv_pch, mark=surv_mark)
    } else plotsurv(xtab[,1],ytab[,1],xtab[,2],ytab[,2],legend=c(xynames),col=surv_col,bg=surv_col, pch=surv_pch,mark=surv_mark);
    #tkconfigure(svplSurvBut,state='normal',foreground='blue',activeforeground='blue');
    #tkconfigure(svPlotBut,state='normal',foreground='blue',activeforeground='blue');
    #tkconfigure(svplQregBut,state='disabled');
    tkentryconfigure(plotMenu, survlabel, foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
    tkconfigure(a,cursor='');
}

plMortSurvFitFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='srv',real_or_fit='fit',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[1], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortSurvBothFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='srv',real_or_fit='both',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[2], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortHazRealFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='haz',real_or_fit='real',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[3], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortHazFitFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='haz',real_or_fit='fit',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[4], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortHazBothFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='haz',real_or_fit='both',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[5], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortDensRealFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='den',real_or_fit='real',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[6], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortDensFitFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='den',real_or_fit='fit',col=surv_col,pch=surv_pch);
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[7], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

plMortDensBothFunc<-function(){
    doGeneralFunc();
    myargs<-list(d=mort.out,what='den',real_or_fit='both',col=surv_col,pch=surv_pch);
    if(surv_percentalive) myargs$fun<-'event';
    do.call(plot.fp,myargs);
    tkentryconfigure(plotMenu, plotmortlabels[8], foreground='black', activeforeground='black');
    tkentryconfigure(fileSaveMenu, saveimglabel, state='normal', foreground='blue', activeforeground='blue');
}

# 	plQSurvFunc<-function(){
# 		doGeneralFunc();
# 		srv<-plotsurv(xtab[,1],ytab[,1],xtab[,2],ytab[,2],pch=NA,lwd=3);
# 		abline(v=scor.out[scor.out$sig=='*','age'],col='green');
# 		tkconfigure(svPlotBut,state='normal',foreground='blue',activeforeground='blue');
# 		legend(x='bottomleft',col=c('black','darkred','green'),lwd=c(3,3,1),legend=c(xynames,'significant'),bg='white');
# 	}

# buttons for some of the above functions
doFixXBut<-tkbutton(frame1,text=paste("Edit",xynames[1]),command=doFixXFunc,font=ar10);
doFixYBut<-tkbutton(frame1,text=paste("Edit",xynames[2]),command=doFixYFunc,font=ar10);
# doSummBut<-tkbutton(frame2,text="Run",command=doSummFunc,font=ar10);
# doScorBut<-tkbutton(frame2,text="Run",command=doScorFunc,font=ar10);
# doQregBut<-tkbutton(frame2,text="Run",command=doQregFunc,font=ar10);
# doMortBut<-tkbutton(frame2,text="Run",command=doMortFunc,font=ar10);
# doAllBut<-tkbutton(frame2,text="Run All",command=doAllFunc,font=ar10);
# plQregBut<-tkbutton(frame2,text="Draw",
# 		    command=plQregFunc,font=ar10,state='disabled');
# plSurvBut<-tkbutton(frame2,text="Draw",command=plSurvFunc,font=ar10);
# 	plQSurvBut<-tkbutton(frame2,text="Draw",command=plQSurvFunc,font=ar10,state='disabled');

delFunc<-function(){
	last<-as.numeric(tclvalue(tcl(textout,'index','end')));
	tcl(textout,'delete','1.0',paste(last,'end',sep='.'));
	assign('cursor',0,envir=penv);
}

# printing functions
prGeneralFunc<-function(d,title='',na='NA',dg=4){
	# doesn't work with mort.out because it's no longer a simple data.frame; use print.fp to extract data`1
	txt<-paste(capture.output(print(d,digits=dg)),collapse='\n');
	if(na!='NA'){txt<-gsub('\\bNA\\b',na,txt);}
	assign('cursor',as.numeric(tclvalue(tcl(textout,'index','end')))-2,envir=penv);
	tkinsert(textout,'end',paste(title,'\n\n'));
	tkinsert(textout,'end',txt);
	tkinsert(textout,'end','\n\n\n');
	tcl(textout,'yview',cursor);
}

makeRefs<-function(){
	myrefs<-c('\nLiterature References',references$survomatic,references$r,references$logrank);
	if(exists('scor.out')){myrefs<-c(myrefs,references$score);}
	if(exists('qreg.out')){myrefs<-c(myrefs,references$quantreg,references$multcomp);}
	if(exists('mort.out')){myrefs<-c(myrefs,references$hazard);}
	return(myrefs);
}

prSummFunc<-function(){
	prGeneralFunc(summ.out,'Summary',na='  ');
	tempcursor<-cursor;
	tkinsert(textout,'end',paste('Log-rank','\n\n'));
	tkinsert(textout,'end',paste('Test statistic:',lrnk.out$stat,',','p =',lrnk.out$pval,
					'based on',lrnk.out$nsim,'permutations.'));
	tkinsert(textout,'end','\n\n\n');
	if(exists('scor.out')){
		prGeneralFunc(scor.out[scor.out$p<0.05,],'Score Test: Significant Quantiles');
	}
	if(exists('qreg.out')){
		prGeneralFunc(qreg.out[!is.na(qreg.out[,5])&qreg.out[,5]<0.05,],
				'Quantile Regression: Significant Quantiles');
	}
	if(exists('mort.out')){
		prGeneralFunc(mort.out,'Mortality Models');
	}
	tkinsert(textout,'end',paste(makeRefs(),collapse='\n'));
	tcl(textout,'yview',tempcursor);
}
prScorFunc<-function(){prGeneralFunc(scor.out,'Score Test');}
prQregFunc<-function(){prGeneralFunc(qreg.out,'Quantile Regression');}

prMortFunc<-function(){prGeneralFunc(mort.out,'Mortality Models');}

prAllFunc<-function(){
	tempcursor<-as.numeric(tclvalue(tcl(textout,'index','end')))-2;
	if(exists('summ.out')){prSummFunc();}
	if(exists('scor.out')){prScorFunc();}
	if(exists('qreg.out')){prQregFunc();}
	if(exists('mort.out')){prMortFunc();}
	tcl(textout,'yview',tempcursor);
}

# buttons for the print functions
prSummBut<-tkbutton(frame2,text="Show Output",command=prSummFunc,font=ar10,
		    state='disabled');
prScorBut<-tkbutton(frame2,text="Show Output",command=prScorFunc,font=ar10,
		    state='disabled');
prQregBut<-tkbutton(frame2,text="Show Output",command=prQregFunc,font=ar10,
		    state='disabled');
prMortBut<-tkbutton(frame2,text="Show Output",command=prMortFunc,font=ar10,
		    state='disabled');
prAllBut<-tkbutton(frame2,text="Show All Output",command=prAllFunc,font=ar10,
		    state='disabled');

# save functions
svGeneralFunc<-function(d,type,ext,filename=NULL,rn=T,cn=NA,q=F,append=F){
	if(is.null(filename)){
		filetypes<-paste('{{',type,'} {',ext,'}}',sep='');
		filename<-tkgetSaveFile(initialdir=tclvalue('path'),filetypes=filetypes);
		filename<-tclvalue(filename);
	}
	if(filename!=""){
		#prGeneralFunc<-function(d,title='',na='NA',dg=4,to=NULL){
		switch(ext,
			'.txt'={write.table(d, file=filename, sep='\t', quote=q, row.names=rn, col.names=cn, append=append);},
			'.rdata'={save(list=d,file=filename);},
			'.Rdata'={save(list=d,file=filename);},
			'.pdf'={dev.copy2pdf(file=filename);},
			'.eps'={dev.copy2eps(file=filename);},
			'.dat'={write(d,file=filename,append=append)}
		)
		tkfocus(a);return(filename);
	} else {tkfocus(a);return(-1)};
}

svSummFunc<-function(){
    if((filename<-svGeneralFunc(NULL,'Tab Delimited Test','.txt'))>0){
	myrefs<-c();
	svGeneralFunc('Summary\n','','.dat',filename=filename,append=T)
	svGeneralFunc(summ.out,'','.txt',filename=filename,append=T);
	svGeneralFunc(paste('\n\nLog-rank\n\nTest statistic:', lrnk.out$stat, ',' , 'p =', lrnk.out$pval, 'based on', lrnk.out$nsim, 'permutations.'), '', '.dat', filename=filename, append=T);
	if(exists('scor.out')){
	    svGeneralFunc('\n\nScore Test: Significant Quantiles\n', '', '.dat', filename=filename, append=T);
	    svGeneralFunc(scor.out[scor.out$p<0.05,],'','.txt',filename=filename, append=T);
	};
	if(exists('qreg.out')){
	    svGeneralFunc('\n\nQuantile Regression: Significant Quantiles\n', '', '.dat', filename=filename, append=T);
	    svGeneralFunc(qreg.out[!is.na(qreg.out[,5])&qreg.out[,5]<0.05,], '', '.txt', filename=filename, append=T);
	}
	if(exists('mort.text.out')){
	    svGeneralFunc('\n\nMortality Models\n', '','.dat',filename=filename,append=T);
	    svGeneralFunc(mort.text.out,'','.txt',filename=filename,rn=F,cn=F,append=T);
	};
	svGeneralFunc(makeRefs(),'','.txt',filename=filename,rn=F,cn=F,q=F,append=T);
	#tkconfigure(svSummBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,savesummlabel,foreground='black',activeforeground='black');
    }
}

svScorFunc<-function(){
    if(svGeneralFunc(scor.out,'Tab Delimited Text','.txt')>0){
	# tkconfigure(svScorBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,savescorelabel,foreground='black',activeforeground='black');
    }
}

svQregFunc<-function(){
    if(svGeneralFunc(qreg.out,'Tab Delimited Text','.txt')>0){
	# tkconfigure(svQregBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,saveqrlabel,foreground='black',activeforeground='black');
    }
}

svMortFunc<-function(){
    if(svGeneralFunc(mort.text.out,'Tab Delimited Text','.txt',rn=F,cn=F)>0){
	# tkconfigure(svMortBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,savemortlabel,foreground='black',activeforeground='black');
    }
}

svAllFunc<-function(){
    doGeneralFunc();
    tosave<-c('xtab','ytab','xynames','surv_percentalive','surv_mark',
		'surv_pch','qregylim','surv_col','logrank_perms');
    if(exists('summ.out')){tosave<-c(tosave,'summ.out');}
    if(exists('lrnk.out')){tosave<-c(tosave,'lrnk.out');}
    if(exists('scor.out')){tosave<-c(tosave,'scor.out');}
    if(exists('qreg.out')){tosave<-c(tosave,'qreg.out');}
    if(exists('qreg.sum')){tosave<-c(tosave,'qreg.sum');}
    if(exists('mort.out')){tosave<-c(tosave,'mort.out');}
    if(exists('mort.summary.out')){tosave<-c(tosave,'mort.summary.out');}
    if(exists('mort.text.out')){tosave<-c(tosave,'mort.text.out');}
    if(svGeneralFunc(tosave,'R Data','.rdata')>0){
	# tkconfigure(svAllBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,savealllabel,foreground='black',activeforeground='black');
    }
}

# 	svplQregFunc<-function(){
# 		if(svGeneralFunc(NULL,'Encapsulated Postscript','.eps')>0){
# 			tkconfigure(svplQregBut,foreground='black',activeforeground='black');
# 		}
# 	}
# 
# 	svplSurvFunc<-function(){
# 		if(svGeneralFunc(NULL,'Encapsulated Postscript','.eps')>0){
# 			tkconfigure(svplSurvBut,foreground='black',activeforeground='black');
# 		}
# 	}

svPlotFunc<-function(){
    if(svGeneralFunc(NULL,'PDF','.pdf')>0){
	# tkconfigure(svPlotBut,foreground='black',activeforeground='black');
	tkentryconfigure(fileSaveMenu,saveimglabel,foreground='black',activeforeground='black');
    }
}
configFunc<-function(){
	b<-tktoplevel();
	tkwm.title(b,'Survomatic Settings');
	surv_percentalive.cb <- tkcheckbutton(b,variable="surv_percentalive");
	surv_mark.txt <- tkentry(b,width=3,textvariable="surv_mark");
	surv_pch1.txt <- tkentry(b,width=3,textvariable="surv_pch1");
	surv_pch2.txt <- tkentry(b,width=3,textvariable="surv_pch2");
	qregylim1.txt <- tkentry(b,width=3,textvariable="qregylim1");
	qregylim2.txt <- tkentry(b,width=3,textvariable="qregylim2");
	surv_col1.txt <- tkentry(b,width=10,textvariable="surv_col1");
	surv_col2.txt <- tkentry(b,width=10,textvariable="surv_col2");
	logrank_perms.txt <- tkentry(b,width=20,textvariable="logrank_perms");
	fontsize.txt <- tkentry(b,width=3,textvariable="pcx");
	configOKfunc<-function(){
		tkconfigure(textout,font=paste('courier',tclvalue("pcx")));
		tkdestroy(b);
	}
	configPCHfunc<-function(){
		plot(1:6,1:6,type='n',axes=F,xlab='',ylab='',
			main='Numeric Codes For Valid Plot Symbols');
		k<- -1;
		for(i in 6:1){
			for(j in 1:5){
				k<-k+1;points(j,i+0.1,pch=ifelse(k<26,k,''),cex=2);
				text(j,i-0.1,ifelse(k<26,k,''));
			}
		}
	}
	configOKbut<-tkbutton(b,text='OK',command=configOKfunc);
	configPCHbut<-tkbutton(b,text='?',command=configPCHfunc);
	tkgrid(tklabel(b,text="Plot survival from 100% to 0%"),surv_percentalive.cb);
	tkgrid(tklabel(b,text="Symbol to use for censored data"),surv_mark.txt);
	tkgrid(tklabel(b,text=paste('Symbols for',
				    xynames[1],'and',xynames[2])),surv_pch1.txt,surv_pch2.txt,
				    configPCHbut);
	tkgrid(tklabel(b,text=paste('Line colors for',
				    xynames[1],'and',xynames[2])),surv_col1.txt,surv_col2.txt);
	tkgrid(tklabel(b,text="Lower and upper limits for the Y-axis of the quantile regression plot"),qregylim1.txt,qregylim2.txt);
	tkgrid(tklabel(b,text="Number of permutations to use for log-rank"),
			logrank_perms.txt,columnspan=2);
	tkgrid(tklabel(b,text="Font size for output"),fontsize.txt,columnspan=2);
	tkgrid(configOKbut);
	doGeneralFunc()
}

# buttons for the save functions
# svSummBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svSummFunc,font=ar10,state='disabled');
# svScorBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svScorFunc,font=ar10,state='disabled');
# svQregBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svQregFunc,font=ar10,state='disabled');
# svMortBut<-tkbutton(frame2,text="Save Output (tab delimited spreadsheet)",command=svMortFunc,font=ar10,state='disabled');
# svAllBut<-tkbutton(frame2,text="Save All Output (.Rdata file)",command=svAllFunc,font=ar10,state='disabled');
# 	svplQregBut<-tkbutton(frame2,text="Save Image (.eps file)",command=svplQregFunc,font=ar10,state='disabled');
# 	svplSurvBut<-tkbutton(frame2,text="Save Image (.eps file)",command=svplSurvFunc,font=ar10,state='disabled');
#svPlotBut<-tkbutton(frame2,text="Save Current Image (PDF)",command=svPlotFunc,font=ar10,state='disabled');
configBut<-tkbutton(frame2,text="Options",command=configFunc,font=ar10);

# begin layout
topMenu<-tkmenu(a);
tkconfigure(a,menu=topMenu);
fileMenu<-tkmenu(topMenu);
fileSaveMenu<-tkmenu(fileMenu);
doMenu<-tkmenu(topMenu);
prMenu<-tkmenu(topMenu);
plotMenu<-tkmenu(topMenu);
tkadd(fileSaveMenu,"command",label=savesummlabel,command=svSummFunc,state='disabled');
tkadd(fileSaveMenu,"command",label=savescorelabel,command=svScorFunc,state='disabled');
tkadd(fileSaveMenu,"command",label=saveqrlabel,command=svQregFunc,state='disabled');
tkadd(fileSaveMenu,"command",label=savemortlabel,command=svMortFunc,state='disabled');
tkadd(fileSaveMenu,"command",label=savealllabel,command=svAllFunc,state='disabled');
tkadd(fileSaveMenu,"command",label=saveimglabel,command=svPlotFunc,state='disabled');
tkadd(fileMenu,"command",label="Load WinModest Data",command=function() loadWMFunc());
tkadd(fileMenu,"command",label="Load R Data",command=function() loadFunc());
tkadd(fileMenu,"cascade",label="Save...",menu=fileSaveMenu);
tkadd(fileMenu,"command",label="Quit",command=quitFunc);
tkadd(doMenu,"command",label=dosummlabel,command=doSummFunc);
tkadd(doMenu,"command",label=doscorelabel,command=doScorFunc);
tkadd(doMenu,"command",label=doqrlabel,command=doQregFunc);
tkadd(doMenu,"command",label=domortlabel,command=doMortFunc);
tkadd(doMenu,"command",label=doalllabel,command=doAllFunc);
tkadd(prMenu,"command",label=prsummlabel,command=prSummFunc,state='disabled');
tkadd(prMenu,"command",label=prscorelabel,command=prScorFunc,state='disabled');
tkadd(prMenu,"command",label=prqrlabel,command=prQregFunc,state='disabled');
tkadd(prMenu,"command",label=prmortlabel,command=prMortFunc,state='disabled');
tkadd(prMenu,"command",label=pralllabel,command=prAllFunc,state='disabled');
tkadd(plotMenu,"command",label=survlabel,command=plSurvFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[1],command=plMortSurvFitFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[2],command=plMortSurvBothFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[3],command=plMortHazRealFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[4],command=plMortHazFitFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[5],command=plMortHazBothFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[6],command=plMortDensRealFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[7],command=plMortDensFitFunc,state='disabled');
tkadd(plotMenu,"command",label=plotmortlabels[8],command=plMortDensBothFunc,state='disabled');
tkadd(plotMenu,"command",label=plotqrlabel,command=plQregFunc,state='disabled');
tkadd(topMenu,"cascade",label="File",menu=fileMenu);
tkadd(topMenu,"cascade",label="Analysis",menu=doMenu);
tkadd(topMenu,"cascade",label="Reports",menu=prMenu);
tkadd(topMenu,"cascade",label="Plot",menu=plotMenu);
# tkgrid(tklabel(frame2,text='Summary Report',font=ar10,anchor='e'),
# 	doSummBut,prSummBut,sticky='we',padx=1,pady=3);
# # tkgrid(tklabel(frame2,text='Survival Curves',font=ar10,anchor='e'),
# # 	tklabel(frame2,text=''),plSurvBut,tklabel(frame2,text=''),sticky='we',padx=1,pady=3);
# tkgrid(tklabel(frame2,text='Score Test',font=ar10,anchor='e'),
# 	doScorBut,prScorBut,sticky='we',padx=1,pady=3);
# # 	tkgrid(tklabel(frame2,text='Surv. Curves w/ Signif. Quantiles',font=ar10,anchor='e'),tklabel(frame2,text=''),plQSurvBut,tklabel(frame2,text=''),sticky='we',padx=1,pady=3);
# tkgrid(tklabel(frame2,text='Quantile Regression',font=ar10,anchor='e'),
# 	doQregBut,prQregBut,sticky='we',padx=1,pady=3);
# # tkgrid(tklabel(frame2,text='Quantile Regression Plot',font=ar10,anchor='e'),
# # 	tklabel(frame2,text=''),plQregBut,tklabel(frame2,text=''),sticky='we',padx=1,pady=3);
# tkgrid(tklabel(frame2,text='Mortality Models',font=ar10,anchor='e'),doMortBut,prMortBut,sticky='we',padx=1,pady=3);
# tkgrid(tklabel(frame2,text='All Tests',font=ar10,anchor='e'),doAllBut,prAllBut,sticky='we',padx=1,pady=3);

# frame1 variables
blank1<-tklabel(frame1,text='');blank2<-tklabel(frame1,text='');
#	blank1left<-tklabel(frame1left,text=' ');
# 	datatxt1<-tklabel(frame1left,text='Age at Death',font='Arial 8',padx=5,anchor='center');
# 	censtxt1<-tklabel(frame1left,text='Censor',font='Arial 8',padx=18,anchor='center');
# 	datatxt2<-tklabel(frame1right,text='Age at Death',font='Arial 8',padx=5,anchor='center');
# 	censtxt2<-tklabel(frame1right,text='Censor',font='Arial 8',padx=18,anchor='center');
tcl("set","name1","Control"); tcl("set","name2","Experimental");
name1txt<-tkentry(frame1,textvariable="name1",font=ar8,bg='white'); 
name2txt<-tkentry(frame1,textvariable="name2",font=ar8,bg='white');
# vard1<-tclArray();vard2<-tclArray();
# vard1[[0,0]]<-'Age at Death'; vard1[[0,1]]<-'Censor';
# vard2[[0,0]]<-'Age at Death'; vard2[[0,1]]<-'Censor';
# for(i in 1:n){vard1[[i,0]]<-'-';vard1[[i,1]]<-1;vard2[[i,0]]<-'-';vard2[[i,1]]<-1;}

# if data passed as r objects from environment
if(!is.null(xtab)&!is.null(ytab)){
#     assign('n',max(dim(xtab)[1],dim(ytab)[1],n),env=penv);
#     array2tcl(xtab,vard1,1:dim(xtab)[1],0:1);
#     array2tcl(ytab,vard2,1:dim(ytab)[1],0:1);
    # tkconfigure(plSurvBut,state='active');
    tkentryconfigure(plotMenu,survlabel,state='normal');
    tcl("set","name1",xynames[1]); tcl("set","name2",xynames[2]);
}

# frame2 variables
tcl("set","save",1);tcl("set","path",getwd());tcl("set","tscor",1); 
tcl("set","tqreg",1);tcl("set","tmort",1);tcl("set","tsumm",1); 
pathl<-tklabel(frame4,text='File Path:',anchor='e',font=ar10);
# set up the input tables for the survival times in frame1
# d1<-tkwidget(frame1,'table',height=13,rows=n,cols=2,titlerows=0,titlecols=0,resizeborders='none',
# 		colwidth=10,titlerows=1,yscrollcommand=function(...) tkset(yscr1,...));
# yscr1<-tkscrollbar(frame1,command=function(...)tkyview(d1,...));
# d2<-tkwidget(frame1,'table',height=13,rows=n,cols=2,titlerows=0,titlecols=0,resizeborders='none',
# 		colwidth=10,titlerows=1,yscrollcommand=function(...) tkset(yscr2,...));
# yscr2<-tkscrollbar(frame1,command=function(...)tkyview(d2,...));
# set up widgets in frame2
pathtxt<-tkentry(frame4,textvariable="path",font=ar8,bg='white');
# 	savechk<-tkcheckbutton(frame2,variable="save", text='Save',font=ar10,anchor='w');
# 	scorchk<-tkcheckbutton(frame2,variable="tscor",text='Score Test',font=ar10,anchor='w');
# 	qregchk<-tkcheckbutton(frame2,variable="tqreg",text='Quant Reg',font=ar10,anchor='w');
# 	mortchk<-tkcheckbutton(frame2,variable="tmort",text='Mortality',font=ar10,anchor='w');
# 	summchk<-tkcheckbutton(frame2,variable="tsumm",text='Summary',font=ar10,anchor='w');
# 	runBut<-tkbutton(frame2,text="Run!",command=runFunc);
#dbBut<-tkbutton(frame4,text="Load WinModest Data",command=function() loadWMFunc(),padx=3);
#loadBut<-tkbutton(frame4,text="Load R Data",command=function() loadFunc(),padx=3);
delBut<-tkbutton(frame4,text='Clear Output Window',command=delFunc,padx=3);
#quitBut<-tkbutton(frame4,text="Quit",command=quitFunc,padx=3);
# set up frame 3
xscr<-tkscrollbar(frame3,command=function(...)tkxview(textout,...),orient='horizontal');
yscr3<-tkscrollbar(frame3,command=function(...)tkyview(textout,...));
textout<-tktext(frame3,bg='white',font=paste('courier',tclvalue("pcx")),wrap='none',width=165,
		xscrollcommand=function(...) tkset(xscr,...),
		yscrollcommand=function(...) tkset(yscr3,...));
# arrange the input tables in frame1
tkgrid(tklabel(frame2,text=''),tklabel(frame2,text=''),configBut,sticky='we',padx=1,pady=3);
tkgrid(doFixXBut,doFixYBut);
tkgrid(name1txt,name2txt);
tkgrid(pathl,pathtxt,delBut);
#tkgrid(delBut,pady=2,padx=5,columnspan=1);
tkwm.title(a,paste('Survomatic v', installed.packages()[installed.packages()[,'Package']=='Survomatic','Version'],sep=''));
if(!is.null(dat)){loadFunc(dat);}
# arrange widgets in frame3
tkgrid(textout,yscr3); tkgrid(xscr); tkgrid.configure(textout,sticky='we');
tkgrid.configure(xscr,sticky='we'); tkgrid.configure(yscr3,sticky='ns');
# arrange frames
tkgrid(frame4,columnspan=4);
tkgrid(frame1,frame2,padx=2,pady=2);
tkgrid(frame3,columnspan=2,sticky='nwe');
tkfocus(a);
browser();
}