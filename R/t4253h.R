t4253h<-function(x,firstpass=T,cycles=1){
	lx<-length(x);
	med4<-c(x[1],apply(embed(x,4),1,median),median(x[(lx-1):lx]),x[lx]);
	z1<-c(x[1],apply(embed(med4,2),1,median));
	z2<-c(z1[1],median(z1[1:3]),apply(embed(z1,5),1,median),median(z1[(lx-2):lx]),z1[lx]);
	z3<-apply(embed(z2,3),1,median);
	z3<-c(median(c(3*z3[1]-2*z3[2],z2[1],z3[1])),z3);
	z3<-c(z3,median(c(3*z3[lx-1]-2*z3[lx-2],z2[lx],z3[lx-1])));
	z4<-0.25*z3[1:(lx-2)]+0.5*z3[2:(lx-1)]+0.25*z3[3:lx];
	z4<-c(z3[1],z4,z3[lx]);
	res<-x-z4;
	if(firstpass){res4<-t4253h(res,FALSE);} else {return(z4);}
	out<-z4+res4;
	if(cycles==1){cat('\n');return(out);} else {cat(cycles,' ');return(t4253h(out,cycles=cycles-1));}
}