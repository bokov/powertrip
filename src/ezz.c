#include <math.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

void testprobs(double *x1,double *x2,double *n1, double *n2, int *m, double *th, double *myans){
	int i;
	printf("x1: ");
	for(i=0;i<= *n1;i++){printf(" %f",x1[i]);} printf("\n");
	printf("x2: ");
	for(i=0;i<= *n2;i++){printf(" %f",x2[i]);} printf("\n");
	printf("n1: %f\n",*n1);
	printf("n2: %f\n",*n2);
	printf("m: %d\n",*m);
	printf("th: %f\n",*th);
}

void zpprob(double *x1, double *x2, double *n1, double *n2, int *ln, double *th, double *ans){
	int i;
	for(i=0;i<= *ln; i++){
		ans[i] = exp(lchoose(*n1,x1[i])+lchoose(*n2,x2[i])+(x1[i]+x2[i])*log(*th)+(*n1+ *n2-x1[i]-x2[i])*log(1- *th));
	}
}

void runzpprobs(double *n1, double *n2, double *th1, double *th2, double *step, double *out, int *rounds, double *zp, double *zptests){
	int i; double j; int k; int l; int m = 0; int m0 = (*n1+1)*(*n2+1); int n; double jout;
	double *x1 = malloc(m0*sizeof(double));
	double *x2 = malloc(m0*sizeof(double));
	int *keep = malloc(m0*sizeof(int));
	for(l=0;l<= *n2;l++){
		for(k=0;k<= *n1;k++){
			x1[m] = k; x2[m] = l;
			m++;
		}
	}
	double *myans = malloc((m+1)*sizeof(double));
	for(i=0; i< *rounds; i++){
		out[i] = 0;
		for(n=0;n<=m;n++){if(labs(zptests[n])>labs(zp[i])){keep[n]=1;}else{keep[n]=0;}}
		for(j= *th1; j<= *th2; j+= *step){
			jout = 0;
			zpprob(x1,x2,n1,n2,&m,&j,myans);
			// below added myans[l] = 0 to try and fix incorrect values when rounds > 1
			// really need to figure out why this was happening even if this fix works
			// not yet tried as of 110625
			for(l=0;l<= m;l++){if(keep[l]>0){jout += myans[l]; myans[l] = 0;}}
			if(out[i] < jout){out[i] = jout;}
		}
		printf(".");
	}
}

