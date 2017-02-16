#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// double bebxf(double *b, int *x){return(exp(*b * *x));}
// double nbebxf(double *b, int *x){return(1);}
// double baebxf(double *a, double *b, double *ebx){return(*a * (*ebx - 1)/ *b);}
// double nbaebxf(double *a, double *b, int *x, double *ebx){return(*a * *x);}
// double ssrvf(double *a, double *c, double *s, double *aebx, int *x){
// 	return(log(exp(- *c * *x)/pow((*s * *aebx + 1),(1/ *s))));
// }
// double nssrvf(double *a, double *c, double *s, double *aebx, int *x){- *aebx - *c * *x;}
// 
// void ofx(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
// {
// 	int n = *size;
// 	double ebx, aebx;
// 	double bebxf(), nbebxf(), baebxf(), nbaebxf(), ssrvf(), nssrvf();
// 	double (*ebxf)(), (*aebxf)(), (*srvf)();
// 	printf("Declarations done\n");
// 	if(*b > 2E-19){ebxf = bebxf; aebxf = baebxf;
// 	printf("B is not zero\n");
// 	} else {printf("B is zero\n"); ebxf = nbebxf; aebxf = nbaebxf;}
// 	if(*s > 6E-10){printf("S is not zero\n");srvf = ssrvf;
// 	} else {printf("S is zero\n");srvf = nssrvf;}
// 	while(--n >= 0){
// 		ebx = (*ebxf)(b,x[n]); aebx = (*aebxf)(a,b,x,&ebx);
// 		printf("exb is %f5.5 and aebx is %f5.5\n",ebx, aebx);
// 		*ans = *ans + (*srvf)(a,c,s,&aebx,x[n]);
// 		printf("so far ans is %f5.5\n",*ans);
// 		*ans = *ans + log(*c + *a * ebx/(*s * aebx + 1)) * censor[n];
// 		printf("now ans is %f5.5\n",*ans);
// 	}
// }

void ofw(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	while(--n >= 0){
		*ans = *ans + log(pow(*a,*b) * (*b) * pow(x[n],(*b -1))) * censor[n];
		*ans = *ans - pow(*a * x[n],*b);
	}
}

void oflm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double ebx, aebx;
	while(--n >= 0){
		if(*b < 2E-19){ebx = 1; aebx = *a * x[n];
		} else {ebx = exp(*b * x[n]); aebx = *a * (ebx -1)/ *b; }
		*ans = *ans + log(*c + *a * ebx/(*s * aebx + 1)) * censor[n];
		if(*s < 6E-10){ *ans = *ans - aebx - *c * x[n];
		} else {*ans = *ans + log(exp(- *c * x[n])/pow((*s * aebx + 1),(1/ *s)));}
	}
}

void ofgm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	double mys = 0;
	oflm(a, b, c, &mys, x, size, censor, ans);
}

void ofg(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	double myc = 0; double mys = 0;
	oflm(a, b, &myc, &mys, x, size, censor, ans);
}

void ofl(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	double myc = 0;
	oflm(a, b, &myc, s, x, size, censor, ans);
}

void gw(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double xpowb1;
	double ab = pow(*a, *b); double abb = *b * ab;
	while ( --n >= 0){
		suma = suma + censor[n] * *b / *a;
		suma = suma - *b * pow(*a * x[n],*b)/ *a;
		xpowb1 = pow(x[n],*b - 1);
		sumb = sumb + censor[n] * pow(x[n],1 - *b) * 
			(xpowb1 * log(x[n]) * abb + log(*a) * xpowb1 * abb + xpowb1 * ab)/abb;
		sumb = sumb - pow(*a * x[n],*b) * log(*a * x[n]);
	}
	ans[0] = suma; ans[1] = sumb;
}


void gg(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb;
	double expBX; double expBXminus1;
	double pow2B = pow(*b,2);
	while ( --n >= 0){
		expBX = exp(*b * x[n]); expBXminus1 = expBX -1;
		suma = suma + censor[n]/ *a;
		suma = suma - expBXminus1/ *b;
		sumb = sumb + x[n] * censor[n];
		sumb = sumb - *a * (expBX * x[n]/ *b - expBXminus1/pow2B);
	}
	ans[0] = suma; ans[1] = sumb;
}
// some values disagree; probably left out a term!
void ggm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; 
	double expBX, expBXminus1divB, AexpBX;
	while ( --n >= 0){
		expBX = exp(*b * x[n]); expBXminus1divB = (expBX - 1)/ *b; 
		AexpBX = *a * expBX;
		suma = suma + censor[n]/(*a + *c/expBX); //ok
		suma = suma - expBXminus1divB; //ok
		sumb = sumb + censor[n] * x[n]/(1+ *c * expBX/ *a); //ok
		sumb = sumb + *a * expBXminus1divB/ *b - x[n] * AexpBX/ *b ;//ok
		sumc = sumc + censor[n]/(*c + AexpBX); //ok
		sumc = sumc - x[n]; //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc;
}
void gl(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma, sumb, sumc, sums;
	double expBX, expBXminus1, prdASexpBXminus1addB;
	double pow2B = pow(*b,2); double pow2S = pow(*s,2); 
	double prdAB = *a * *b; double prdAS = *a * *s; double prdABS = prdAB * *s;
	while ( --n >= 0){
		expBX = exp(*b * x[n]);
		expBXminus1 = expBX - 1;
		prdASexpBXminus1addB = prdAS * expBXminus1 + *b;
		suma = suma + censor[n] * *b/(*a * prdASexpBXminus1addB); //ok
		suma = suma - expBXminus1/(prdASexpBXminus1addB); //ok
		sumb = sumb + censor[n] * (prdAS * expBX + (pow2B - prdABS) *x[n] - prdAS)/(*b * prdASexpBXminus1addB); //ok
		sumb = sumb - ((prdAB * x[n] - *a) * expBX + *a)/(*b * prdASexpBXminus1addB); //ok
		sums = sums - censor[n] * *a * expBXminus1/prdASexpBXminus1addB; //ok
		sums = sums + (prdASexpBXminus1addB * log(prdASexpBXminus1addB/ *b) - prdAS * expBXminus1)/(pow2S * prdASexpBXminus1addB); //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sums;
}

void glm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; double sums;
	double exp1B = exp(1) * *b; 
	double pow2S = pow(*s,2); double pow2B = pow(*b,2); 
	double prdAS = *a * *s; double prdABS = prdAS * *b;
	double prdACS = prdAS * *c; double prdAB = *a * *b;
	double prdBC = *b * *c;
	double expBX, prdASexpBX, prdACSexpBX, prdABexpBX, BminusprdAS;
	while ( --n >= 0){
		expBX = exp(*b * x[n]); prdASexpBX = prdAS * expBX;
		prdACSexpBX = prdACS * expBX; prdABexpBX = prdAB * expBX;
		BminusprdAS = *b - prdAS;

		suma = suma + censor[n] * (pow2B * expBX)/((prdASexpBX + BminusprdAS)*(prdACSexpBX + prdABexpBX - prdACS + prdBC)); //ok
		suma = suma - expBX/(prdASexpBX + exp1B); //ok
		sumb = sumb + censor[n] * (*a * expBX * (prdASexpBX - prdABS * x[n] + pow2B * x[n] - prdAS))/((prdASexpBX + BminusprdAS)*(prdACSexpBX +prdABexpBX - prdACS + prdBC)); //ok
		sumb = sumb - ((prdAB * x[n] - *a) * expBX)/(prdABS * expBX + exp1B * *b); //ok
		sumc = sumc + censor[n] * (prdASexpBX + BminusprdAS)/(prdACSexpBX + prdABexpBX - prdACS + prdBC); //ok
		sumc = sumc - x[n];
		sums = sums - censor[n] * (*a * prdAB * (expBX-1) * expBX)/((prdASexpBX + BminusprdAS)*(prdACSexpBX + prdABexpBX - prdACS + prdBC)); //ok
		sums = sums + ((prdASexpBX + exp1B)*log(((prdASexpBX + exp1B))/exp1B)- prdASexpBX)/(pow2S * prdASexpBX + exp1B * pow2S); //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc; ans[3] = sums;
}

/* Need to rename these functions to gmsrv, lmsrv, wsrv */
double gmsrv(double *a, double *b, double *c, double *s, int x)
{return(exp(- (*a * (exp(*b * x) - 1)/ *b) - *c * x));}
double lmsrv(double *a, double *b, double *c, double *s, int x)
{return(exp(- *c * x)*pow((1+ *s * *a * (exp(*b * x)-1)/ *b),(-1/ *s)));}
double lmsrv2(double *a, double *b, double *c, double *s, int x)
{double asb = *a * *s/ *b; double expbx = exp(*b * x); double expcx = exp(- *c * x); 
	double denom = pow(asb * expbx - asb + 1, 1/ *s);
return(exp(- *c * x)/pow((asb * exp(*b * x) - asb + 1),1/ *s));}
double wsrv(double *a, double *b, double *c, double *s, int x)
{double ax = *a * x; double axb = pow(ax,*b); double expaxb = exp(-axb);return(expaxb);}

void srvshp(double *a, double *b, double *c, double *s, int *size, int *model, int *x, double *ans){
	int n = *size; double surv;
	double (*sf)(double *,double *,double *,double *,int);
	sf = &wsrv;
	switch(*model){
		case 1: sf = &wsrv; break;
		case 2: sf = &lmsrv; break;
		case 3: sf = &gmsrv; break;
		case 4: sf = &lmsrv2; break;
		default: printf("Warning: no case selected!\n");	
	}
	while(--n >= 0){
		ans[n] = (*sf)(a,b,c,s,x[n]);
	}
}


void simsurv(double *a, double *b, double *c, double *s, int *size, int *model, double *ans)
{
	int n = *size; int i;
	double surv, rn;
	double (*sf)(double *,double *,double *,double *,int);
	sf = &wsrv;
	switch(*model){
		case 1: sf = &wsrv; break;
		case 2: sf = &lmsrv; break;
		case 3: sf = &gmsrv; break;
		case 4: sf = &lmsrv2; break;
		default: printf("Warning: no case selected!\n");	
	}
	while(--n >= 0){
		i = 1; surv = 1; rn = 0;
		while(surv > rn){
			surv = (*sf)(a,b,c,s,i+1)/(*sf)(a,b,c,s,i);
			rn = (double)rand()/((double)(RAND_MAX)+(double)(1));
			ans[n] = i; i++;
		}
	}
}

void zptest(double *x1, double *x2, double *n1, double *n2, double *n, int *ln, double *that1, double *that2, double *ttil, double *zp)
{
	int i; double denom;
	for(i=0;i< *ln; i++){
		that1[i] = x1[i]/ *n1; 
		that2[i] = x2[i]/ *n2;
		ttil[i] = (x1[i]+ x2[i])/ *n;
		zp[i] = (that1[i] - that2[i])/sqrt(ttil[i] * (1 - ttil[i]) * (1/ *n1 + 1/ *n2));
	}
}
