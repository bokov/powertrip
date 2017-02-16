#include <math.h>
#include <stdio.h>

void ofgm2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size; double ebx = 0;
	while ( --n >= 0){
		ebx = exp(*b * x[n]);
		*ans = *ans + log(*c + *a * ebx) * censor[n];
		*ans = *ans -(*c * x[n]) - *a * (ebx - 1) / *b;
	}
}

void ofg2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	*c = 0;
	ofgm2(a, b, c, s, x, size, censor, ans);
}

/* The objective function for the Logistic and Logistic-Makeham models (if c==0) */
void oflm2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size; double ebx = 0; double saebxb = 0;
	while ( --n >= 0){
		ebx = exp(*b * x[n]); saebxb = 1 + *s * *a * (ebx - 1)/ *b;
		*ans = *ans + log((*c + *a * ebx/saebxb)) * censor[n];
		*ans = *ans + log(exp(- *c * x[n]) * pow(saebxb,(- 1/ *s)));
	}
}

void ofl2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	*c = 0;
	oflm2(a, b, c, s, x, size, censor, ans);
}

void gg2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb;
	while ( --n >= 0){
		suma = suma + censor[n]/ *a;
		suma = suma -((exp(*b * x[n]) - 1)/ *b);
		sumb = sumb + x[n] * censor[n];
		sumb = sumb -(*a * (exp(*b * x[n]) * x[n])/ *b - *a * (exp(*b * x[n]) - 1)/pow( *b ,2));
	}
	ans[0] = suma; ans[1] = sumb;
}
// some values disagree; probably left out a term!
void ggm2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; double ebx = 0; double aebx = 0;
	while ( --n >= 0){
		ebx = exp(*b * x[n]); aebx = *a * ebx;
		suma = suma + 1/(*a + *c/ebx) * censor[n];
		suma = suma -((ebx - 1)/ *b);
		sumb = sumb + x[n]/(1+ *c/aebx) * censor[n];
//		sumb = sumb + *a * (exp(*b * x[n]) * x[n])/(*c + *a * exp(*b * x[n])) * censor[n];
		sumb = sumb -(aebx * x[n]/ *b - *a * (ebx - 1)/pow(*b,2));
//		sumb = sumb -(*a * (exp(*b * x[n]) * x[n])/ *b - *a * (exp(*b * x[n]) - 1)/pow( *b ,2));
		sumc = sumc + (1/(*c + aebx)) * censor[n];
//		sumc = sumc + (1/(*c + *a * exp(*b * x[n]))) * censor[n];
		sumc = sumc - x[n];
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc;
}
void gl2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; double sums; double ebx = 0;
	while ( --n >= 0){
		ebx = exp(*b * x[n]);
		//suma = suma + ((1 - *a * *s * (ebx - 1)/(*b * (1 + *s * *a * (ebx - 1)/ *b)))) * censor[n];
		suma = suma + (ebx/(1 + *s * *a * (ebx - 1)/ *b) - *a * ebx * (*s * (ebx - 1)/ *b)/pow((1 + *s * *a * (ebx - 1)/ *b),2))/(*a * ebx/(1 + *s * *a * (ebx - 1)/ *b)) * censor[n];
		suma = suma + pow((1 + *s * *a * (ebx - 1)/ *b),((-1/ *s) - 1)) * ((-1/ *s) * (*s * (ebx - 1)/ *b))/pow((1 + *s * *a * (ebx - 1)/ *b),(-1/ *s));
		sumb = sumb + (*a * (ebx * x[n])/(1 + *s * *a * (ebx - 1)/ *b) - *a * ebx * (*s * *a * (ebx * x[n])/ *b - *s * *a * (ebx - 1)/pow(*b,2))/pow((1 + *s * *a * (ebx - 1)/ *b),2))/(*a * ebx/(1 + *s * *a * (ebx - 1)/ *b)) * censor[n];
		sumb = sumb + pow((1 + *s * *a * (ebx - 1)/ *b),((-1/ *s) - 1)) * ((-1/ *s) * (*s * *a * (ebx * x[n])/ *b - *s * *a * (ebx - 1)/pow(*b,2)))/pow((1 + *s * *a * (ebx - 1)/ *b),(-1/ *s));
		sums = sums - (*a * ebx * (*a * (ebx - 1)/ *b)/pow((1 + *s * *a * (ebx - 1)/ *b),2)/(*a * ebx/(1 + *s * *a * (ebx - 1)/ *b))) * censor[n];
		sums = sums + (pow((1 + *s * *a * (ebx - 1)/ *b),((-1/ *s) - 1)) * ((-1/ *s) * (*a * (ebx - 1)/ *b)) + pow((1 + *s * *a * (ebx - 1)/ *b),(-1/ *s)) * (log((1 + *s * *a * (ebx - 1)/ *b)) * (1/pow(*s,2))))/pow((1 + *s * *a * (ebx - 1)/ *b),(-1/ *s));
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sums;
}
void glm2(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; double sums;
	while ( --n >= 0){
		suma = suma + (exp(*b * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b) - *a * exp(*b * x[n]) * (*s * (exp(*b * x[n]) - 1)/ *b)/pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),2))/(*c + *a * exp(*b * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b)) * censor[n];
		suma = suma + exp(- *c * x[n]) * (pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),((-1/ *s) - 1)) * ((-1/ *s) * (*s * (exp(*b * x[n]) - 1)/ *b)))/(exp(- *c * x[n]) * pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(- 1/ *s)));
		sumb = sumb + (*a * (exp(*b * x[n]) * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b) - *a * exp(*b * x[n]) * (*s * *a * (exp(*b * x[n]) * x[n])/ *b - *s * *a * (exp(*b * x[n]) - 1)/pow( *b ,2))/pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),2))/(*c + *a * exp(*b * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b)) * censor[n];
		sumb = sumb + exp(- *c * x[n]) * (pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s) - 1) * ((-1/ *s) * (*s * *a * (exp(*b * x[n]) * x[n])/ *b - *s * *a * (exp(*b * x[n]) - 1)/pow(*b,2))))/(exp(- *c * x[n]) * pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s)));
		sumc = sumc + 1/(*c + *a * exp(*b * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b)) * censor[n];
		sumc = sumc - (exp(- *c * x[n]) * x[n] * pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s))/(exp(- *c * x[n]) * pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s))));
		sums = sums -(*a * exp(*b * x[n]) * (*a * (exp(*b * x[n]) - 1)/ *b)/pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),2)/(*c + *a * exp(*b * x[n])/(1 + *s * *a * (exp(*b * x[n]) - 1)/ *b))) * censor[n];
		sums = sums + exp(- *c * x[n]) * (pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),((-1/ *s) - 1)) * ((-1/ *s) * (*a * (exp(*b * x[n]) - 1)/ *b)) + pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s)) * (log((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b)) * (1/ pow(*s,2))))/(exp(- *c * x[n]) * pow((1 + *s * *a * (exp(*b * x[n]) - 1)/ *b),(-1/ *s)));
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc; ans[3] = sums;
}
/****************** optimized versions *****************/
