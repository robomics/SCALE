#include <stdlib.h>
#include <math.h>

void utmvMul(int *i,int *j,float *x,long m,double *v,int k,double *res, int threads, double **space);

double ppNormVector(long m,int *ii,int *jj,float *xx, double *b,int k, int threads, double **space) {
	double s1,s2;
	double *one = (double *) malloc(k*sizeof(double));
	double *u = (double *) malloc(k*sizeof(double));
	double *v = (double *) malloc(k*sizeof(double));
	int i;
	long p;
	for (p=0;p<k;p++) if (b[p] == 0) b[p] = NAN;
	for (p=0;p<k;p++) one[p] = isnan(b[p])?0:1;
	utmvMul(ii,jj,xx,m,one,k,u,threads,space);
	s1 = 0;
	for (p=0;p<k;p++) s1 += u[p]*one[p];
	for (p=0;p<k;p++) v[p] = isnan(b[p])?0:b[p];
	utmvMul(ii,jj,xx,m,v,k,u,threads,space);
	s2 = 0;
	for (p=0;p<k;p++) s2 += u[p]*v[p];
	double s = sqrt(s2/s1);
	for (p=0;p<k;p++) if (!isnan(b[p])) b[p] = s/b[p];

	return(s1);
}




