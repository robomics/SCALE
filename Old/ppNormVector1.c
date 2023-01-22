#include <stdlib.h>
#include <math.h>

void utmvMul(int *i,int *j,float *x,long m,float *v,int k,float *res, int threads, float **space);

double ppNormVector(long m,int *ii,int *jj,float *xx, float *b,int k, int threads, float **space) {
	double s1,s2;
	float *one = (float *) malloc(k*sizeof(float));
	float *u = (float *) malloc(k*sizeof(float));
	float *v = (float *) malloc(k*sizeof(float));
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
	float s = sqrt(s2/s1);
	for (p=0;p<k;p++) if (!isnan(b[p])) b[p] = s/b[p];

	return(s1);
}




