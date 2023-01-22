#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>

struct th2 {
	int *i;
	int *j;
	float *x;
	float *v;
	float *res;
	long m;
};

void *Mul(void *threadid) {
	struct th2 *a = (struct th2 *) threadid;
	int *i = a->i;
	int *j = a->j;
	float *x = a->x;
	float *v = a->v;
	float *res = a->res;
	long m = a->m;
	long p;
//printf("m = %ld; i = %ld; j = %ld; x = %ld; v = %ld; res = %ld\n",m,(long) i,(long) j,(long) x,(long) v,(long) res);
	for (p=0;p<m;p++) {
               res[i[p]] += x[p]*v[j[p]];
               res[j[p]] += x[p]*v[i[p]];
	}
//	pthread_exit(NULL);
}

void utmvMul(int *i,int *j,float *x,long m,float *v,int k,float *res) {
	long m1 = (long) m/2;
	long m2 = m - m1;
	long p;
	struct th2 a1,a2;
	float *res1 = (float *) malloc(k*sizeof(float));
	float *res2 = (float *) malloc(k*sizeof(float));
	for (p=0;p<k;p++) res1[p] = 0;
	for (p=0;p<k;p++) res2[p] = 0;
	a1.m = m1;
	a1.i = i;
	a1.j = j;
	a1.x = x;
	a1.v = v;
	a1.res = res1;
	a2.m = m2;
	a2.i = i+m1;
	a2.j = j+m1;
	a2.x = x+m1;
	a2.v = v;
	a2.res = res2;

	pthread_t threads[2];
	int rc1 = pthread_create(&threads[0], NULL, Mul, (void *) &a1);
	int rc2 = pthread_create(&threads[1], NULL, Mul, (void *) &a2);
	pthread_join(threads[0], NULL);
	pthread_join(threads[1], NULL);
	for (p=0;p<k;p++) res[p] = res1[p]+res2[p];
	free(res1);
	free(res2);
}

