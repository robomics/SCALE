#include <real.h>

void utmvMul(int *i,int *j,real *x,long m,real *v,int k,real *res) {
	long p;
	for (p=0;p<k;p++) res[p] = 0;
	for (p=0;p<m;p++) {
               res[i[p]] += x[p]*v[j[p]];
               res[j[p]] += x[p]*v[i[p]];
	}
}

