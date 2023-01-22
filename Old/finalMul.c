void utmvMul(int *i,int *j,float *x,long m,float *v,int k,float *res) {
	long p;
	for (p=0;p<k;p++) res[p] = 0;
	for (p=0;p<m;p++) {
               res[i[p]] += x[p]*v[j[p]];
               res[j[p]] += x[p]*v[i[p]];
	}
}

