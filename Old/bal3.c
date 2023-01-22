#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <real.h>

void utmvMul(int *i,int *j,real *x,long m,real *v,int k,real *res);

//	this function is used to sort rows sums
int cmpfunc (const void * a, const void * b) {
   if ( *(real*)a < *(real*)b ) return -1;
   if ( *(real*)a > *(real*)b ) return 1;
   return(0);
}

/********************************************************************************************************************
*
*	m is the number of nonzero elements

*	i and j are 0-based arrays each containing the row and column indices of the nonzero bins

*	x is an array containing the nonzero matrix entries

*
*	i, j, and x define the upper triangle of the (squarre symmetric) matrix
*

*	on exit b will hold the scaling vector, i.e. by multiplying the rows and columns of the original matrix 
	by b we get the scaled matrix; 

*	on exit report contains the relative error after each iteration

*	verb indicates whether the function needs to output yhe progress; 1 means report, 0 means run silent

*	tol is the desired relative error

*	perc is the percentage of low rows sums to be excluded (i.e. 0.01, 0.02, etc.)

*	maxiter is the maximum number of iterations allowed

*	del and trial are for determining that the convergence is too slow and early termination (before maxiter iteration): if 
	the relative error decreased by less than del for trials consecuitive iterations the call is terminated and -iter is 
	returned (where iter is the number of iterations); calling function can check the return value and know whether convergence 
	was reached
***********************************************************************************************************************/

int bal(long m,int *i,int *j,real *x, real *d, real report[][2],int verb, double tol,double *pppp, int maxiter, real del)
{  
	long p;
	int k,n,n0;
	real high, low;
//	will be allocated so need to be freed
	int *bad, lind, hind;
	real *r, *r0, *s, *one;
	double perc = *pppp;

//	find the matrix dimensions
	k = 0;
	for (p=0; p<m;p++) if (j[p] > k) k=j[p];
	k++;
	real *current = (real *) malloc(k*sizeof(real));
	r = (real *) malloc(k*sizeof(real));
	r0 = (real *) malloc(k*sizeof(real));
	real *rr = (real *) malloc(k*sizeof(real));
	bad = (int *) malloc(k*sizeof(int));
	one = (real *) malloc(k*sizeof(real));
	s = (real *) malloc(k*sizeof(real));

	for (p=0;p<k;p++) one[p] = 1.0;
	for (p=0;p<k;p++) bad[p] = 0;

//	find rows sums
	utmvMul(i,j,x,m,one,k,rr);
//	find relevant percentiles
	for (p=0; p<k;p++) r[p] = rr[p];
	for (p=0; p<k;p++) r0[p] = r[p];
	qsort(r0,k,sizeof(real),cmpfunc);
	n0 = 0;
	for (p=0;p<k;p++) if (r0[p] == 0) n0++;
	lind = n0-1 + (int)(((real)(k-n0))*perc+0.5);

	if (lind < 0) lind = 0;
	low = r0[lind];

//	find the "bad" rows and exclude them
	for (p=0;p<k;p++) if (r[p] <= low)  bad[p] = 1;

	for (p=0;p<k;p++) d[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
	
//	start iterations
	real ber = 10.0*(1.0+tol);
	real err = 10.0*(1.0+tol);
	int iter=0;
	int fail;
	for (p=0;p<k;p++) current[p] = d[p];
	while((ber > tol || err > 5.0*tol) && iter < maxiter) {
		fail = 1;
		iter++;
		for (p=0;p<k;p++) if (bad[p] == 1) r[p] = 1.0;
		for (p=0;p<k;p++) s[p] = 1.0/sqrt(r[p]);
		for (p=0;p<k;p++) d[p] *= s[p];
	
		utmvMul(i,j,x,m,d,k,r);
		for (p=0;p<k;p++) r[p] *= d[p];
		for (p=0;p<k;p++) if (bad[p] == 1) r[p] = 1.0;

//	calculate the current relative error
		ber=0;
		for (p=0;p<k;p++) {
			if (bad[p] == 1) continue;
			if (fabs(d[p]-current[p]) > ber) ber = fabs(d[p]-current[p]);
		}
		err = 0;
		for (p=0;p<k;p++) {
			if (bad[p] == 1) continue;
			if (err < fabs(r[p] - 1.0)) err = fabs(r[p] - 1.0);
		}
		report[iter-1][0]=ber;
		report[iter-1][1]=err;
		if (verb) printf("%d: %30.15lf %30.15lf\n",iter,ber,err);

		for (p=0;p<k;p++) current[p] = d[p];
		if (iter < 6) continue;
		for (p=2;p<6;p++) {
			if ((report[iter-p+1][0]*(1.0+del) < report[iter-p][0]) && (report[iter-p+1][1]*(1.0+del) < report[iter-p][1])) fail = 0;
		}
		if (fail == 1) {
			perc = 1.5*perc;
			if (perc < 0.005) perc = 0.005;
			if (verb) printf("new p = %g\n",perc);
			lind = n0-1 + (int)(((real)(k-n0))*perc+0.5);
			if (lind < 0) lind = 0;
			low = r0[lind];
			for (p=0;p<k;p++) if (rr[p] <= low)  bad[p] = 1;
			if (report[iter-1][0] < report[iter-5][0]) for (p=0;p<k;p++) d[p] = 1.0 - bad[p];
			else for (p=0;p<k;p++) d[p] *= (1.0 - bad[p]);
			for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
//			utmvMul(i,j,x,m,d,k,r);
			ber = 10.0*(1.0+tol);
			err = 10.0*(1.0+tol);
			iter=0;
			fail = 1;
			for (p=0;p<k;p++) current[p] = d[p];
		}
	}

	report[maxiter+1][0]=ber;
	report[maxiter+1][1]=err;
	(*pppp) = perc;

	for (p=0;p<k;p++) if (bad[p] == 1) d[p] = NAN;
	free(bad);
	free(one);
	free(current);
	free(r);
	free(s);
	free(r0);
	
	return(iter);	
}

