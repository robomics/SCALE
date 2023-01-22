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

*	z is the "target" vector, i.e. we want rows (and columns) sums to be equal to z

*	on exit b will hold the scaling vector, i.e. by multiplying the rows and columns of the original matrix 
	by b we get the scaled matrix; 

*	on exit report contains the relative error after each iteration

*
*	bekow are arguments having default values
*

*	verb indicates whether the function needs to output yhe progress; 1 means report, 0 means run silent

*	tol is the desired relative error

*	perc is the percentage of low rows sums to be excluded (i.e. 0.01, 0.02, etc.)

*	perc1 is the percentage of low and high values of z to ignore

*	maxiter is the maximum number of iterations allowed

*	del and trial are for determining that the convergence is too slow and early termination (before maxiter iteration): if 
	the relative error decreased by less than del for trials consecuitive iterations the call is terminated and -iter is 
	returned (where iter is the number of iterations); calling function can check the return value and know whether convergence 
	was reached
*
*	Note that making any optional argument negative causes the default value to be used
*
***********************************************************************************************************************/

int scale(long m,int *i,int *j,real *x, real *z,real *b, real *report,int verb, real tol,real perc,real perc1, int maxiter, real del,int trials,int zerodiag) 
{  
	long p;
	int k,n;
	real high, low;
//	will be allocated so need to be freed
	int *bad, lind, hind;
	real *row, *col, *dr, *dc, *r0, *s, *one, *zz;	

//	find the matrix dimensions
	k = 0;
	for (p=0; p<m;p++) if (j[p] > k) k=j[p];
	k++;
	real *current = (real *) malloc(k*sizeof(real));
	row = (real *) malloc(k*sizeof(real));
	col = (real *) malloc(k*sizeof(real));
	dr = (real *) malloc(k*sizeof(real));
	dc = (real *) malloc(k*sizeof(real));
	r0 = (real *) malloc(k*sizeof(real));
	bad = (int *) malloc(k*sizeof(int));
	int *bad1 = (int *) malloc(k*sizeof(int));
	one = (real *) malloc(k*sizeof(real));
	s = (real *) malloc(k*sizeof(real));
	zz = (real *) malloc(k*sizeof(real));
	int l = 0;
	for (p=0;p<k;p++) {
		if (isnan(z[p])) continue; 
		if (z[p] > 0) zz[l++] = z[p];
	}
	qsort(zz,l,sizeof(real),cmpfunc);
	lind = (int)(((real) l)*perc1+0.5);
	hind = (int)(((real) l)*(1.0-perc1)+0.5);
	if (lind < 0) lind = 0;
	if (hind >= l) hind = l-1;
	low = zz[lind];
	high = zz[hind];
	free(zz);
	for (p=0;p<k;p++) if (z[p] > 0 && (z[p] < low || z[p] > high)) z[p] = NAN;

	for (p=0;p<k;p++) one[p] = 1.0;
        for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;

	if (zerodiag == 1) {
		for (p=0;p<k;p++) bad[p] = 1;
		for (p=0;p<m;p++) if (i[p] == j[p]) bad[i[p]] = 0;
	} else {
		for (p=0;p<k;p++) bad[p] = 0;
	}

//	find rows sums
	utmvMul(i,j,x,m,one,k,row);

//	find relevant percentiles
	for (p=0; p<k;p++) r0[p] = row[p];
	qsort(r0,k,sizeof(real),cmpfunc);
	n = 0;
	for (p=0;p<k;p++) if (r0[p] == 0) n++;
	lind = n-1 + (int)(((real)(k-n))*perc+0.5);
        hind = n-1 + (int)(((real)(k-n))*(1.0-0.1*perc)+0.5);

	if (lind < 0) lind = 0;
        if (hind >= k) hind = k-1;
	low = r0[lind];
	high = r0[hind];
	free(r0);

//	find the "bad" rows and exclude them
	for (p=0;p<k;p++) {
		if (((row[p] < low || row[p] > high) && z[p] > 0) || isnan(z[p])) {
//		if ((row[p] < low  && z[p] > 0) || isnan(z[p])) {
			bad[p] = 1;
			z[p] = 1.0;
		}
	}

	for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;
	for (p=0;p<k;p++) bad1[p] = 1.0 - one[p];
	
//	start iterations
//	row is the current rows sum; s is the correction vector to be applied to rows and columns
	real ber = 10.0*(1.0+tol);
	real err = 10.0*(1.0+tol);
	int iter=0;
	int stuck=0;
	for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
	while((ber > tol || err > 5.0*tol) && iter < maxiter) {
		iter++;
		for (p=0;p<k;p++) if (bad1[p] == 1) row[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/row[p];
		for (p=0;p<k;p++) dr[p] *= s[p];
	
		utmvMul(i,j,x,m,dr,k,col);
		for (p=0;p<k;p++) col[p] *= dc[p];
		for (p=0;p<k;p++) if (bad1[p] == 1) col[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/col[p];
		for (p=0;p<k;p++) dc[p] *= s[p];

		utmvMul(i,j,x,m,dc,k,row);
		for (p=0;p<k;p++) row[p] *= dr[p];

		for (p=0;p<k;p++) b[p] = sqrt(dr[p]*dc[p]);
	
//	calculate the current relative error
		ber=0;
		for (p=0;p<k;p++) {
			if (bad1[p] == 1) continue;
			if (fabs(b[p]-current[p]) > ber) ber = fabs(b[p]-current[p]);
		}
		report[iter-1]=ber;
		if (verb) printf("%d: %30.15lf\n",iter,ber);
		if (iter % 10 == 0) {
			utmvMul(i,j,x,m,b,k,col);
			err = 0;
			for (p=0;p<k;p++) {
				if (bad1[p] == 1) continue;
				if (err < fabs(col[p]*b[p] - z[p])) err = fabs(col[p]*b[p] - z[p]);
			}
			if (verb) printf("the error is %30.15f\n",err);
		}

		for (p=0;p<k;p++) current[p] = b[p];
		if (iter < trials+2) continue;
		if (ber > (1.0-del)*report[iter-2]) stuck++; 
		else stuck = 0;
		if (stuck >= trials) break;
	}

	utmvMul(i,j,x,m,b,k,col);
	err = 0;
	for (p=0;p<k;p++) {
		if (bad1[p] == 1) continue;
		if (err < fabs(col[p]*b[p] - z[p])) err = fabs(col[p]*b[p] - z[p]);
	}
	if (verb) printf("the error is %30.15f\n",err);
	report[maxiter+1] = ber;
	report[maxiter+2] = err;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	free(bad);
	free(bad1);
	free(one);
	free(current);
	free(row);
	free(col);
	free(dr);
	free(dc);
	free(s);
	
	if (ber > tol || err > 5.0*tol) return(-iter);	
	else return(iter);	
}

