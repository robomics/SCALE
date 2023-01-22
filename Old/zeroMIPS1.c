#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void utmvMul(int c, int **i,int **j,float **x,int *m,float *v,int k,float *es);

//	this function is used to sort rows sums
int cmpfunc (const void * a, const void * b) {
   if ( *(float*)a < *(float*)b ) return -1;
   if ( *(float*)a > *(float*)b ) return 1;
   return(0);
}

/********************************************************************************************************************
*
*	This function allows more that 2^31 - 1 nonzero entries. It acceptc a list of c arrays where array i contains m[i] elements
*

*	c is the number of arrays

*	m is array containing the number of elements of the c arrays

*	i and j are lists of c 0-based arrays each containing the row and column indices of the nonzero bins

*	x is a list of c arrays containing the nonzero matrix entries

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

int scale(int c, int *m,int **i,int **j,float **x, float *z,float *b, float *report,
int verb, float tol,float perc,float perc1, int maxiter, float del,int trials,int zerodiag) 
{  
	int k,p,ic,n;
	float high, low;
//	will be allocated so need to be freed
	int *bad, lind, hind;
	float *row, *col, *dr, *dc, *r0, *s, *one, *zz;	

//	find the matrix dimensions
	k = 0;
	for (ic=0;ic<c;ic++) for (p=0; p<m[ic];p++) if (j[ic][p] > k) k=j[ic][p];
	k++;
	float *current = (float *) malloc(k*sizeof(float));
	row = (float *) malloc(k*sizeof(float));
	col = (float *) malloc(k*sizeof(float));
	dr = (float *) malloc(k*sizeof(float));
	dc = (float *) malloc(k*sizeof(float));
	r0 = (float *) malloc(k*sizeof(float));
	bad = (int *) malloc(k*sizeof(int));
	int *bad1 = (int *) malloc(k*sizeof(int));
	one = (float *) malloc(k*sizeof(float));
	s = (float *) malloc(k*sizeof(float));
	zz = (float *) malloc(k*sizeof(float));
	int l = 0;
	for (p=0;p<k;p++) {
		if (isnan(z[p])) continue; 
		if (z[p] > 0) zz[l++] = z[p];
	}
	qsort(zz,l,sizeof(float),cmpfunc);
	lind = (int)(((float) l)*perc1+0.5);
	hind = (int)(((float) l)*(1.0-perc1)+0.5);
	if (lind < 0) lind = 0;
	if (hind >= l) hind = l-1;
	low = zz[lind];
	high = zz[hind];
	free(zz);
	for (p=0;p<k;p++) if (z[p] > 0 && (z[p] < low || z[p] > high)) z[p] = NAN;

        for (p=0;p<k;p++) one[p] = 1.0;
        for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;

	for (p=0;p<k;p++) bad[p] = 0;
	if (zerodiag == 1) for (ic=0;ic<c;ic++) for (p=0;p<m[ic];p++) if (i[ic][p] == j[ic][p] && x[ic][p] == 0) bad[i[ic][p]] = 1;

//	find rows sums
	utmvMul(c,i,j,x,m,one,k,row);

//	find relevant percentiles
	for (p=0; p<k;p++) r0[p] = row[p];
	qsort(r0,k,sizeof(float),cmpfunc);
	n = 0;
	for (p=0;p<k;p++) if (r0[p] == 0) n++;
	lind = n-1 + (int)(((float)(k-n))*perc+0.5);
        hind = n-1 + (int)(((float)(k-n))*(1.0-0.1*perc)+0.5);

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
	float ber = 10.0*(1.0+tol);
	float err = 10.0*(1.0+tol);
	int iter=0;
	int stuck=0;
	for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
	while((ber > tol || err > 5.0*tol) && iter++ < maxiter) {
		for (p=0;p<k;p++) if (bad1[p] == 1) row[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/row[p];
		for (p=0;p<k;p++) dr[p] *= s[p];
	
		utmvMul(c,i,j,x,m,dr,k,col);
		for (p=0;p<k;p++) col[p] *= dc[p];
		for (p=0;p<k;p++) if (bad1[p] == 1) col[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/col[p];
		for (p=0;p<k;p++) dc[p] *= s[p];

		utmvMul(c,i,j,x,m,dc,k,row);
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
			utmvMul(c,i,j,x,m,b,k,col);
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

	utmvMul(c,i,j,x,m,b,k,col);
	err = 0;
	for (p=0;p<k;p++) {
		if (bad1[p] == 1) continue;
		if (err < fabs(col[p]*b[p] - z[p])) err = fabs(col[p]*b[p] - z[p]);
	}
	printf("the error is %30.15f\n",err);

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
	
	if (ber > tol || err > 5*tol) return(-iter);	
	else return(iter);	
}

