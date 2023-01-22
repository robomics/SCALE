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
***********************************************************************************************************************/

int scale(long m,int *i,int *j,real *x, real *z,real *b, real *report,int verb, double tol,double *pppp,double *pppp1, int maxiter, int zerodiag, double del, double dp, double dp1)
{  
	long p;
	int k,n0;
	real low, zhigh, zlow;
//	will be allocated so need to be freed
	int *bad, lind, hind;
	real *row, *col, *r0, *s, *one, *zz, *row0, *dr, *dc;	

//	find the matrix dimensions
	k = 0;
	for (p=0; p<m;p++) if (j[p] > k) k=j[p];
	k++;
	real *current = (real *) malloc(k*sizeof(real));
	row = (real *) malloc(k*sizeof(real));
	row0 = (real *) malloc(k*sizeof(real));
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

        double perc = *pppp;
        double perc1 = *pppp1;

	for (p=0;p<k;p++) {
		if (isnan(z[p])) continue; 
		if (z[p] > 0) zz[l++] = z[p];
	}
	qsort(zz,l,sizeof(real),cmpfunc);
	lind = (int)(l*perc1+0.5);
	hind = (int)(l*(1.0-perc1)+0.5);
	if (lind < 0) lind = 0;
	if (hind >= l) hind = l-1;
	zlow = zz[lind];
	zhigh = zz[hind];
	for (p=0;p<k;p++) if (z[p] > 0 && (z[p] < zlow || z[p] > zhigh)) z[p] = NAN;

	for (p=0;p<k;p++) one[p] = 1.0;
        for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;

	if (zerodiag == 1) {
		for (p=0;p<k;p++) bad[p] = 1;
		for (p=0;p<m;p++) if (i[p] == j[p]) bad[i[p]] = 0;
	} else {
		for (p=0;p<k;p++) bad[p] = 0;
	}

//	find rows sums
	utmvMul(i,j,x,m,one,k,row0);
//	find relevant percentiles
	for (p=0; p<k;p++) row[p] = row0[p];
	n0 = 0;
        for (p=0; p<k;p++) if (row0[p] > 0) r0[n0++] = row0[p];
        qsort(r0,n0,sizeof(real),cmpfunc);
	lind =  (int)(n0*perc+0.5);
        if (lind < 0) lind = 0;
        low = r0[lind];
//	find the "bad" rows and exclude them
	for (p=0;p<k;p++) {
		if ((row[p] <= low  && z[p] > 0) || isnan(z[p])) {
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
	for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
      	int fail;
        int nerr = 0;
        double errors[10000];

	while((ber > tol || err > 5.0*tol) && iter < maxiter) {
		iter++;
		fail = 1;
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
			errors[nerr++] = err;
			if (verb) printf("the error is %30.15f\n",err);
		}

		for (p=0;p<k;p++) current[p] = b[p];

            	if ((ber < tol) && (nerr < 2 || (nerr >= 2 && errors[nerr-1] < 0.5*errors[nerr-2])))  continue;
                if (iter > 5) {
                        for (p=2;p<=5;p++) if (report[iter-p+1]*(1.0+del) < report[iter-p]) fail = 0;
                        if (nerr >= 2 && errors[nerr-1] > 0.75*errors[nerr-2]) fail = 1;
                        if (fail == 1) {
                                perc += dp;
                                perc1 += dp1;
                                nerr = 0;
				if (verb == 1) printf("perc = %g and perc1 = %g\n",perc,perc1);
                                lind =  (int)(n0*perc+0.5);
                                low = r0[lind];
				lind = (int)(l*perc1+0.5);
			        hind = (int)(l*(1.0-perc1)+0.5);
        			if (lind < 0) lind = 0;
        			if (hind >= l) hind = l-1;
        			zlow = zz[lind];
        			zhigh = zz[hind];
        			for (p=0;p<k;p++) if (z[p] > 0 && (z[p] < zlow || z[p] > zhigh)) z[p] = NAN;
				for (p=0;p<k;p++) {
					if ((row0[p] <= low  && z[p] > 0) || isnan(z[p])) {
						bad[p] = 1;
						bad1[p] = 1;
						one[p] = 0;
						z[p] = 1.0;
                			}
        			}


                                ber = 10.0*(1.0+tol);
                                err = 10.0*(1.0+tol);
                                if (report[iter-1] > report[iter-5]) {
                                        for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
                                        for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
                                        for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
                                        for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
                                        for (p=0;p<k;p++) row[p] = row0[p];
                                }
                                else {
                                        for (p=0;p<k;p++) dr[p] *= (1.0-bad[p]);
                                        for (p=0;p<k;p++) dc[p] *= (1.0-bad[p]);
                                }
                                iter=0;

                        }
                }
		if (perc > 0.2 || perc1 > 0.1) break;
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
	*pppp = perc;
	*pppp1 = perc1;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

