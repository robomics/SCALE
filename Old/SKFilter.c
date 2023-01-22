#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);

//	this function is used to sort rows sums
int cmpfunc (const void * a, const void * b) {
   if ( *(double*)a < *(double*)b ) return -1;
   if ( *(double*)a > *(double*)b ) return 1;
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

*	b - on exit will hold the scaling vector, i.e. by multiplying the rows and columns of the original matrix 
	by b we get the scaled matrix; 

*	report and allIters - on exit contain the error after each iteration and the accumulative number of iterations

*	tol is the desired relative error

*	maxiter is the maximum number of iterations allowed for one set of perc and perc1 values

*	totalIt is the maximum total number of iterations allowed; on exit contains the total number of iterations

*	del - if the relative error decreased by less than del for 4 consecuitive iterations we conclude that there is no convergence or it is too slow and then perc and perc1 are increased by dp and dp1 respectively
*
***********************************************************************************************************************/

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol, int maxiter, double del, int *totIter, int threads, unsigned int k, int width, double *pperc)
{  
	long p;
	unsigned int n0, lind;
	double low;
	int *bad;
	double *row, *col, *r0, *s, *one, *row0, *dr, *dc;	
	bool div = true;
	double perc = *pperc;

//	now we can allocate the required arrays
	double *current = (double *) malloc(k*sizeof(double));
	row = (double *) malloc(k*sizeof(double));
	row0 = (double *) malloc(k*sizeof(double));
	col = (double *) malloc(k*sizeof(double));
	dr = (double *) malloc(k*sizeof(double));
	dc = (double *) malloc(k*sizeof(double));
	r0 = (double *) malloc(k*sizeof(double));
	bad = (int *) malloc(k*sizeof(int));
	one = (double *) malloc(k*sizeof(double));
	s = (double *) malloc(k*sizeof(double));
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0; p<threads;p++) space[p] = (double *) malloc(k*sizeof(double));

	for (p=0;p<k;p++) bad[p] = 0;
	for (p=0;p<k;p++) one[p] = 1.0;

//	find the number of nonzeros in each row
	double *medians = (double *) malloc(k*sizeof(double));
	for (p=0;p<k;p++) medians[p] = 0;
	for (p=0;p<m;p++) if (j[p] - i[p] <= width) {
		medians[i[p]] += x[p];
		medians[j[p]] += x[p];
	}
//	find the cutoff
	n0 = 0;
        for (p=0; p<k;p++) if (medians[p] > 0) r0[n0++] = medians[p];
        qsort(r0,n0,sizeof(double),cmpfunc);
	double med = r0[(int) n0/2];
	low = med*perc;
printf("n0 = %d; median = %g; low = %g\n",n0,med,low);
	for (p=0;p<k;p++) if (medians[p] < low) bad[p] = 1;
	for (p=0;p<k;p++) one[p] = 1.0-bad[p];
//	if bad[p] is 1 we remove row p
	for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
	
//	start iterations
//	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
	double ber = 10.0*(1.0+tol);
	double err = 10.0*(1.0+tol);
	int iter=0;
	for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
	int all_iters = 0;

	utmvMul(i,j,x,m,dc,k,row,threads,space);
	for (p=0;p<k;p++) row[p] *= dr[p];
	while(ber > tol  && iter < maxiter) {
		iter++;
		all_iters++;
		for (p=0;p<k;p++) dr[p] = (bad[p] == 0)? dr[p]/row[p]:0;

//printf("perc = %g; iter = %d; all_iters = %d; ber = %g\n",perc,iter,all_iters,ber);
//	find column sums and update rows scaling vector
		utmvMul(i,j,x,m,dr,k,col,threads,space);
		for (p=0;p<k;p++) col[p] *= dc[p];
		for (p=0;p<k;p++) dc[p] = (bad[p] == 0)? dc[p]/col[p]:0;

//	find row sums and update columns scaling vector
		utmvMul(i,j,x,m,dc,k,row,threads,space);
		for (p=0;p<k;p++) row[p] *= dr[p];
//	calculate current scaling vector
		for (p=0;p<k;p++) b[p] = sqrt(dr[p]*dc[p]);
//	calculate the current error
		ber=0;
		for (p=0;p<k;p++) {
			if (bad[p] == 1) continue;
			if (fabs((double) ((b[p]-current[p])/(b[p]+current[p]))) > ber) ber = fabs((double) ((b[p]-current[p])/(b[p]+current[p])));
		}
		report[all_iters-1]=ber;
		allIters[all_iters-1]=iter;

		for (p=0;p<k;p++) current[p] = b[p];

            	if (ber < tol) {
			div = false;
			break;
		}

                if (iter <=  5) continue;
//	check whether convergence rate is satisfactory
		if ((report[all_iters-1]*(1.0+del) < report[all_iters-6]) && (iter < maxiter)) continue;
//	the scaling vector does not seen to converge well enough or the row sums erros do not decrease enough - 
//	increase perc and perc1 and continue

		perc *= 1.1;
		double junk = med*perc;
		if ( junk < low+1) {
			junk = low+1;
			perc = junk /med;
		}
		low = junk;
		for (p=0;p<k;p++) if (medians[p] < low) bad[p] = 1;
		for (p=0;p<k;p++) one[p] = 1.0-bad[p];
//		if bad[p] is 1 we remove row p
		for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
		for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
		iter = 0;
printf("current perc = %g and low = %g\n",perc,low);

		if (all_iters > *totIter) break;
	}

//	find the final error in row sums
	utmvMul(i,j,x,m,b,k,col,threads,space);
	err = 0;
	for (p=0;p<k;p++) {
		if (bad[p] == 1) continue;
		if (err < fabs((double) (col[p]*b[p] - one[p]))) err = fabs((double) (col[p]*b[p] - one[p]));
	}
	report[all_iters+1] = ber;
	report[all_iters+2] = err;
	*totIter = all_iters;
	*pperc = perc;
	
	if (div) {
		for (p=0;p<k;p++) b[p] = NAN;
		return(-iter);
	}
	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

