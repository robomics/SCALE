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

*	z is the "target" vector, i.e. we want rows (and columns) sums to be equal to z

*	b - on exit will hold the scaling vector, i.e. by multiplying the rows and columns of the original matrix 
	by b we get the scaled matrix; 

*	report and allIters - on exit contain the error after each iteration and the accumulative number of iterations

*	tol is the desired relative error

*	zerodiag: if 1, remove rows with 0 at the diagonal; if 0 do not do this

*	perc is the proportion of low rows sums to be excluded; on exit contains the actual value used (pppp is a pointer to perc)

*	dp - by how much to increase perc each time

*	perc1 is the proportion of low and high values of z to be excluded; on exit contains the actual value used (pppp1 is a pointer to perc1)

*	dp1 - by how much to increase perc1 each time

*	maxiter is the maximum number of iterations allowed for one set of perc and perc1 values

*	totalIt is the maximum total number of iterations allowed; on exit contains the total number of iterations

*	del - if the relative error decreased by less than del for 4 consecuitive iterations we conclude that there is no convergence or it is too slow and then perc and perc1 are increased by dp and dp1 respectively
*
***********************************************************************************************************************/

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,double *pppp, int maxiter, double del, double dp, int *totIter, int threads, unsigned int k)
{  
	long p;
	int n0;
	double low;
	int *bad, lind;
	double *row, *r0, *one, *row0;

//	now we can allocate the required arrays
	double *current = (double *) malloc(k*sizeof(double));
	row = (double *) malloc(k*sizeof(double));
	row0 = (double *) malloc(k*sizeof(double));
	r0 = (double *) malloc(k*sizeof(double));
	bad = (int *) malloc(k*sizeof(int));
	one = (double *) malloc(k*sizeof(double));
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0; p<threads;p++) space[p] = (double *) malloc(k*sizeof(double));

        double perc = *pppp;

	for (p=0;p<k;p++) one[p] = 1.0;
	for (p=0;p<k;p++) bad[p] = 0;

//	find rows sums
	double *nz = (double *) malloc(k*sizeof(double));
	for (p=0;p<k;p++) nz[p] = 0;
	for (p=0;p<m;p++) {
		nz[i[p]]+=1.0;
		if (i[p] != j[p]) nz[j[p]]+=1.0;
	}
//	utmvMul(i,j,x,m,one,k,row0,threads,space);
//	find relevant percentiles
	n0 = 0;
        for (p=0; p<k;p++) if (nz[p] > 0) r0[n0++] = nz[p];
        qsort(r0,n0,sizeof(double),cmpfunc);
	lind =  (int)(n0*perc+0.5);
        if (lind < 0) lind = 0;
        low = r0[lind];
//	find rows which are identically 0 and the perc proportion of rows with the lowest row sums and exclude them
	for (p=0;p<k;p++) if (nz[p] < low) bad[p] = 1;
//	if bad[p] is 1 we remove row p
	utmvMul(i,j,x,m,one,k,row0,threads,space);
	for (p=0; p<k;p++) row[p] = row0[p];

	for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) if (bad[p] == 1) row[p] = 1.0;
	for (p=0;p<k;p++) b[p] = one[p]/sqrt(row[p]);

//	start iterations
//	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
	double ber = 10.0*(1.0+tol);
	double err = 10.0*(1.0+tol);
	int iter=0;
	for (p=0;p<k;p++) current[p] = b[p];
      	int fail;   // checks whether the convergence rate seems to be good enough; 0 if yes, 1 if no
        int nerr = 0;
        double errors[10000];
	int all_iters = 0;

	while((ber > tol || err > 5.0*tol) && iter < maxiter) {
		iter++;
		all_iters++;
		fail = 1;

		utmvMul(i,j,x,m,b,k,row,threads,space);
		for (p=0;p<k;p++) row[p] *= b[p];
		for (p=0;p<k;p++) if (bad[p] == 1) row[p] = 1.0;
		for (p=0;p<k;p++) b[p] *= one[p]/sqrt(row[p]);

//	calculate the current error
		ber=0;
		for (p=0;p<k;p++) {
			if (bad[p] == 1) continue;
			if (fabs((double) (b[p]-current[p])) > ber) ber = fabs((double) (b[p]-current[p]));
		}
		report[all_iters-1]=ber;
		allIters[all_iters-1]=iter;

//	since calculating the error in row sums requires matrix-vector multiplication we are are doing this every 10 
//	iterations
		if (iter % 10 == 0) {
			utmvMul(i,j,x,m,b,k,row,threads,space);
			err = 0;
			for (p=0;p<k;p++) {
				if (bad[p] == 1) continue;
				if (err < fabs((double) (row[p]*b[p] - one[p]))) err = fabs((double) (row[p]*b[p] - one[p]));
			}
			errors[nerr++] = err;
		}

		for (p=0;p<k;p++) current[p] = b[p];

//	check whether convergence rate is satisfactory
//	if less than 5 iterations (so less than 5 errors) and less than 2 row sums errors, there is nothing to chek
            	if ((ber < tol) && (nerr < 2 || (nerr >= 2 && errors[nerr-1] < 0.5*errors[nerr-2])))  continue;
//	otherwise check
                if (iter > 5) {
                        for (p=1;p<=5;p++) if (report[all_iters-p]*(1.0+del) < report[all_iters-p-1]) fail = 0;
                        if (nerr >= 2 && errors[nerr-1] > 0.75*errors[nerr-2]) fail = 1;
			if (iter >= maxiter) fail = 1;
//	the scaling vector does not seen to converge well enough or the row sums erros do not decrease enough - 
//	increase perc and perc1 and continue
                        if (fail == 1) {
                                perc += dp;
                                nerr = 0;
                                lind =  (int)(n0*perc+0.5);
                                low = r0[lind];
				for (p=0;p<k;p++) {
					if (nz[p] < low) {
						bad[p] = 1;
						one[p] = 0;
                			}
        			}

                                ber = 10.0*(1.0+tol);
                                err = 10.0*(1.0+tol);
//	if the current error is larger than 5 iteration ago start from scratch, otherwise continue from the current 
//	position 
                                if (report[all_iters-1] > report[all_iters-6]) {
                                        for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
                                        for (p=0;p<k;p++) row0[p] *= one[p];
                                        for (p=0;p<k;p++) row[p] = row0[p];
                                        for (p=0;p<k;p++) current[p] = b[p];
                                }
                                else for (p=0;p<k;p++) b[p] *= (1.0-bad[p]);
                                iter=0;

                        }
                }

//	if perc or perc1 reached upper bound or the total number of iterationbs is too high, exit
		if (perc > 0.2) break;
		if (all_iters > *totIter) break;
	}

//	find the final error in row sums
	utmvMul(i,j,x,m,b,k,row,threads,space);
	err = 0;
	for (p=0;p<k;p++) {
		if (bad[p] == 1) continue;
		if (err < fabs((double) (row[p]*b[p] - one[p]))) err = fabs((double) (row[p]*b[p] - one[p]));
	}
	report[all_iters+1] = ber;
	report[all_iters+2] = err;
	*pppp = perc;
	*totIter = all_iters;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

