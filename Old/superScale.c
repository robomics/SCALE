#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);

//	this function is used to sort rows sums; may use C++ std::sort instead in the future
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

*	perc is the proportion of low rows sums to be excluded; on exit contains the actual value used (pppp is a pointer to perc)

*	maxiter is the maximum number of iterations allowed for one set of perc and perc1 values

*	totalIt is the maximum total number of iterations allowed; on exit contains the total number of iterations

*	del - if the relative error decreased by less than del  fter 5 iterations we conclude that there is no convergence or it is too slow and then perc is increased accordingly
*
***********************************************************************************************************************/

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,double *pppp, int maxiter, double del, int *totIter, int threads, unsigned int k)
{  
	long p;
	unsigned int n0, lind;
	double low;
	int *bad;
	double *row, *col, *r0, *s, *one, *dr, *dc;	

//	now we can allocate the required arrays
	double *current = (double *) malloc(k*sizeof(double));
	row = (double *) malloc(k*sizeof(double));
	col = (double *) malloc(k*sizeof(double));
	dr = (double *) malloc(k*sizeof(double));
	dc = (double *) malloc(k*sizeof(double));
	r0 = (double *) malloc(k*sizeof(double));
	bad = (int *) malloc(k*sizeof(int));
	one = (double *) malloc(k*sizeof(double));
	s = (double *) malloc(k*sizeof(double));
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0; p<threads;p++) space[p] = (double *) malloc(k*sizeof(double));

//	for convenience and in order not to change the existing code
        double perc = *pppp;

//      New percentage policy   
        double conv_perc, div_perc;
        bool conv = false, div = false;
	double low_conv = 1000, low_div = 0;
	double *b_conv = (double *) malloc(k*sizeof(double));
	double *b0 = (double *) malloc(k*sizeof(double));
	int *bad_conv = (int *) malloc(k*sizeof(int));
	double ber_conv;
	bool yes = true;

	for (p=0;p<k;p++) bad[p] = 0;
	for (p=0;p<k;p++) one[p] = 1.0;

//	find the number of nonzeros in each row
	double *nz = (double *) malloc(k*sizeof(double));
	for (p=0;p<k;p++) nz[p] = 0;
	for (p=0;p<m;p++) {
		nz[i[p]]+=1.0;
		if (i[p] != j[p]) nz[j[p]]+=1.0;
	}
	n0 = 0;
        for (p=0; p<k;p++) if (nz[p] > 0) r0[n0++] = nz[p];
        qsort(r0,n0,sizeof(double),cmpfunc);
	double abc = n0*perc+0.5;
        low = 1.0;
	int junk = (int) n0*0.2;
	double bound = r0[junk];

//	find rows which are identically 0 and the perc proportion of rows with the lowest row sums and exclude them
	for (p=0;p<k;p++) if (nz[p] < low) bad[p] = 1;
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
		double temp1;
		int numBad=0;

		for (p=0;p<k;p++) {
			if (bad[p] == 1) continue;
			temp1 = fabs((b[p]-current[p])/(b[p]+current[p]));
			if (temp1 > ber) ber = temp1;
			if (temp1 > tol) numBad++;
		}
		report[all_iters-1]=ber;
		allIters[all_iters-1]=iter;

		for (p=0;p<k;p++) b0[p] = current[p];
		for (p=0;p<k;p++) current[p] = b[p];

            	if (ber < tol) {
			printf("conv: %g %d %d\n",low,iter,all_iters);
			yes = true;
			if (low == 1) break;
                        conv = true;
			for (p=0;p<k;p++) b_conv[p] = b[p];
			for (p=0;p<k;p++) bad_conv[p] = bad[p];
			ber_conv = ber;
			low_conv = low;
                        if (div) {
				if (low_conv - low_div <= 1) break;
                                low = (int) (low_conv + low_div)/2;
                        }
                        else low = (int) low_conv/2;

                        for (p=0;p<k;p++) {
                                one[p] = 1.0;
                                bad[p] = 0;
                        }
                        for (p=0;p<k;p++) {
                                if (nz[p] < low) {
                                        bad[p] = 1;
                                        one[p] = 0;
                                }
                        }
                        iter = 0;
                        ber = 10.0;
			for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
			for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
			utmvMul(i,j,x,m,dc,k,row,threads,space);
			for (p=0;p<k;p++) row[p] *= dr[p];
			continue;
		}

                if (iter <=  5) continue;
//	check whether convergence rate is satisfactory
		if ((report[all_iters-1]*(1.0+del) < report[all_iters-6]) && (iter < maxiter)) continue;
//	the scaling vector does not seen to converge well enough - increase perc and continue

		div = true;
		printf("div: %g %g %d %d\n",low,((double) numBad)/n0,iter,all_iters);
		low_div = low;
                if (conv) {
			if (low_conv - low_div <= 1) {
				for (p=0;p<k;p++) b[p] = b_conv[p];
				for (p=0;p<k;p++) bad[p] = bad_conv[p];
				ber = ber_conv;
				break;
			}
			else if (((double) numBad)/n0 < 1.0e-5 && yes) {
				for (p=0;p<k;p++) {
					if (bad[p] == 1) continue;
					temp1 = fabs((b[p]-b0[p])/(b[p]+b0[p]));
					if (temp1 > tol) {
						bad[p] = 1;
						one[p] = 0;
					}
				}
				yes = false;
				goto next;
			}
			else {
				low = (int) (low_div + low_conv)/2;
				yes = true;
			}
                }
		else if (((double) numBad)/n0 < 1.0e-5 && yes) {
			for (p=0;p<k;p++) {
				if (bad[p] == 1) continue;
				temp1 = fabs((b[p]-b0[p])/(b[p]+b0[p]));
				if (temp1 > tol) {
					bad[p] = 1;
					one[p] = 0;
				}
			}
			yes = false;
			goto next;
		}
                else {
			low = 2*low;
			yes = true;
		}

                for (p=0;p<k;p++) {
                        bad[p] = 0;
                        one[p] = 1.0;
                }
                for (p=0;p<k;p++) {
                        if (nz[p] < low) {
                                bad[p] = 1.0;
                                one[p] = 0;
                        }
                }
next:		ber = 10.0;
                iter=0;
		for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
		for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
		utmvMul(i,j,x,m,dc,k,row,threads,space);
		for (p=0;p<k;p++) row[p] *= dr[p];

//	if perc reached upper bound or the total number of iterationbs is too high, exit
		if (low > bound) break;
		if (all_iters > *totIter) break;
	}

//	find the final error in row sums
	utmvMul(i,j,x,m,b,k,col,threads,space);
	err = 0;
	for (p=0;p<k;p++) {
		if (bad[p] == 1) continue;
		if (err < fabs((double) (col[p]*b[p] - one[p]))) err = fabs((double) (col[p]*b[p] - one[p]));
	}
	junk = 0;
	while (r0[junk] < low) junk++;
	perc = ((double) junk)/n0;
	report[all_iters+1] = ber;
	report[all_iters+2] = err;
	*pppp = perc;
	*totIter = all_iters;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

