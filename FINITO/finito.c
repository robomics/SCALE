#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

*	low is the number of nonzeros so that rows with few nonzeros ae excluded; on exit rcontains the actual value used (llow is a pointer to low)

*	maxiter is the maximum number of iterations allowed for one set of perc and perc1 values

*	totalIt is the maximum total number of iterations allowed; on exit contains the total number of iterations
***********************************************************************************************************************/

//      this function is used to sort the number of nonzeros in rows
int cmpfunc (const void * a, const void * b) {
   if ( *(int*)a < *(int*)b ) return -1;
   if ( *(int*)a > *(int*)b ) return 1;
   return(0);
}

// multithreaded matrix-vector multiplication
void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,int *llow, int maxiter, double del, int *totIter, int threads, unsigned int k, double *pppp, int width)
{  
	long p;
	unsigned int n0, lind;
	int *bad;
	double *row, *col, *s, *one, *dr, *dc;	

//	now we can allocate the required arrays
	double *current = (double *) malloc(k*sizeof(double));
	row = (double *) malloc(k*sizeof(double));
	col = (double *) malloc(k*sizeof(double));
	dr = (double *) malloc(k*sizeof(double));
	dc = (double *) malloc(k*sizeof(double));
	bad = (int *) malloc(k*sizeof(int));
	int *bad0 = (int *) malloc(k*sizeof(int));
	one = (double *) malloc(k*sizeof(double));
	s = (double *) malloc(k*sizeof(double));
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0; p<threads;p++) space[p] = (double *) malloc(k*sizeof(double));

//	for convenience and in order not to change the existing code
        int low, low0;
	double perc = *pppp;

//      New percentage policy   
        bool conv = false, div = false;
	double low_conv = 1000, low_div = 0;
	double *b_conv = (double *) malloc(k*sizeof(double));
	double *b0 = (double *) malloc(k*sizeof(double));
	int *bad_conv = (int *) malloc(k*sizeof(int));
	double ber_conv;
	bool yes = true;
	double erez = 1.0e-4;

	for (p=0;p<k;p++) bad[p] = 0;
	for (p=0;p<k;p++) bad0[p] = 0;
	for (p=0;p<k;p++) one[p] = 1.0;

//	find the number of nonzeros in each row
	int *nz = (int *) malloc(k*sizeof(int));
	int *nz0 = (int *) malloc(k*sizeof(int));
	float *diag = (float *) malloc(k*sizeof(float));
	float *diag0 = (float *) malloc(k*sizeof(float));
	for (p=0;p<k;p++) nz[p] = 0;
	for (p=0;p<m;p++) {
		nz[i[p]]+=1.0;
		if (i[p] != j[p]) nz[j[p]]+=1.0;
	}
	n0 = 0;
        for (p=0; p<k;p++) if (nz[p] > 0) nz0[n0++] = nz[p];
	qsort(nz0,n0,sizeof(int),cmpfunc);
	int junk = (int) n0*0.2;
	int bound = nz0[junk];
	junk = (int) n0*perc;
	low = nz0[junk];
	low0 = low;
	int lowMed;

	if (width > 0) {
		for (p=0;p<k;p++) diag[p] = 0;
        	for (p=0;p<m;p++) if (j[p] - i[p] <= width) {
        	        diag[i[p]] += x[p];
        	        diag[j[p]] += x[p];
        	}
        	int nd0 = 0;
        	for (p=0; p<k;p++) if (diag[p] > 0) diag0[nd0++] = diag[p];
        	qsort(diag0,nd0,sizeof(float),cmpfunc);
        	double med = diag0[(int) nd0/2];
		lowMed = (int) (med*0.01);
		for (p=0;p<k;p++) if (diag[p] < lowMed) bad0[p] = 1;
//printf("lowMed = %d; nd0 = %d; med = %g\n",lowMed, nd0,med);
//		FILE *out = fopen("c.diag0","w");
//		for (p=0;p<nd0;p++) fprintf(out,"%g\n",diag0[p]);
//		fclose(out);
	}

//	find rows which are identically 0 and the perc proportion of rows with the lowest row sums and exclude them
	for (p=0;p<k;p++) if (nz[p] < low) bad[p] = 1;
	for (p=0;p<k;p++) bad[p] |= bad0[p];
	for (p=0;p<k;p++) one[p] = 1.0-bad[p];
//	if bad[p] is 1 we remove row p
	for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
	
//	start iterations
//	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
	double ber = 10.0*(1.0+tol);
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

// if converged
            	if (ber < tol) {
//printf("conv: %d %d %d\n",low,iter,all_iters);
			ber_conv = ber;
			low_conv = low;
			yes = true;
			if (low <= low0) break;
                        conv = true;
			for (p=0;p<k;p++) b_conv[p] = b[p];
			for (p=0;p<k;p++) bad_conv[p] = bad[p];
//  did it diverge before?
                        if (div) {
				if (low_conv - low_div <= 1) break;
                                low = (int) (low_conv + low_div)/2;
                        }
//  just halve low
                        else low = (int) low_conv/2;

                        for (p=0;p<k;p++) {
                                one[p] = 1.0;
                                bad[p] = 0;
                        }
                        for (p=0;p<k;p++) if (nz[p] < low) bad[p] = 1;
			for (p=0;p<k;p++) {
				bad[p] |= bad0[p]; 
				one[p] = 1-bad[p];
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

//  diverged
//printf("div: %d %g %d %d\n",low,((double) numBad)/n0,iter,all_iters);
		div = true;
		low_div = low;
//  did it converge before? If it converged for low+1 and diverged for low, use the last converged norm vector 
                if (conv) {
			if (low_conv - low_div <= 1) {
				for (p=0;p<k;p++) b[p] = b_conv[p];
				for (p=0;p<k;p++) bad[p] = bad_conv[p];
				ber = ber_conv;
				break;
			}
//  if it almost converged (only a very small fraction of errors is above tol) remove bad rows and try again
//  with the same low (Erez's trick)			
			else if (((double) numBad)/n0 < erez && yes) {
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
//  if neither, set the new value for low and continue
			else {
				low = (int) (low_div + low_conv)/2;
				yes = true;
			}
                }
//  have never converged before
//  Erez's trick
		else if (((double) numBad)/n0 < erez && yes) {
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
//  otherwise double the low and continue
                else {
			low = 2*low;
			yes = true;
		}

                for (p=0;p<k;p++) bad[p] = 0;
                for (p=0;p<k;p++) if (nz[p] < low) bad[p] = 1.0;
		for (p=0;p<k;p++) {
			bad[p] |= bad0[p];
			one[p] = 1-bad[p];
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
//  find maximum row sums error (all "good" rows should sum to 1)
	double err = 0;
	for (p=0;p<k;p++) {
		if (bad[p] == 1) continue;
		if (err < fabs((double) (col[p]*b[p] - one[p]))) err = fabs((double) (col[p]*b[p] - one[p]));
	}
//  to report statistics
	junk = 0;
	while (nz0[junk] < low_conv) junk++;
	perc = ((double) junk)/n0;
	if (perc > 0.2) for (p=0;p<k;p++) b[p] = NAN;
	report[all_iters+1] = ber;
	report[all_iters+2] = err;
	*llow = low_conv;
	*pppp = perc;
	*totIter = all_iters;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

