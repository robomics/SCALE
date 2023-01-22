#include <stdlib.h>
#include <math.h>

void utmvMul(int *i,int *j,float *x,long m,float *v,int k,float *res, int threads, float **space);

//	this function is used to sort rows sums
int cmpfunc (const void * a, const void * b) {
   if ( *(float*)a < *(float*)b ) return -1;
   if ( *(float*)a > *(float*)b ) return 1;
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

int scale(long m,int *i,int *j,float *x, float *z,float *b, float *report,int *allIters, float tol,float *pppp,float *pppp1, int maxiter, int zerodiag, float del, float dp, float dp1,int *totIter, int threads)
{  
	long p;
	int n0;
	float low, zhigh, zlow;
	int *bad, lind, hind;
	float *row, *col, *r0, *s, *one, *zz, *row0, *dr, *dc;	

//	find the matrix dimensions
	int k = 0;
	for (p=0; p<m;p++) if (j[p] > k) k=j[p];
	k++;

//	now we can allocate the required arrays
	float *current = (float *) malloc(k*sizeof(float));
	row = (float *) malloc(k*sizeof(float));
	row0 = (float *) malloc(k*sizeof(float));
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
	float **space = (float **) malloc(threads*sizeof(float *));
	for (p=0; p<threads;p++) space[p] = (float *) malloc(k*sizeof(float));

//	for convenience and in order not to change the existing code
        float perc = *pppp;
        float perc1 = *pppp1;

//	treat possible NaNs in the scaling vector and remove perc1 proportion of highest and lowest positive entries
	for (p=0;p<k;p++) {
		if (isnan(z[p])) continue; 
		if (z[p] > 0) zz[l++] = z[p];
	}
	qsort(zz,l,sizeof(float),cmpfunc);
	lind = (int)(l*perc1+0.5);
	hind = (int)(l*(1.0-perc1)+0.5);
	if (lind < 0) lind = 0;
	if (hind >= l) hind = l-1;
	zlow = zz[lind];
	zhigh = zz[hind];
	for (p=0;p<k;p++) if (z[p] > 0 && (z[p] < zlow || z[p] > zhigh)) z[p] = NAN;

	for (p=0;p<k;p++) one[p] = 1.0;
        for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;

//	if zerodiag is 1 remove rows with 0 on the diagonal
	if (zerodiag == 1) {
		for (p=0;p<k;p++) bad[p] = 1;
		for (p=0;p<m;p++) if (i[p] == j[p]) bad[i[p]] = 0;
	} else {
		for (p=0;p<k;p++) bad[p] = 0;
	}

//	find rows sums
	utmvMul(i,j,x,m,one,k,row0,threads,space);
//	find relevant percentiles
	for (p=0; p<k;p++) row[p] = row0[p];
	n0 = 0;
        for (p=0; p<k;p++) if (row0[p] > 0) r0[n0++] = row0[p];
        qsort(r0,n0,sizeof(float),cmpfunc);
	lind =  (int)(n0*perc+0.5);
        if (lind < 0) lind = 0;
        low = r0[lind];
//	find rows which are identically 0 and the perc proportion of rows with the lowest row sums and exclude them
	for (p=0;p<k;p++) {
		if ((row[p] <= low  && z[p] > 0) || isnan(z[p])) {
			bad[p] = 1;
			z[p] = 1.0;
		}
	}
//	if bad[p] is 1 we remove row p


	for (p=0;p<k;p++) dr[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) dc[p] = 1.0 - bad[p];
	for (p=0;p<k;p++) one[p] = 1.0 - bad[p];
	
//	treat separately rows for which z[p] = 0
	for (p=0;p<k;p++) if (z[p] == 0) one[p] = 0;
	for (p=0;p<k;p++) bad1[p] = 1.0 - one[p];
	
//	start iterations
//	row is the current rows sum; dr and dc are the current rows and columns scaling vectors
	float ber = 10.0*(1.0+tol);
	float err = 10.0*(1.0+tol);
	int iter=0;
	for (p=0;p<k;p++) current[p] = sqrt(dr[p]*dc[p]);
      	int fail;   // checks whether the convergence rate seems to be good enough; 0 if yes, 1 if no
        int nerr = 0;
        float errors[10000];
	int all_iters = 0;

	while((ber > tol || err > 5.0*tol) && iter < maxiter) {
		iter++;
		all_iters++;
		fail = 1;
		for (p=0;p<k;p++) if (bad1[p] == 1) row[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/row[p];
		for (p=0;p<k;p++) dr[p] *= s[p];

//	find column sums and update rows scaling vector
		utmvMul(i,j,x,m,dr,k,col,threads,space);
		for (p=0;p<k;p++) col[p] *= dc[p];
		for (p=0;p<k;p++) if (bad1[p] == 1) col[p] = 1.0;
		for (p=0;p<k;p++) s[p] = z[p]/col[p];
		for (p=0;p<k;p++) dc[p] *= s[p];

//	find row sums and update columns scaling vector
		utmvMul(i,j,x,m,dc,k,row,threads,space);
		for (p=0;p<k;p++) row[p] *= dr[p];
//	calculate current scaling vector
		for (p=0;p<k;p++) b[p] = sqrt(dr[p]*dc[p]);

//	calculate the current error
		ber=0;
		for (p=0;p<k;p++) {
			if (bad1[p] == 1) continue;
			if (fabs((double) (b[p]-current[p])) > ber) ber = fabs((double) (b[p]-current[p]));
		}
		report[all_iters-1]=ber;
		allIters[all_iters-1]=iter;

//	since calculating the error in row sums requires matrix-vector multiplication we are are doing this every 10 
//	iterations
		if (iter % 10 == 0) {
			utmvMul(i,j,x,m,b,k,col,threads,space);
			err = 0;
			for (p=0;p<k;p++) {
				if (bad1[p] == 1) continue;
				if (err < fabs((double) (col[p]*b[p] - z[p]))) err = fabs((double) (col[p]*b[p] - z[p]));
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
//	the scaling vector does not seen to converge well enough or the row sums erros do not decrease enough - 
//	increase perc and perc1 and continue
                        if (fail == 1) {
                                perc += dp;
                                perc1 += dp1;
                                nerr = 0;
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
//	if the current error is larger than 5 iteration ago start from scratch, otherwise continue from the current 
//	position 
                                if (report[all_iters-1] > report[all_iters-6]) {
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

//	if perc or perc1 reached upper bound or the total number of iterationbs is too high, exit
		if (perc > 0.2 || perc1 > 0.1) break;
		if (all_iters > *totIter) break;
	}

//	find the final error in row sums
	if (iter % 10 == 0) {
		utmvMul(i,j,x,m,b,k,col,threads,space);
		err = 0;
		for (p=0;p<k;p++) {
			if (bad1[p] == 1) continue;
			if (err < fabs((double) (col[p]*b[p] - z[p]))) err = fabs((double) (col[p]*b[p] - z[p]));
		}
	}
	report[all_iters+1] = ber;
	report[all_iters+2] = err;
	*pppp = perc;
	*pppp1 = perc1;
	*totIter = all_iters;

	for (p=0;p<k;p++) if (bad[p] == 1) b[p] = NAN;
	
	return(iter);	
}

