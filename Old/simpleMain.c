#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>
#include <real.h>

int scale(long m,int *i,int *j,real *x, real *z,real *b, real *report, int verb,real tol,real perc,real perc1, int maxiter,real del,int trials,int zerodiag);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-m memarray][-p percent][-q percent1][-v verbose][-s silent][-t tol][-d delta][-f trials][-a attempts] [-z remove_zero_diag] <infile> <vecor_file> <outfile>\n", argv0);
  fprintf(stderr, "  <infile>: matrix file in sparse upperc triangular notation\n");
  fprintf(stderr, "  <vecor_file>: vecor to scale to, all 1s for balaned\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}


int main(int argc, char *argv[]) {
  long m0,m1,m,k,p;
  int q,n,ibad,iter;
  int *i, *j;
  real *x;
  real *z, *z0;
  struct timeb t0,t1,start,end;
  int opt;
  
  // parameters
  real tol=5.0e-4;
  m1 = (long) 7e8; 
  real perc = 1.0e-2;
  real perc1=0.25e-2;
  int verb=0;
  int silent=1;
  int maxiter=300;
  real delta=0.01;
  int trials=5;
  int attempts=3;
  int zerodiag=1;
  
  while ((opt = getopt(argc, argv, "m:p:q:v:s:t:d:f:I:a:z:h")) != -1) {
    switch (opt) {
    case 'm':
      m1 = (long) atof(optarg);
      break;
    case 'p':
      perc = (real) atof(optarg);
      break;
    case 'q':
      perc1 = (real) atof(optarg);
      break;
    case 't':
      tol = (real) atof(optarg);
      break;
    case 'd':
      delta = (real) atof(optarg);
      break;
    case 'f':
      trials = atoi(optarg);
      break;
    case 'v':
      verb = atoi(optarg);
      break;
    case 's':
      silent = atoi(optarg);
      break;
    case 'I':
      maxiter=atoi(optarg);
      break;
    case 'a':
      attempts=atoi(optarg);
      break;
    case 'z':
      zerodiag=atoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  if (argc - optind < 3) {
    usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  FILE *fin = fopen(argv[optind++],"r");
  if (fin==NULL) {
    fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
    exit(EXIT_FAILURE);
  }
  else if (verb) printf("Reading from %s", argv[optind-1]);
  
  FILE *finV = fopen(argv[optind++],"r");
  if (finV==NULL) {
    fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
    exit(EXIT_FAILURE);
  }
  else if (verb) printf(" and %s", argv[optind-1]);
  

  FILE *fout = fopen(argv[optind++],"w");
  if (fout==NULL) {
    fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
    exit(EXIT_FAILURE);
  }
  else if (verb) printf(" and writing to %s\n", argv[optind-1]);
  
  ftime(&start);

  i = (int *) malloc(m1*sizeof(int));
  if (i == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(int));
    exit(EXIT_FAILURE);
  }
  j = (int *) malloc(m1*sizeof(int));
  if (j == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(int));
    exit(EXIT_FAILURE);
  }
  x = (real *) malloc(m1*sizeof(real));
  if (x == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(real));
    exit(EXIT_FAILURE);
  }
  m = 0; 
  m0 = m1;
  ftime(&t0);
  while(fscanf(fin,"%d %d %f",&i[m],&j[m],&x[m]) == 3) {
    m++;
    if (m == m1) {
      m1 += m0;
      i = (int *) realloc(i,m1*sizeof(int));
      if (i == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(int));
        exit(EXIT_FAILURE);
      }
      j = (int *) realloc(j,m1*sizeof(int));
      if (j == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(int));
        exit(EXIT_FAILURE);
      }
      x = (real *) realloc(x,m1*sizeof(real));
      if (x == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(real));
        exit(EXIT_FAILURE);
      }
    }
  }
   
  for (p=0;p<m;p++) {
    i[p]--;
    j[p]--;
  }
  ftime(&t1);
  fclose(fin);
  printf("finished reading\n");
  printf("took %ld seconds to read %ld records\n",(long) (t1.time - start.time),m);

  k = 0;
  for (p=0; p<m;p++) if (j[p] > k) k=j[p];
  k++;

  for (p=0; p<m;p++) if (j[p] == i[p]) x[p] *= 0.5;

  z = (real *) malloc((2*k+1000)*sizeof(real));
  z0 = (real *) malloc((2*k+1000)*sizeof(real));

  n = 0;
  while(fscanf(finV,"%f",&z0[n]) == 1) n++;
  if (n < k) {
	for (p=n;p<k;p++) z0[p] = NAN;
	n = k;
  }
  for (p=0;p<n;p++) z[p] = z0[p];
  real *b = (real *) malloc((2*k+1)*sizeof(real));
  for (p=0;p<n;p++) b[p] = NAN;
  fclose(finV);

  real *report = (real *) malloc((maxiter+3)*sizeof(real));
  
  ftime(&t1);
  iter = scale(m,i,j,x, z,b,report,1-silent,tol,perc,perc1,maxiter,delta,trials,zerodiag);

  ftime(&end);
  printf("took %ld seconds\n",(long) (end.time - start.time));
  printf("iterations took %15.10lf seconds\n",((real) (end.time - t1.time)) + 0.001*(end.millitm - t1.millitm) );
  if (silent) for (p=0;p<abs(iter);p++) printf("%d: %30.15lf\n",(int) (p+1),report[p]);
  printf("final error in scaling vector is %g and in row sums is %g\n",report[maxiter+1],report[maxiter+2]);
  int count=0;
  while (iter < 0 && count++ < attempts) {
    printf("Did not converge!!!\n");
    perc=1.5*perc;
    perc1=1.5*perc1;
    printf("new perc = %lg and new perc1 = %lg\n",perc,perc1);
    for (p=0;p<n;p++) z[p] = z0[p];
    ftime(&t1);
    iter = scale(m,i,j,x, z,b,report,1-silent,tol,perc,perc1,maxiter,delta,trials,zerodiag);
    ftime(&end);
    printf("iterations took %15.10lf seconds\n",((real) (end.time - t1.time)) + 0.001*(end.millitm - t1.millitm) );
    if (silent) for (p=0;p<abs(iter);p++) printf("%d: %30.15lf\n",(int) (p+1),report[p]);
    printf("final error in scaling vector is %g and in row sums is %g\n",report[maxiter+1],report[maxiter+2]);
  }
  
  for (p=0;p<n;p++) fprintf(fout,"%30.15lf\n",b[p]);
  fclose(fout);
  if (iter < 0) printf("still did not converge!!!\n");
  return(iter);
}

