#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>

int scale(int c, int *m,int **i,int **j,float **x, float *z,float *b, float *report, int verb,float tol,float perc,float perc1, int maxiter,float del,int trials,int zerodiag);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-m memarray][-p percent][-q percent1][-v verbose][-s silent][-C max_array_count][-t tol][-d delta][-f trials][-a attempts] [-z remoce_zero_diagonal] <infile> <vector_file> <outfile>\n", argv0);
  fprintf(stderr, "  <infile>: matrix file in sparse upper triangular notation\n");
  fprintf(stderr, "  <vector_file>: vector to scale to, all 1s for balanced\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}


int main(int argc, char *argv[]) {
  int *m,k,p,q,c,n,ibad,iter;
  int **i, **j;
  float **x;
  float *z, *z0;
  struct timeb t0,t1,start,end;
  int opt;
  
  // parameters
  int maxC = 100;
  float tol=5.0e-4;
  int m1 = (int) 7e8; 
  float perc = 1.0e-2;
  float perc1=0.25e-2;
  int verb=0;
  int silent=1;
  int maxiter=300;
  float delta=0.01;
  int trials=5;
  int attempts=3;
  int zerodiag=1;
  
  while ((opt = getopt(argc, argv, "m:p:q:v:s:C:t:d:f:I:a:z:h")) != -1) {
    switch (opt) {
    case 'm':
      m1 = (int) atof(optarg);
      break;
    case 'p':
      perc = atof(optarg);
      break;
    case 'q':
      perc1 = atof(optarg);
      break;
    case 'C':
      maxC = atoi(optarg);
      break;
    case 't':
      tol = atof(optarg);
      break;
    case 'd':
      delta = atof(optarg);
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

  m = (int *) malloc(maxC*sizeof(int));
  
  c = 0;
  i = (int **) malloc(maxC*sizeof(int *));
  j = (int **) malloc(maxC*sizeof(int *));
  x = (float **) malloc(maxC*sizeof(float *));
  while(c < maxC) {
    ftime(&t0);
    i[c] = (int *) malloc(m1*sizeof(int));
    j[c] = (int *) malloc(m1*sizeof(int));
    x[c] = (float *) malloc(m1*sizeof(float));
    k = 0;
    while(fscanf(fin,"%d %d %f",&i[c][k],&j[c][k],&x[c][k]) == 3) {
      k++;
      if (k == m1) break;
    }
    if (k == 0) break;	
    for (p=0;p<k;p++) {
      i[c][p]--;
      j[c][p]--;
    }
    m[c++] = k;
    ftime(&t1);
    printf("took %ld seconds to read %d records; round %d\n",(long) (t1.time - t0.time),k,c);
    fflush(stdout);
    if (k < m1) break;
  }
  fclose(fin);
  printf("finished reading\n");
  printf("took %ld seconds to read %ld records\n",(long) (t1.time - start.time),((long) m1)*(c-1)+k);

  k = 0;
  for (int ic=0;ic<c;ic++) for (int p=0; p<m[ic];p++) if (j[ic][p] > k) k=j[ic][p];
  k++;

  for (int ic=0;ic<c;ic++) for (int p=0; p<m[ic];p++) if (j[ic][p] == i[ic][p]) x[ic][p] *= 0.5;

  z = (float *) malloc((2*k+1000)*sizeof(float));
  z0 = (float *) malloc((2*k+1000)*sizeof(float));

  n = 0;
  while(fscanf(finV,"%f",&z0[n]) == 1) n++;
  if (n < k) {
	for (int p=n;p<k;p++) z0[p] = NAN;
	n = k;
  }
  for (p=0;p<n;p++) z[p] = z0[p];
  float *b = (float *) malloc((2*k+1)*sizeof(float));
  for (p=0;p<n;p++) b[p] = NAN;
  fclose(finV);

  float *report = (float *) malloc((maxiter+1)*sizeof(float));
  
  ftime(&t1);
  iter = scale(c, m,i,j,x, z,b,report,1-silent,tol,perc,perc1,maxiter,delta,trials,zerodiag);

  ftime(&end);
  printf("took %ld seconds\n",(long) (end.time - start.time));
  printf("iterations took %15.10lf seconds\n",((float) (end.time - t1.time)) + 0.001*(end.millitm - t1.millitm) );
  if (silent) for (p=0;p<abs(iter);p++) printf("%d: %30.15lf\n",p+1,report[p]);
  int count=0;
  while (iter < 0 && count++ < attempts) {
    printf("Did not converge!!!\n");
    perc=1.5*perc;
    perc1=1.5*perc1;
    printf("new perc = %lg and new perc1 = %lg\n",perc,perc1);
    for (p=0;p<n;p++) z[p] = z0[p];
    ftime(&t1);
    //		iter = scale(c, m,i,j,x, z,b, tol,perc,perc1,maxiter, report,verb);
    iter = scale(c, m,i,j,x, z,b,report,1-silent,tol,perc,perc1,maxiter,delta,trials,zerodiag);
    ftime(&end);
    printf("iterations took %15.10lf seconds\n",((float) (end.time - t1.time)) + 0.001*(end.millitm - t1.millitm) );
    if (silent) for (p=0;p<abs(iter);p++) printf("%d: %30.15lf\n",p+1,report[p]);
  }
  
  for (p=0;p<n;p++) fprintf(fout,"%30.15lf\n",b[p]);
  fclose(fout);
  if (iter < 0) printf("still did not converge!!!\n");
  return(iter);
}

