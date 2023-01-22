#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/timeb.h>
#include <time.h>

int scale(long m,int *i,int *j,float *x, float *z,float *b, float *report, int *iters, float tol,float *pppp,float *pppp1, int maxiter, int zerodiag, float del, float dp, float dp1, int *totalIt, int threads);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-m memarray][-p percent][-P delta_p][-q percent1]i[-Q delta_q][-v verbose][-t tol][-I max_iterations][-z remove_zero_diag][-d delta][-A All_iterations][-T threads] <infile> <vecor_file> <outfile>\n", argv0);
  fprintf(stderr, "  <infile>: matrix file in sparse upperc triangular notation\n");
  fprintf(stderr, "  <vecor_file>: vecor to scale to, all 1s for balaned\n");
}


int main(int argc, char *argv[]) {
  long m0,m1,m,k,p;
  int q,n,ibad,iter;
  int *i, *j;
  float *x;
  float *z, *z0;
  struct timeb t0,t1,start,end;
  int opt;
  
  // parameters
  float tol=5.0e-4;
  float del=2.0e-2;
  m1 = (long) 7e8; 
  float perc = 1.0e-2;
  float dp = 5.0e-3;
  float perc1=0.25e-2;
  float dp1 = 1.0e-3;
  int verb=1;
  int maxiter=100;
  int zerodiag=0;
  int All_iter = 200;
  int threads = 1;
  
  while ((opt = getopt(argc, argv, "m:p:P:q:Q:v:t:d:I:z:A:T:h")) != -1) {
    switch (opt) {
    case 'm':
      m1 = (long) atof(optarg);
      break;
    case 'p':
      perc = atof(optarg);
      break;
    case 'P':
      dp = atof(optarg);
      break;
    case 'q':
      perc1 = atof(optarg);
      break;
    case 'Q':
      dp1 = atof(optarg);
      break;
    case 't':
      tol = atof(optarg);
      break;
    case 'd':
      del = atof(optarg);
      break;
    case 'v':
      verb = atoi(optarg);
      break;
    case 'A':
      All_iter=atoi(optarg);
      break;
    case 'I':
      maxiter=atoi(optarg);
      break;
    case 'z':
      zerodiag=atoi(optarg);
      break;
    case 'T':
      threads=atoi(optarg);
      break;
    case 'h':
      usage(argv[0]);
      exit(EXIT_SUCCESS);
    default:
      usage(argv[0]);
      exit(EXIT_FAILURE);
    }
  }
  
  if (argc - optind < 2) {
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
  else if (verb) printf(" and %s\n", argv[optind-1]);
  fflush(stdout);
  
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
  x = (float *) malloc(m1*sizeof(float));
  if (x == NULL) {
    fprintf(stderr,"Failed to allocate %ld bytes. Exiting\n",m1*sizeof(float));
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
      x = (float *) realloc(x,m1*sizeof(float));
      if (x == NULL) {
        fprintf(stderr,"Failed to allocate additional %ld bytes. Exiting\n",m0*sizeof(float));
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

  z = (float *) malloc((2*k+1000)*sizeof(float));
  z0 = (float *) malloc((2*k+1000)*sizeof(float));

  n = 0;
  while(fscanf(finV,"%f",&z0[n]) == 1) n++;
  if (n < k) {
	for (p=n;p<k;p++) z[p] = NAN;
	n = k;
  }
  float *b = (float *) malloc((2*k+1)*sizeof(float));
  for (p=0;p<n;p++) b[p] = NAN;
  fclose(finV);

  int totalIt = All_iter;
  float *report = (float *) malloc((All_iter+3+100)*sizeof(float));
  int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));
  
  char answer;
  double p0 = perc, p10 = perc1;
  while (true) {
	for (p=0;p<n;p++) z[p] = z0[p];
	perc = p0;
	perc1 = p10;
	printf("Enter the number of threads\n");
	scanf("%d",&threads);
	ftime(&t0);
	iter = scale(m,i,j,x, z,b, report,allIters, tol,&perc,&perc1, maxiter, zerodiag, del, dp, dp1, &totalIt,threads);
	ftime(&t1);
  	printf("iterations took %15.10lf seconds\n",((float) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
  	printf("total %d iterations; final perc = %g and perc1 = %g\n",totalIt,perc,perc1);
  	printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);
	printf("Continue?\n");
	scanf("%s",&answer);
	if ((answer != 'y') && (answer != 'Y')) break;
  }

  return(iter);
}

