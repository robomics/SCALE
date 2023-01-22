#include <cstring>
#include <iostream>
#include <fstream>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,double *pppp, int maxiter, double del, double dp, int *totIter, int threads, unsigned int k);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-m num. records][-p percent][-P delta_p][-v verbose][-t tol][-I max_iterations][-d delta][-A All_iterations][-T threads] <hicfile> <length> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <length>: chromosome length\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

int main(int argc, char *argv[]) {
	time_t t0,t1;

	long mm = 2.0e9;
	double tol=5.0e-4;
	double del=2.0e-2;
	double perc = 1.0e-2;
	double dp = 5.0e-3;
	int maxiter=100;
	int All_iter = 200;
	int threads = 1;
	int verb = 2;
	int opt;

	while ((opt = getopt(argc, argv, "m:p:P:v:t:d:I:A:T:h")) != -1) {
		switch (opt) {
				case 'm':
					mm = strtoll(optarg,NULL,10);
					break;
				case 'P':
					dp = atof(optarg);
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

	if (argc - optind < 4) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	char *in_name = argv[optind++];
	if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	int length = atoi(argv[optind++]);
	int binsize = atoi(argv[optind++]);
	char *out_name = argv[optind++];

	FILE *fin = fopen(in_name,"r");
	if (fin==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	}
	else if (verb) printf("Reading matrix from %s\n", argv[optind-1]);

	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	}
	else if (verb) printf("Writing norm vector to %s\n", argv[optind-1]);

	unsigned int *ii = (unsigned int *) malloc(mm*sizeof(int));
	unsigned int *jj = (unsigned int *) malloc(mm*sizeof(int));
	float *xx = (float *) malloc(mm*sizeof(float));
	long p;
	long m;
	unsigned int k = ceil(length/binsize);

	p = 0;
	m = 0;

	time(&t0);
	while(fscanf(fin,"%d %d %f",&ii[m],&jj[m],&xx[m]) == 3) m++;
	time(&t1);
	if (verb) printf("took %ld seconds for %ld records\n",t1 - t0,m);
	for (p=0;p<m;p++) {
        	ii[p]/=binsize;
        	jj[p]/=binsize;
	}

	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

	double *b = (double *) malloc(k*sizeof(double));
	for (p=0;p<k;p++) b[p] = NAN;
	int totalIt = All_iter;
	double *report = (double *) malloc((All_iter+3+100)*sizeof(double));
	int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));

	int iter;
	time(&t0);
	iter = balance(m,ii,jj,xx, b, report,allIters, tol,&perc, maxiter, del, dp, &totalIt, threads, k);
	time(&t1);

	if (verb) printf("iterations took %ld seconds\n",t1 - t0);
	if (verb > 1) for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
	if (verb) printf("total %d iterations; final perc = %g\n",totalIt,perc);
	if (verb) printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);

	for (int p=0;p<k;p++) {
		if (!std::isnan(b[p])) fprintf(fout,"%20.10f\n",b[p]);
		else fprintf(fout,"%20s\n","NaN");
	}
	fclose(fout);
	return(iter);
}

