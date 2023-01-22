#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;

int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,double *pppp, int maxiter, double del, double dp, int *totIter, int threads, unsigned int k);

//int getMatrix(string fname, int binsize, string norm, int **i, int **j, double **x, long *m, vector<std::string> &CH,vector<int> &CHL);
unsigned int getMatrix(string fname, int binsize, string norm, string ob, bool interOnly, unsigned int **i, unsigned int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL);

double ppNormVector(long m,unsigned int *ii,unsigned int *jj,float *xx, double *b,unsigned int k, int threads, double **space);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-i (inter_only)][-p percent][-P delta_p][-v verbose][-t tol][-I max_iterations][-d delta][-A All_iterations][-T threads] <hicfile> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	string ob("observed");
	time_t t0,t1;
	string norm_name("GW_BAL");

	double tol=5.0e-4;
	double del=2.0e-2;
	double perc = 1.0e-2;
	double dp = 5.0e-3;
	int maxiter=100;
	int All_iter = 200;
	int threads = 1;
	int verb = 1;
	bool interOnly = false;
	int opt;

	while ((opt = getopt(argc, argv, "ip:P:v:t:d:I:A:T:h")) != -1) {
		switch (opt) {
				case 'i':
					interOnly = true;
					break;
				case 'p':
					perc = atof(optarg);
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

	if (argc - optind < 3) {
		usage(argv[0]);
		exit(EXIT_FAILURE);
	}

	string fname = argv[optind++];
	if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	int binsize = atoi(argv[optind++]);
	char *out_name = argv[optind++];

	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else if (verb) printf("Writing norm vector to %s\n", argv[optind-1]);

	unsigned int *ii;
	unsigned int *jj;
	float *xx;
	long p;
	long m;
	unsigned int k;
	vector<std::string> chroms;
	vector<int> chrLen;

	time(&t0);
	k = getMatrix(fname, binsize, norm, ob, interOnly, &ii, &jj, &xx, &m,chroms,chrLen);
	time(&t1);
	if (k == 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}
	if (verb) printf("took %ld seconds for %ld records\n",t1 - t0,m);

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

//	if (verb) printf("iterations took %15.10lf seconds\n",((double) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
	if (verb) printf("iterations took %ld seconds\n",t1 - t0);
	for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
	if (verb) printf("total %d iterations; final perc = %g\n",totalIt,perc);
	if (verb) printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);

	for (int j=0;j< k; j++) if (!isnan(b[j])) fprintf(fout,"%20.10f\n",b[j]);
		else fprintf(fout,"%20s\n","NaN");

/*
	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0;p<threads; p++) space[p] = (double *) malloc(k*sizeof(double));
	double sum = ppNormVector(m,ii,jj,xx,b,k,threads,space);

	int n = 0;
	for (int i=0; i<chroms.size(); i++) {
		fprintf(fout,"vector %s %s %d BP\n",const_cast<char*> (norm_name.c_str()),const_cast<char*> (chroms.at(i).c_str()),binsize);
		for (int j=0;j< ((int) ceil(chrLen.at(i)/((double) binsize)));j++) {
			if (!isnan(b[n])) fprintf(fout,"%20.10f\n",b[n++]);
			else {
				fprintf(fout,"%20s\n","NaN");
				n++;
			}
		}
	}
*/
	fclose(fout);
	return(iter);
}

