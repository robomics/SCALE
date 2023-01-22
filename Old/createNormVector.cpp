#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int scale(long m,int *i,int *j,float *x, float *z,float *b, float *report, int *iters, float tol,float *pppp,float *pppp1, int maxiter, int zerodiag, float del, float dp, float dp1, int *totalIt, int threads, int k);

int getMatrix(string fname, int binsize, string norm, int **i, int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL);

double ppNormVector(long m,int *ii,int *jj,float *xx, float *b,int k, int threads, float **space);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-p percent][-P delta_p][-q percent1][-Q delta_q][-v verbose][-t tol][-I max_iterations][-z remove_zero_diag][-d delta][-A All_iterations][-T threads] <hicfile> <resolution> <outfile> <[vector_file]>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
  fprintf(stderr, "  [optional] <vecor_file>: vecor to scale to (otherwise balancing)\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	struct timeb t0,t1;

	float tol=5.0e-4;
	float del=2.0e-2;
	float perc = 1.0e-2;
	float dp = 5.0e-3;
	float perc1 = 0.25e-2;
	float dp1 = 1.0e-3;
	int maxiter=100;
	int zerodiag=0;
	int All_iter = 200;
	int threads = 1;
	int verb = 1;
	int opt;

	while ((opt = getopt(argc, argv, "m:p:P:q:Q:v:t:d:I:z:A:T:h")) != -1) {
		switch (opt) {
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

	FILE *fscal;
	bool scaling = false;
	if (argc - optind > 0) scaling = true;
	if (scaling) {
		char *scal_name = argv[optind++];
		fscal = fopen(scal_name,"r");
		if (fscal==NULL) {
			fprintf(stderr, "Error! File %s cannot be opened for reading\n", argv[optind-1]);
			exit(EXIT_FAILURE);
	  	}
	  	else if (verb) printf("Reading target vector from %s\n", argv[optind-1]);
	}
	printf("\n");

	int *ii;
	int *jj;
	float *xx;
	long p;
	long m;
	int k;
	vector<std::string> chroms;
	vector<int> chrLen;

	ftime(&t0);
	k = getMatrix(fname, binsize, norm, &ii, &jj, &xx, &m,chroms,chrLen);
	ftime(&t1);
	if (k < 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}
	if (verb) printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),m);

	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

	float *z = (float *) malloc(k*sizeof(float));
	string norm_name;
	if (scaling) {
		norm_name = "GW_SCAL";
		int n = 0;
		while(fscanf(fscal,"%f",&z[n]) == 1 && n < k) n++;
		if (n < k) for (p=n;p<k;p++) z[p] = NAN;
		fclose(fscal);
	}
	else {
		norm_name = "GW_BAL";
		for (p=0;p<k;p++) z[p] = 1.0;
		perc1 = 0;
		dp1 = 0;	
	}

	float *b = (float *) malloc(k*sizeof(float));
	for (p=0;p<k;p++) b[p] = NAN;
	int totalIt = All_iter;
	float *report = (float *) malloc((All_iter+3+100)*sizeof(float));
	int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));

	int iter;
	ftime(&t0);
	iter = scale(m,ii,jj,xx, z,b, report,allIters, tol,&perc,&perc1, maxiter, zerodiag, del, dp, dp1, &totalIt,threads,k);
	ftime(&t1);

	if (verb) printf("iterations took %15.10lf seconds\n",((float) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
	for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
	if (verb) printf("total %d iterations; final perc = %g and perc1 = %g\n",totalIt,perc,perc1);
	if (verb) printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);

	float **space = (float **) malloc(threads*sizeof(float *));
	for (p=0;p<threads; p++) space[p] = (float *) malloc(k*sizeof(float));
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
	fclose(fout);
	return(iter);
}

