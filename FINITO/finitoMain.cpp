#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>

#include "scale/scale.hpp"

using namespace std;

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-e (oe)][-D (diag)][-p perc][-v verbose][-t tol][-I max_iterations][-A All_iterations][-T threads] <hicfile> <chromosome> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chr>: chromosome\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	time_t t0,t1;

	string ob("observed");
	string norm_name("SCALE");
	double tol=1.0e-4;
	double del=5.0e-2;
	double perc = 0.01;
	int low;
	int maxiter=1000;
	int All_iter = 2000;
	int threads = 1;
	int verb = 1;
	bool diag = false;
	int opt;

	while ((opt = getopt(argc, argv, "ep:v:t:Dd:I:A:T:h")) != -1) {
		switch (opt) {
				case 'e':
					ob = "oe";
                                        break;
				case 'D':
					diag = true;
                                        break;
				case 'p':
					perc = atof(optarg);
					break;
				case 'd':
					del = atof(optarg);
					break;
				case 't':
					tol = atof(optarg);
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

	string fname = argv[optind++];
	if (verb) printf("Reading hic file from %s\n", argv[optind-1]);
	string chr = argv[optind++];
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
	int width = -1;
	if (diag) {
		width = (int) (10000/binsize);
		if (width < 2) width = 2;
	}

	time(&t0);
	k = getMatrix(fname, binsize, norm, &ii, &jj, &xx, &m, ob, chr);
	time(&t1);
	if (k < 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}
	else if (k == 0) {
		cerr << "Error! chromosome " << chr << " is not found in " << fname << endl;
		exit(EXIT_FAILURE);
	}
	if (verb) printf("took %ld seconds for %ld records\n",t1 - t0,m);

	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;


	std::vector<double> b(k);
	for (p=0;p<k;p++) b[p] = NAN;
	int totalIt = All_iter;
	double *report = (double *) malloc((All_iter+3+100)*sizeof(double));
	int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));

	int iter;
	time(&t0);
	iter =  balance(m,ii,jj,xx, b.data(), report,allIters, tol,&low, maxiter, del, &totalIt, threads, k, &perc,width);
	time(&t1);



  for (const auto&w: b) {
    printf("%f\n", w);
  }

	if (verb) printf("iterations took %ld seconds\n",t1 - t0);
	if (verb > 1) for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
	if (verb) printf("total %d iterations; final perc = %g; final cutoff = %d\n",totalIt,perc,low);
	if (verb) printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);

	double **space = (double **) malloc(threads*sizeof(double *));
	for (p=0;p<threads; p++) space[p] = (double *) malloc(k*sizeof(double));
	double sum = ppNormVector(m,ii,jj,xx,b.data(),k,threads,space);

	fprintf(fout,"vector %s %s %d BP\n",const_cast<char*> (norm_name.c_str()),const_cast<char*> (chr.c_str()),binsize);
	for (int p=0;p<k;p++) {
		if (!std::isnan(b[p])) fprintf(fout,"%20.10f\n",b[p]);
		else fprintf(fout,"%20s\n","NaN");
	}
	fclose(fout);
	return(iter);
}

