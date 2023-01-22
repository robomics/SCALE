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

int getHiCInfo(string fname, vector<std::string> &CH,vector<int> &CHL);

double ppNormVector(long m,int *ii,int *jj,float *xx, float *b,int k, int threads, float **space);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-o oe] [-p percent][-P delta_p][-v verbose][-t tol][-I max_iterations][-z remove_zero_diag][-d delta][-A All_iterations][-T threads] <hicfile> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	string unit("BP");
	string matrix("observed");
	string norm_name("INTRA_BAL");
	struct timeb t0,t1,start,end;

	float tol=5.0e-4;
	float del=2.0e-2;
	float perc0 = 1.0e-2;
	float dp = 5.0e-3;
	float perc1 = 0.25e-2;
	float dp1 = 1.0e-3;
	int maxiter=100;
	int zerodiag=0;
	int All_iter = 200;
	int threads = 1;
	int verb = 2;
	int opt;

	ftime(&start);

	while ((opt = getopt(argc, argv, "o:m:p:P:v:t:d:I:z:A:T:h")) != -1) {
		switch (opt) {
				case 'o':
					if (strcmp(optarg,"oe") == 0) matrix = "oe";
					else {
						usage(argv[0]);
						exit(EXIT_FAILURE);
					}
					break;
				case 'p':
					perc0 = atof(optarg);
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
	if (verb) printf("\nReading hic file from %s\n", argv[optind-1]);
	int binsize = atoi(argv[optind++]);
	char *out_name = argv[optind++];

	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	  else if (verb) printf("Writing norm vector to %s\n\n", argv[optind-1]);

	int *ii;
	int *jj;
	float *xx;
	long p;
	long m;
	int k;
	int i;
	vector<std::string> chroms;
	vector<int> chrLen;

	i = getHiCInfo(fname,chroms,chrLen);

	if (i < 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}

	vector<contactRecord> records;
	int iter0 = 0;
	float **space = (float **) malloc(threads*sizeof(float *));

//	loop accross all chromosomes
        for (i=0; i<chroms.size(); i++) {
		ftime(&t0);
        	records = straw(matrix,norm, fname, chroms.at(i), chroms.at(i), unit, binsize);
        	m = records.size();
		ii = (int *) malloc(m*sizeof(int));
		jj = (int *) malloc(m*sizeof(int));
		xx = (float *) malloc(m*sizeof(float));
                for (p=0; p<m; ) {
                	ii[p] = records[p].binX/binsize;
                        jj[p] = records[p].binY/binsize;
                        xx[p] = (float) records[p].counts;
                        if (!isnan(xx[p])) p++;
                }
        	records.clear();
        	records.shrink_to_fit();
		ftime(&t1);
		if (verb > 1) printf("%s: took %ld seconds for %ld records\n",const_cast<char*> (chroms.at(i).c_str()),(long) (t1.time - t0.time),m);
		m = p;
		for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

		k = (int) ceil(chrLen.at(i)/((double) binsize));
		float *z = (float *) malloc(k*sizeof(float));
		for (p=0;p<k;p++) z[p] = 1.0;
		perc1 = 0;
		dp1 = 0;	

		float *b = (float *) malloc(k*sizeof(float));
		for (p=0;p<k;p++) b[p] = NAN;
		int totalIt = All_iter;
		float *report = (float *) malloc((All_iter+3+100)*sizeof(float));
		int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));

		float perc=perc0;
		int iter;
		ftime(&t0);
		iter = scale(m,ii,jj,xx, z,b, report,allIters, tol,&perc,&perc1, maxiter, zerodiag, del, dp, dp1, &totalIt,threads,k);
		ftime(&t1);
		iter0 += iter;
		if (verb > 1) printf("iterations took %15.10lf seconds\n",((float) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
		if (verb > 2) for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
		if (verb > 1) {
			printf("total %d iterations; final perc = %g and perc1 = %g\n",totalIt,perc,perc1);
			printf("final error in scaling vector is %g and in row sums is %g\n\n",report[totalIt+1],report[totalIt+2]);
		}
		for (p=0;p<threads; p++) space[p] = (float *) malloc(k*sizeof(float));
		if (report[totalIt+1] >= tol || report[totalIt+2] >= 5.0*tol) for (p=0;p<k;p++) b[p] = NAN;
		else double sum = ppNormVector(m,ii,jj,xx,b,k,threads,space);
		fprintf(fout,"vector %s %s %d BP\n",const_cast<char*> (norm_name.c_str()),const_cast<char*> (chroms.at(i).c_str()),binsize);
		for (int j=0;j<k;j++) {
			if (!isnan(b[j])) fprintf(fout,"%20.10f\n",b[j]);
			else fprintf(fout,"%20s\n","NaN");
		}
		free(ii);
		free(jj);
		free(xx);
		for (p=0;p<threads; p++) free(space[p]);
	}
	fclose(fout);
	ftime(&end);
	if (verb) printf("\n**************    all togeter took %15.10lf seconds\n",((float) (end.time - start.time)) + 0.001*(end.millitm - start.millitm));
	return(iter0);
}

