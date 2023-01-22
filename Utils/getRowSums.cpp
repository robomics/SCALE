#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;

unsigned int getMatrix(string fname, int binsize, string norm, unsigned int **i, unsigned int **j, float **x, long *m, string ob, string chr);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s <hicfile> <chromosome> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <chr>: chromosome\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: normalization vector output  file\n");
}

int main(int argc, char *argv[]) {
	string norm("NONE");
	time_t t0,t1;

	string ob("observed");
	int opt;

	while ((opt = getopt(argc, argv, "ep:v:t:d:I:A:T:h")) != -1) {
		switch (opt) {
				case 'e':
					ob = "oe";
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
	printf("Reading hic file from %s\n", argv[optind-1]);
	string chr = argv[optind++];
	int binsize = atoi(argv[optind++]);
	char *out_name = argv[optind++];

	FILE *fout = fopen(out_name,"w");
	if (fout==NULL) {
		fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
		exit(EXIT_FAILURE);
	  }
	else printf("Writing norm vector to %s\n", argv[optind-1]);

	unsigned int *i;
	unsigned int *j;
	float *x;
	long p;
	long m;
	unsigned int k;

	time(&t0);
	k = getMatrix(fname, binsize, norm, &i, &j, &x, &m, ob, chr);
	time(&t1);
	if (k < 0) {
		cerr << "Error! File " << fname << " cannot be opened for reading" << endl;
		exit(EXIT_FAILURE);
	}
	else if (k == 0) {
		cerr << "Error! chromosome " << chr << " is not found in " << fname << endl;
		exit(EXIT_FAILURE);
	}
	printf("took %ld seconds for %ld records\n",t1 - t0,m);

	double *b = (double *) malloc(k*sizeof(double));
	int *nz = (int *) malloc(k*sizeof(int));
	for (p=0;p<k;p++) {
		b[p] = 0;
		nz[p] = 0;
	}

	for (p=0; p<m;p++) {
		b[i[p]] += x[p];
		nz[i[p]]++;
		if (j[p] != i[p]) {
			b[j[p]] += x[p];
			nz[j[p]]++;
		}
	}

	for (int p=0;p<k;p++) fprintf(fout,"%g %d\n",b[p],nz[p]);
	fclose(fout);
	return(0);
}

