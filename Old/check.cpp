#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <ctime>
#include <cmath>
#include <unistd.h>
using namespace std;

unsigned int getMatrix(string fname, int binsize, string norm, string ob, bool interOnly, unsigned int **i, unsigned int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL);

int main(int argc, char *argv[]) {
	string norm("NONE");
	string ob("observed");

	bool interOnly = false;
	int opt;

	while ((opt = getopt(argc, argv, "i")) != -1) {
		switch (opt) {
				case 'i':
					interOnly = true;
					break;
				default:
					exit(EXIT_FAILURE);
		}
	}

	string fname = argv[1];
	int binsize = atoi(argv[2]);

	unsigned int *ii;
	unsigned int *jj;
	float *xx;
	long p;
	long m;
	unsigned int k;
	vector<std::string> chroms;
	vector<int> chrLen;

	k = getMatrix(fname, binsize, norm, ob, interOnly, &ii, &jj, &xx, &m,chroms,chrLen);
	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;

	double *b = (double *) malloc(k*sizeof(double));
	double *y = (double *) malloc(k*sizeof(double));
	int i;

	for (int f=3;f<argc;f++) {
		FILE *fvec = fopen(argv[f],"r");
		for (i=0;i<k;i++) fscanf(fvec,"%lf",&b[i]);
		for (i=0;i<k;i++) if (std::isnan(b[i])) b[i] = 0;
		for (i=0;i<k;i++) y[i] = 0;
		for (p=0; p<m;p++) {
			y[ii[p]] += xx[p]*b[jj[p]];
			y[jj[p]] += xx[p]*b[ii[p]];
		}
		for (i=0;i<k;i++) y[i] *= b[i];
		double l,u;
		i=0;
		while (b[i]==0) i++;
		l = y[i];
		u = y[i];
		for (;i<k;i++) {
			if (b[i] == 0) continue;
			if (y[i] < l) l = y[i];
			if (y[i] > u) u = y[i];
		}
		printf("%s %lg %lg %lg\n",argv[f],l,u,u-l);
		fclose(fvec);
	}
	return 0;
}

