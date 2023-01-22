#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

long realSize(string fname, int binsize, bool intra);

int main(int argc, char *argv[]) {
	string norm("NONE");
	struct timeb t0,t1;

	string fname = argv[1];
	int binsize = atoi(argv[2]);
	bool intra = (strcmp(argv[3],"true")==1);

	long m;

	ftime(&t0);
	m = realSize(fname, binsize,intra);
	ftime(&t1);
	printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),m);

	return(0);
}

