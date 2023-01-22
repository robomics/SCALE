#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int main(int argc, char *argv[]) {
	string fname = argv[1];
	int binsize = atoi(argv[2]);
	bool interOnly = false;
	if (atoi(argv[3]) == 1) interOnly = true;
	long nonZer = getNumRecordsForFile(fname,binsize,interOnly);
	cout << nonZer << " nonzeros" << endl;
	return 0;
}
