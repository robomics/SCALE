#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int main(int argc, char *argv[]) {
	int version;
	string norm("NONE");
	string unit("BP");
	string ALL("ALL");
	string MT("MT");
	ifstream fin;
	struct timeb t0,t1;

	string fname = argv[1];
	fin.open(fname, fstream::in);
	int binsize = atoi(argv[2]);
	char *out_name = argv[3];
	FILE *fout = fopen(out_name,"w");

	string str;
	string one = "1";
	ftime(&t0);
	getline(fin, str, '\0' );
//	cout << str << "\n";
	fin.read((char*)&version, sizeof(int));
	if (version < 6) {
		cerr << "Version " << version << " no longer supported" << endl;
		 exit(1);
	}
//	cout << "version = " << version << "\n";
	long master;
	fin.read((char*)&master, sizeof(long));
	string genome;
	getline(fin, genome, '\0' );
	int nattributes;
	fin.read((char*)&nattributes, sizeof(int));
//	cout << master << " " << genome << " " << nattributes << "\n";

// reading and ignoring attribute-value dictionary
	for (int i=0; i<nattributes; i++) {
		string key, value;
		getline(fin, key, '\0');
		getline(fin, value, '\0');
//		cout << key << " " << value << "\n";
//		cout << key << "\n";
	}
	int nChrs;
	fin.read((char*)&nChrs, sizeof(int));
	vector<std::string> chroms;
	vector<int> chrLen;
// chromosome map for finding matrix
	for (int i=0; i<nChrs; i++) {
		string name;
		int length;
		getline(fin, name, '\0');
		fin.read((char*)&length, sizeof(int));
//		if (name != ALL && name != MT) {
		if (name != ALL && name != "All") {
				chroms.insert(chroms.end(),name);
				chrLen.insert(chrLen.end(),length);
			}
    	}
//	for (int i=0; i<chroms.size(); i++) cout << chroms.at(i) << " - " << chrLen.at(i) << "\n";
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=1; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	fin.close();
//	fin.open(fname, fstream::in);
//	fin.seekg(0, ios::beg);
	long nonZer = 0;
	ftime(&t0);
	int ii,jj;
	float xx;
	long pos=0;
	for (int j=1; j<chroms.size();j++) {
		vector<contactRecord> records = straw(norm, fname, one, chroms.at(j), unit, binsize);
		long length=records.size();
		for (long k=0; k<length; k++) {
			ii = records[k].binX/binsize; 
			jj = offset[j] + records[k].binY/binsize; 
			xx = (float) records[k].counts;
			pos++;
//			fprintf(fout,"%d %d %g\n",ii,jj,xx);	
		}
	}
	ftime(&t1);
	printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),pos);
	fclose(fout);
	return 0;
}
