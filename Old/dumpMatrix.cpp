#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

int main(int argc, char *argv[]) {
	int version;
	string norm("NONE");
	string unit("BP");
	ifstream fin;

	fin.open(argv[1], fstream::in);
	int binsize = atoi(argv[2]);

	string str;
	getline(fin, str, '\0' );
	fin.read((char*)&version, sizeof(int));
	if (version < 6) {
		cerr << "Version " << version << " no longer supported" << endl;
		 exit(1);
	}
	long master;
	fin.read((char*)&master, sizeof(long));
	string genome;
	getline(fin, genome, '\0' );
	int nattributes;
	fin.read((char*)&nattributes, sizeof(int));

// reading and ignoring attribute-value dictionary
	for (int i=0; i<nattributes; i++) {
		string key, value;
		getline(fin, key, '\0');
		getline(fin, value, '\0');
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
		if (name != "ALL" && name != "All") {
				chroms.insert(chroms.end(),name);
				chrLen.insert(chrLen.end(),length);
			}
    	}
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	fin.close();
	int ii,jj;
	float xx;
	vector<contactRecord> records;
	FILE *fout = fopen(argv[3],"w");
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			records = straw(norm, argv[1], chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii = offset[i] + records[k].binX/binsize; 
				jj = offset[j] + records[k].binY/binsize; 
				xx = (float) records[k].counts;
				if (!isnan(xx)) fprintf(fout,"%d %d %g\n",ii,jj,xx);
			}
		        records.clear();
       			records.shrink_to_fit();
		}
	}
	fin.close();
	return(0);
}
