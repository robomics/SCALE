#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
using namespace std;

int main(int argc, char *argv[]) {
	int version;
	string fname = argv[1];
	int binsize = atoi(argv[2]);
	string norm = "NONE";
	string unit = "BP";
	ifstream fin;
	fin.open(fname, fstream::in);
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
	vector<string> chroms;
	vector<int> chrLen;
	int length;
	string name;
	for (int i=0; i<nChrs; i++) {
		getline(fin, name, '\0');
		fin.read((char*)&length, sizeof(int));
		if (name != "ALL") {
			chroms.insert(chroms.end(),name);
			chrLen.insert(chrLen.end(),length);
		}
    	}
	int numChroms = chroms.size();
//	for (int i=0; i<numChroms; i++) cout << chroms.at(i) << " - " << chrLen.at(i) << "\n";
	fin.close();
	for (int i=0; i<numChroms; i++) {
		for (int j=i; j<numChroms;j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			vector<contactRecord> records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			int len = records.size();
			if (num != len) cout << chroms.at(i) << "-" << chroms.at(j) << ": " << num << " " << len << endl;
		}
		cout << "done " << chroms.at(i) << endl;
	}
	return 0;
}

