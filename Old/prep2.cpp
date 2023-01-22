#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
using namespace std;

int main(int argc, char *argv[]) {
	int version;
//	string fname("LCL_mega_62B_1k_30.hic");
	string fname(argv[1]);
	int binsize = atoi(argv[2]);
	string norm("NONE");
	string unit("BP");
	string bad("ALL");
	ifstream fin;
	struct timeb t0,t1;
	fin.open(fname, fstream::in);
	string str;
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
	cout << master << " " << genome << " " << nattributes << "\n";

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
		if (name != bad) {
				chroms.insert(chroms.end(),name);
				chrLen.insert(chrLen.end(),length);
			}
    	}
//	for (int i=0; i<chroms.size(); i++) cout << chroms.at(i) << " - " << chrLen.at(i) << "\n";
	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	for (int i=0; i<chroms.size(); i++) cout << chroms.at(i) << " - " << chrLen.at(i) << " - " << offset[i] << "\n";
	fin.close();
//	fin.open(fname, fstream::in);
//	fin.seekg(0, ios::beg);
	long nonZer = 0;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			nonZer += num;
//			cout << chroms[i] << "-" << chroms[j] << ": " << num << "\n";
		}
	}
	ftime(&t1);
	printf("Total %ld records\n",nonZer);
	printf("took %ld seconds\n",(long) (t1.time - t0.time));
	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	float *xx = (float *) malloc(nonZer*sizeof(float));
	long pos = 0;
	ftime(&t0);
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			vector<contactRecord> records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			if (length != (long) num) cout << chroms.at(i) << " " << chroms.at(j) << " " << num << " " << length << "\n";		
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (float) records[k].counts;
				pos++;
			}
//			cout << chroms.at(i) << "-" << chroms.at(j) << ": " << pos << "\n";
		}
		cout << "done " << chroms.at(i) << endl;
	}
	ftime(&t1);
	printf("took %ld seconds\n",(long) (t1.time - t0.time));
	cout << "Enter anything" << endl;
	string answer;
	cin >> answer;
	cout << answer << endl;
	return 0;
}

