#include <cstring>
#include <iostream>
#include <fstream>
using namespace std;

int main(int argc, char *argv[]) {
	int version;
//	string fname("LCL_mega_62B_1k_30.hic");
	string fname(argv[1]);
	ifstream fin;
	fin.open(fname, fstream::in);
	string str;
	getline(fin, str, '\0' );
	cout << str << "\n";
	fin.read((char*)&version, sizeof(int));
	if (version < 6) {
		cerr << "Version " << version << " no longer supported" << endl;
		 exit(1);
	}
	cout << "version = " << version << "\n";
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
		cout << key << "\n";
	}
	int nChrs;
	fin.read((char*)&nChrs, sizeof(int));
// chromosome map for finding matrix
	for (int i=0; i<nChrs; i++) {
		string name;
		int length;
		getline(fin, name, '\0');
		fin.read((char*)&length, sizeof(int));
		cout << name << " - " << length << "\n";
    	}
//	fin.close();
//	fin.open(fname, fstream::in);
	fin.seekg(0, ios::beg);
	getline(fin, str, '\0' );
	cout << str << "\n";
	return 0;
}

