#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
using namespace std;

int main(int argc, char *argv[]) {
  string fname = argv[1];
  int binsize = atoi(argv[2]);
  ifstream fin;
  long master;
  fin.open(fname, fstream::in);
  map<string, chromosome> chromosomeMap = readHeader(fin, master);

  
  for (const auto &chrEntry1 : chromosomeMap) {
    if (chrEntry1.first != "ALL" && chrEntry1.first != "All") {
      chromosome chr1 = chrEntry1.second;
      for (const auto &chrEntry2 : chromosomeMap) {
	if (chrEntry2.first != "ALL" && chrEntry2.first != "All") {
	  chromosome chr2 = chrEntry2.second;
	  if (chr1.index <= chr2.index) {
	    cerr << "calling getSize " <<chr1.name << " " << chr2.name << " " << endl;
	    int num = getSize("NONE", fname, chr1.name, chr2.name, "BP", binsize);
	    cerr << "num=" << num << endl;
	    vector<contactRecord> records = straw("NONE", fname, chr1.name, chr2.name, "BP", binsize);
	    int len = records.size();
	    if (num != len) {
	      cout << chr1.name << "-" << chr2.name << ": " << num << " " << len << "\n";
	    }	    
	    else cout << "good" << " " << chr1.name << "-" << chr2.name << endl;
	  } 
	}
      }
    }
  }
  return 0;
}

