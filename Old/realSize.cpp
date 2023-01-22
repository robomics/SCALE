#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);


long realSize(string fname, int binsize, bool intra) {
	string norm("NONE");
	string unit("BP");
	ifstream fin;
        string ob("observed");


    long master = 0L;
    map<string, chromosome> chromosomeMap;
    string genomeID;
    int version = 0;
    long nviPosition = 0;
    long nviLength = 0;
    long totalFileSize;

	int nChrs;
	vector<std::string> chroms;
	vector<int> chrLen;
	fin.open(fname, fstream::in);
    chromosomeMap = readHeader(fin, master, genomeID, nChrs, version, nviPosition, nviLength);
    fin.close();
    for (map<string,chromosome>::iterator itr =  chromosomeMap.begin(); itr !=  chromosomeMap.end(); ++itr) {
         string chrom = itr->second.name;
         int length = itr->second.length;
	 if (chrom != "ALL" && chrom != "All") {
		chroms.insert(chroms.end(),chrom);
		chrLen.insert(chrLen.end(),length);
	}
    }

	long pos = 0;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size(); i++) {
		int upto = chroms.size();
		if (intra) upto = i+1;
		for (int j=i; j<chroms.size();j++) {
			records = straw(ob,norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			pos += length;
			records.clear();
			records.shrink_to_fit();
		}
	}
	fin.close();

	return pos;
}
