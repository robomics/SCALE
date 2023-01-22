#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);

int getHiCInfo(string fname, vector<std::string> &CH,vector<int> &CHL) {

	 ifstream fin;
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

	CH = chroms;
	CHL = chrLen;
	return nChrs;
}
