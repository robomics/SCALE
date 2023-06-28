#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>


#include "scale/scale.hpp"

using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);


unsigned int getMatrix(string fname, int binsize, string norm, unsigned int **i, unsigned int **j, float **x, long *m, string ob, string chr) {
	string norm1("NONE");
	string unit("BP");
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
    bool found = false;
    int length;
    for (map<string,chromosome>::iterator itr =  chromosomeMap.begin(); itr !=  chromosomeMap.end(); ++itr) {
         string chrom = itr->second.name;
	 if (chrom != chr) continue;
         length = itr->second.length;
	 found = true;
	 break;
    }
    if (!found) return 0;

	vector<contactRecord> records;
	records = straw(ob,norm1, fname, chr, chr, unit, binsize);
	long pos=0;
	unsigned int *ii = (unsigned int *) malloc(records.size()*sizeof(int)); 
	unsigned int *jj = (unsigned int *) malloc(records.size()*sizeof(int)); 
	float *xx = (float *) malloc(records.size()*sizeof(float)); 
	for (long k=0; k<records.size(); k++) {
		ii[pos] = records[k].binX/binsize; 
		jj[pos] = records[k].binY/binsize; 
		xx[pos] = (float) records[k].counts;
		if (!isnan(xx[pos])) pos++;
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	return (int) ceil(length/((double) binsize));
}

Interactions getSingleMatrix(const std::string& fname, int binsize, const std::string& chr) {


        long master = 0L;
        map<string, chromosome> chromosomeMap;
        string genomeID;
        int version = 0;
        long nviPosition = 0;
        long nviLength = 0;

        int nChrs;
        vector<std::string> chroms;
        vector<int> chrLen;
        std::ifstream fin(fname, fstream::in);
        chromosomeMap = readHeader(fin, master, genomeID, nChrs, version, nviPosition, nviLength);
        fin.close();
        bool found = false;
        int length;
        for (const auto& itr : chromosomeMap) {
                string chrom = itr.second.name;
                if (chrom != chr) continue;
                length = itr.second.length;
                found = true;
                break;
        }
        if (!found) {
                throw std::runtime_error("chromosome not found");
        }

        Interactions interactions{};
        vector<contactRecord> records;
        records = straw("observed", "NONE", fname, chr, chr, "BP", binsize);

        interactions.chrom1 = {chr, length};
        interactions.chrom2 = interactions.chrom1;
        interactions.ii.reserve(records.size());
        interactions.jj.reserve(records.size());
        interactions.xx.reserve(records.size());

        for (auto &record : records) {
                if (!isnan(record.counts)) {
                  interactions.ii.push_back(record.binX / binsize);
                  interactions.jj.push_back(record.binY / binsize);
                  interactions.xx.push_back(record.counts);
                }
        }

        return interactions;

}
