#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);


int getMatrix(string fname, int binsize, string norm, int **i, int **j, float **x, long *m, vector<std::string> &CH,vector<int> &CHL) {
	string norm1("NONE");
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

	int offset[chroms.size()];
	int current_offset = 0;
	for (int i=0; i<chroms.size(); i++) {
		offset[i] = current_offset;
		current_offset += (int) ceil(chrLen.at(i)/((double) binsize));
	}
	long nonZer = 21e9;
/* getSize!!!
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			nonZer += num;
		}
	}
*/

	int *ii = (int *) malloc(nonZer*sizeof(int));
	int *jj = (int *) malloc(nonZer*sizeof(int));
	float *xx = (float *) malloc(nonZer*sizeof(float));
	long pos = 0;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size()-1; i++) {
		for (int j=i+1; j<chroms.size();j++) {
			records = straw(ob,norm1, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (float) records[k].counts;
				if (!isnan(xx[pos])) pos++;
			}
		}
	}
	fin.close();
	records.clear();
	records.shrink_to_fit();

	*i = ii;
	*j = jj;
	*x = xx;
	*m = pos;
	CH = chroms;
	CHL = chrLen;
	return current_offset;
}
