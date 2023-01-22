#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <cmath>
#include <unistd.h>
using namespace std;

map<string, chromosome> readHeader(istream &fin, long &masterIndexPosition, string &genomeID, int &numChromosomes,
                                   int &version, long &nviPosition, long &nviLength);

static void usage(const char *argv0)
{
  fprintf(stderr, "Usage: %s [-n norm] <hicfile> <resolution> <outfile>\n", argv0);
  fprintf(stderr, "  <hicfile>: hic file\n");
  fprintf(stderr, "  <resolution>: resolution in bp\n");
  fprintf(stderr, "  <outfile>: output file\n");
}

int main(int argc, char *argv[]) {
        string norm("NONE");
        string ob("observed");
        string unit("BP");
	ifstream fin;
        int opt;

        while ((opt = getopt(argc, argv, "n:h")) != -1) {
                switch (opt) {
                                case 'n':
					norm=optarg;
                                        break;
                                case 'h':
                                        usage(argv[0]);
                                        exit(EXIT_SUCCESS);
                                default:
                                        usage(argv[0]);
                                        exit(EXIT_FAILURE);
                }
        }

        if (argc - optind < 3) {
                usage(argv[0]);
                exit(EXIT_FAILURE);
        }

	string fname = argv[optind++];
        int binsize = atoi(argv[optind++]);
        char *out_name = argv[optind++];

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
        if (!fin) {
                cerr << "File " << fname << " cannot be opened for reading" << endl;
                exit(EXIT_FAILURE);
        }

        FILE *fout = fopen(out_name,"w");
        if (fout==NULL) {
                fprintf(stderr, "Error! File %s cannot be opened for writing\n", argv[optind-1]);
                exit(EXIT_FAILURE);
          }

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

	long ii,jj;
	float xx;
	vector<contactRecord> records;
	for (int i=0; i<chroms.size(); i++) {
		for (int j=i; j<chroms.size();j++) {
			records = straw(ob,norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
			for (long k=0; k<length; k++) {
				ii = offset[i] + records[k].binX/binsize; 
				jj = offset[j] + records[k].binY/binsize; 
				xx = (float) records[k].counts;
				fprintf(fout,"%ld %ld %f\n",ii,jj,xx);
			}
			records.clear();
			records.shrink_to_fit();
		}
	}
	fin.close();
	fclose(fout);
	return 0;
}
