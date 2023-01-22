#include <cstring>
#include <iostream>
#include <fstream>
#include <straw.h>
#include <sys/timeb.h>
#include <time.h>
#include <cmath>
using namespace std;

int scale(long m,int *i,int *j,float *x, float *z,float *b, float *report, int *iters, float tol,float *pppp, int maxiter, float del, float dp, int *totalIt, int threads, int k);

int main(int argc, char *argv[]) {
	int version;
	string fname(argv[1]);
	int binsize = atoi(argv[2]);
	int threads = 1;
	if (argc >= 4) threads = atoi(argv[3]);
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
	cout << "k = " << current_offset << endl;
//	for (int i=0; i<chroms.size(); i++) cout << chroms.at(i) << " - " << chrLen.at(i) << " - " << offset[i] << "\n";
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
//			int num = getSize(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			vector<contactRecord> records = straw(norm, fname, chroms.at(i), chroms.at(j), unit, binsize);
			long length=records.size();
//			if (length != (long) num) cout << chroms.at(i) << " " << chroms.at(j) << " " << num << " " << length << "\n";		
			for (long k=0; k<length; k++) {
				ii[pos] = offset[i] + records[k].binX/binsize; 
				jj[pos] = offset[j] + records[k].binY/binsize; 
				xx[pos] = (float) records[k].counts;
				if (!isnan(xx[pos])) pos++;
			}
//			cout << chroms.at(i) << "-" << chroms.at(j) << ": " << pos << "\n";
		}
		cout << "done " << chroms.at(i) << endl;
	}
	ftime(&t1);
	printf("took %ld seconds for %ld records\n",(long) (t1.time - t0.time),pos);

	float tol=5.0e-4;
	float del=2.0e-2;
	float perc = 1.0e-2;
	float dp = 5.0e-3;
	int maxiter=100;
	int All_iter = 200;
	long p;
	long m = nonZer;
	int k = current_offset;
	for (p=0; p<m;p++) if (jj[p] == ii[p]) xx[p] *= 0.5;
	float *z = (float *) malloc(k*sizeof(float));
	for (p=0;p<k;p++) z[p] = 1.0;
	float *b = (float *) malloc(k*sizeof(float));
	for (p=0;p<k;p++) b[p] = NAN;
	int totalIt = All_iter;
	float *report = (float *) malloc((All_iter+3+100)*sizeof(float));
	int *allIters = (int *) malloc((All_iter+3+100)*sizeof(int));

	int iter;
	ftime(&t0);
	iter = scale(m,ii,jj,xx, z,b, report,allIters, tol,&perc, maxiter, del, dp, &totalIt,threads,current_offset);
	ftime(&t1);

	printf("iterations took %15.10lf seconds\n",((float) (t1.time - t0.time)) + 0.001*(t1.millitm - t0.millitm));
//	for (p=0;p<totalIt;p++) printf("%d: %30.15lf\n",allIters[p],report[p]);
	printf("total %d iterations; final perc = %g\n",totalIt,perc);
	printf("final error in scaling vector is %g and in row sums is %g\n",report[totalIt+1],report[totalIt+2]);
//	for (p=0;p<k;p++) fprintf(fout,"%30.15lf\n",b[p]);
//	fclose(fout);
	return(iter);
}

