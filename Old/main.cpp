/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <iostream>
#include <string>
#include "straw.h"
using namespace std;

int main(int argc, char *argv[])
{
    string matrixType = "observed";
    string norm = "NONE";
    string unit = "BP";
    string fname = argv[1];
    string chr1loc = argv[2];
    string chr2loc = argv[3];
    string size = argv[4];
    int binsize = stoi(size);
    vector<contactRecord> records;
    records = straw(matrixType, norm, fname, chr1loc, chr2loc, unit, binsize);
    size_t length = records.size();
    printf("NR = %ld\n",length);
/*
 * for (int i = 0; i < length; i++) {
        printf("%d\t%d\t%.14g\n", records[i].binX, records[i].binY, records[i].counts);
    }
*/
}
