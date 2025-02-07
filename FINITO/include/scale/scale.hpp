// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>


int balance(long m,unsigned int *i,unsigned int *j,float *x, double *b, double *report,int *allIters, double tol,int *llow, int maxiter, double del, int *totIter, int threads, unsigned int k, double *pppp, int width);

unsigned int getMatrix(std::string fname, int binsize, std::string norm, unsigned int **i, unsigned int **j, float **x, long *m, std::string ob, std::string chr);

double ppNormVector(long m,unsigned int *ii,unsigned int *jj,float *xx, double *b,unsigned int k, int threads, double **space);


// multithreaded matrix-vector multiplication
void utmvMul(unsigned int *i,unsigned int *j,float *x,long m,double *v,unsigned int k,double *res, int threads, double **space);
