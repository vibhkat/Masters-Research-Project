#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include<complex>
#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<cctype>
#include"stdio.h"
#include"stdlib.h"
#include<fstream>
#include<cmath>
#include<cstring>
#include<ctime>
using namespace std;
int AinverseB(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
int QR_factorization(vector<vector<double> > &,vector<vector<double> > &);
int vectormultiplication(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
int vector_transpose(vector<vector<double> > &,vector<vector<double> > &);
int blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int & ,vector<double> &);
int multi_blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int &,vector<double> &);
int block_inverse(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,double );
int complex_vector_multiplication(vector<vector<complex<double> > > ,vector<vector<complex<double> > > ,vector<vector<complex<double> > > &);
int multi_Arnoldi(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,int, double );
int original_plus_MOR_giving_X_big(vector<vector<double> > &,double,double,double,string);
int SVD(vector<vector<double> > &,vector<double> &,vector<vector<double> >);
int retrun_X_MATRIX(vector<vector<double> > & ,vector<vector<double> >,int );
int original_plus_MOR_sending_X_MATRIX(vector<vector<double> > ,double ,double ,double );

extern "C" int zgesv_( int *n , int *nrhs , complex<double> *a , int *lda , int *ipiv ,complex<double> *b , int *ldb , int *info  );
extern "C" int dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
extern  "C" int dgeev_( char* jobvl, char* jobvr, int* n, double* a,int* lda,  double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info );
extern "C" void dgesvd_(char* jobu,char* jobvt,int* m,int* n,double* a,int* lda,double* s,double* u,int* ldu,double* vt,int* ldvt,double* work,int* lwork,int* info );
int file_write(string , vector<vector<double> > );

#define pi 3.14
#endif
