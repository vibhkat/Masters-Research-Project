#ifndef PRINT_H
#define PRINT_H
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
using namespace std;
int print(vector<vector<double> >&);
int printcomplexvector(vector<complex<double> >&);
int print_eigenvalues( char* desc, int n, double* wr, double* wi );
int print_matrix( char* desc, int m, int n, complex<double> *a, int lda );
int print_int_vector( char* desc, int n, int* a );
int print_complex_vector(vector<vector<complex<double> > >);
#endif
