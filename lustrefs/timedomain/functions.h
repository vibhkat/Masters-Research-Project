
#ifndef FUNCTIONS_H
#define FUNCTIONS_H

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
int print_2D(std::vector<std::vector<double> >);
int print_1D(std::vector<double> );
int vectoraddition(std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > &);
int vectormultiplication(std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > &);
int AinverseB(std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > &);
extern "C" int dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
int original_solver(std::string,std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > , double , int ,std::vector<double> , std::vector<double> );
int QR_factorization(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &);
int vector_transpose(std::vector<std::vector<double> > ,std::vector<std::vector<double> > &);
int blockArnoldi(std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > ,std::vector<std::vector<double> > &,int );
int reduced_MATRIX_generator(std::vector<std::vector<double> >,std::vector<std::vector<double> >,std::vector<std::vector<double> >,std::vector<std::vector<double> >,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &);
int reduced_solver(std::string,std::vector<std::vector<double> >,std::vector<std::vector<double> >,std::vector<std::vector<double> >,std::vector<std::vector<double> >,double,int,std::vector<double>,std::vector<double>,std::vector<std::vector<double> > );




#endif
