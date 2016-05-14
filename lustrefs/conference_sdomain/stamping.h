#ifndef STAMPING_H
#define STAMPING_H
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
int addzeroes(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,char, double ,double ,int &);
int insertvaluesRG(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesC(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesL(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesJ(vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
int insertvaluesV(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
int insertvaluesZ(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
int insertvaluesE(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
int insertvaluesH(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double,double,double);
int read_nodes(string, vector<double>& );
int populate_matrices(string , vector<vector<double> >& , vector<vector<double> >& , vector<vector<double> >& , vector<complex<double> >&);
#endif
