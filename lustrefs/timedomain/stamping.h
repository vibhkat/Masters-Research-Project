#ifndef STAMPING_H
#define STAMPING_H

#include <vector>
#include <complex>
#include <string>
#include <iostream>
#include <fstream>
int read_nodes(std::string ,std::vector<double>& );
int populate_matrices(std::string ,std::vector<std::vector<double> >&,std::vector<std::vector<double> >& ,std::vector<std::vector<double> >&,std::vector<std::vector<double> >& ,double ,double );
int addzeroes(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,char , double,double ,int &);
int insertvaluesRG(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double,double,double);
int insertvaluesC(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double,double,double);
int insertvaluesL(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double,double,double);
int insertvaluesZ(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double,double,double,double,double);
int insertvaluesE(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double,double,double,double,double);
int insertvaluesV(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double ,double,std::vector<double> ,std::vector<double> ,int &, double, double );
int insertvaluesJ(std::vector<std::vector<double> > &,std::vector<std::vector<double> > &,double ,double ,std::vector<double>,std::vector<double> ,int &, double, double);

struct netlist_t
{
std::vector<std::string>name;
std::vector<double>values;
};


#endif
