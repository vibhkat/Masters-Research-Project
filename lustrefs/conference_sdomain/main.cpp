/*First enter RGC irrespective of the order,then enter L,the dependent sources and then dependent sources */
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
#include"print.h"
#include"stamping.h"
#include"functions.h"
#define pi 3.14
#define frequency_max 10e9
#define frequency_min 1e3
#define points 1000 




int main(int argc, char *argv[])
{      
	double f_max=frequency_max;
	double f_min=frequency_min;
	double pts=points;
	int count_s=0;
	vector<vector<double> > X_big;


	for(int num=1;num<argc;++num)
		//	char cont_1='y';
		//	while(cont_1=='y')
	{
		//if(cont_1=='y')
		++count_s;
		string input_file=argv[num];
		cout<<" Entering the file : "<<input_file<<endl;
		original_plus_MOR_giving_X_big(X_big,f_max,f_min,points,input_file);	
		//	cout<<" DO YOU WANT TO ENTER THE FILES FOR OTHER CIRCUIT (Y/N)----: ";
		//	cin>> cont_1;
	}// main while loop for different cicuits
	file_write("SERIAL_X_BIG.txt",X_big);
	cout<<"X_big written in SERIAL_X_BIG.txt"<<endl;
	cout<<" number of circuits ran are : "<<count_s<<endl;
	cout<<"the number of rows in X_big are : "<<X_big.size()<<" the number of columns in X_big are : "<<(X_big.at(0)).size()<<endl;
	vector<vector<double> > X_svd;
	vector<double> sigma;
	cout<<"performing SVD................"<<endl;
	double t1=clock();
	SVD(X_svd,sigma,X_big);
	file_write("SERIAL_X_SVD.txt",X_svd);
	double t2=clock();
	cout<<" TIME taken to compute SVD is : "<<(t2-t1)/double(CLOCKS_PER_SEC) << " seconds" <<endl;

	cout<<"the number of rows in X_svd are : "<<X_svd.size()<<" the number of columns in X_svd are : "<<(X_svd.at(0)).size()<<endl;
	cout<<"the number of elements in SIGMA are : "<<sigma.size()<<endl;

	//writing the ratio in a file	
	fstream write("sigma.txt",fstream::out);
	if(write.is_open())
	{
		for(int i=0;i<sigma.size();++i)
		{
			write<<i<<" "<<sigma[i]<<" "<<sigma[i]/sigma[0]<<endl;
		}
		write.close();
	}
	else cout<<"enable to open file"<<endl;

	char answer='y';
	while(answer=='y')
	{

		double thr;
		int count_a=0;
		cout<<"Enter the threshold for sigma matrix :";
		cin>>thr;
		for(int i=1;i<sigma.size();++i)
		{
			if(sigma[i]/sigma[0]>thr || sigma[i]/sigma[0]==thr)
			{
				++count_a;
			}
		}
		cout<<" the no of columns are : "<<count_a<<endl;
		vector<vector<double> > X_MATRIX;
		retrun_X_MATRIX(X_MATRIX,X_svd,count_a);
		cout<<"the number of rows in X_MATRIX are : "<<X_MATRIX.size()<<" the number of columns in X_MATRIX are : "<<(X_MATRIX.at(0)).size()<<endl;
		original_plus_MOR_sending_X_MATRIX(X_MATRIX,f_max,f_min,points);
		X_MATRIX.clear();
		cout<<"do you want to continue for other vlaue of threshold (y/n) :";
		cin>>answer;
	}
	cout<<"end of program"<<endl;
	return 0;
}



















