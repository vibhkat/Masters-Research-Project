#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<cctype>
#include<cmath>
#include"stdio.h"
#include"stdlib.h"
#include<sstream>
#include<cstring>
#include"stamping.h"
#include"functions.h"
#define points 1000
#define time_stop 10e-9
#define time_start 0
int main()
{
	double delta_t=(time_stop-time_start)/points;
	std::cout<<"delta_t :"<<delta_t<<std::endl;
	std::vector<std::vector<double> >G_MATRIX;
	std::vector<std::vector<double> >C_MATRIX;
	std::vector<std::vector<double> >B_MATRIX,u_MATRIX;

	std::vector<double> input_node,output_node;



	read_nodes("./text_files/input_terminal.txt", input_node);
	read_nodes("./text_files/output_terminal.txt", output_node);
	std::cout<<std::endl<<"input nodes :"<<std::endl;
	for(int i=0;i<input_node.size();++i)
	{
		std::cout<<input_node[i]<<" ";
	}
	std::cout<<std::endl<<"output nodes :"<<std::endl;  
	for(int i=0;i<output_node.size();++i)                         
	{
		std::cout<<output_node[i]<<" ";                            
	}

	std::cout<<std::endl<<"Generating the individual Matrices"<<std::endl;
	populate_matrices("./text_files/time_multi_10.txt",G_MATRIX,C_MATRIX, B_MATRIX,u_MATRIX,points,delta_t);
	std::cout<<" THE number of variable are : "<<G_MATRIX.size()<<std::endl;
	std::cout<<std::endl<<"DONE Generating the individual Matrices"<<std::endl;
	std::cout<<"Solving the original model "<<std::endl;


	//solving the original model

	original_solver("./text_files/original_10.txt",G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,delta_t,points,input_node,output_node);
	std::cout<<" DONE solving the original mode"<<std::endl;
	//computing orthonormal matric and also solving reduced model       
	char answer='y';
	while(answer=='y')
	{
		std::vector<std::vector<double> > X_MATRIX,G_cap,C_cap,B_cap;
		int q_val;
		std::cout<<"enter the value of q : ";
		std::cin>>q_val;
		blockArnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX,q_val);
		std::cout<<"the number of rows in X_MATRIX are: "<<X_MATRIX.size()<<" the number of columns are :"<<(X_MATRIX[0]).size()<<std::endl;
		//generating the reduced matrices
		std::cout<<"generating the reduced Matrices of orginal matrices"<<std::endl;
		reduced_MATRIX_generator(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX,G_cap,C_cap,B_cap);		
		std::cout<<"done generating the reduced Matrices of original matrices"<<std::endl;
		std::cout<<"the number of X reduced are :"<<G_cap.size()<<std::endl; 
		std::cout<<"Solving the reduced model "<<std::endl;
		

std::string X_str_col;
std::ostringstream ss1;
ss1<<(X_MATRIX.at(0)).size();
X_str_col=ss1.str();
std::string MOR_filename="./text_files/MOR("+X_str_col+").txt";
reduced_solver(MOR_filename,G_cap,C_cap,B_cap,u_MATRIX,delta_t,points,input_node,output_node,X_MATRIX);	
		std::cout<<"done solving the reduced model"<<std::endl;


		std::cout<<"do you want to continue for new value of Q (y/n)---- -:";
		std::cin>>answer;


	}



	return 0;
}

