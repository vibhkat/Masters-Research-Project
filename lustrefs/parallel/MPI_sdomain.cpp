
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
#include<mpi.h>
#include<iomanip>
using namespace std;
void print(vector<vector<double> >&);
void printcomplexvector(vector<complex<double> >&);
void addzeroes(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,char, double ,double ,int &);
void insertvaluesRG(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesC(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesL(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesJ(vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
void insertvaluesV(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
void insertvaluesZ(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
void insertvaluesE(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
void insertvaluesH(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double,double,double);

#define pi 3.14
#define frequency_max 10e9
#define frequency_min 1e3
#define points 1000 
extern "C" void zgesv_( int *n , int *nrhs , complex<double> *a , int *lda , int *ipiv ,complex<double> *b , int *ldb , int *info  );
void print_matrix( char* desc, int m, int n, complex<double> *a, int lda );
void print_int_vector( char* desc, int n, int* a );
extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void AinverseB(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
void QR_factorization(vector<vector<double> > &,vector<vector<double> > &);
void vectormultiplication(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
void vector_transpose(vector<vector<double> > &,vector<vector<double> > &);
void blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int & ,vector<double> &);
extern  "C" void dgeev_( char* jobvl, char* jobvr, int* n, double* a,int* lda,  double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info );
int multi_blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int &,vector<double> &);
void print_eigenvalues( char* desc, int n, double* wr, double* wi );
int block_inverse(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,double );
int complex_vector_multiplication(vector<vector<complex<double> > > ,vector<vector<complex<double> > > ,vector<vector<complex<double> > > &);
int print_complex_vector(vector<vector<complex<double> > >);
int multi_Arnoldi(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,int, double );
int SVD(vector<vector<double> > &,vector<double> &,vector<vector<double> >);
extern "C" void dgesvd_(char* jobu,char* jobvt,int* m,int* n,double* a,int* lda,double* s,double* u,int* ldu,double* vt,int* ldvt,double* work,int* lwork,int* info );
int retrun_X_MATRIX(vector<vector<double> > & ,vector<vector<double> >,int );
int original_plus_MOR_sending_X_MATRIX(vector<vector<double> >,double,double, double, string,int, string);
int read_nodes(string, vector<double>& );
int populate_matrices(string , vector<vector<double> >& , vector<vector<double> >& , vector<vector<double> >& , vector<complex<double> >&);

int main(int argc, char *argv[])
{       



	int numbernodes, mynode,len ;
	char hostname[MPI_MAX_PROCESSOR_NAME];

	MPI_Status status;
	MPI_Request request, request2;
	int rc=MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) 
	{
		cout<<"Error starting MPI program. Terminating.\n";
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
	MPI_Comm_size(MPI_COMM_WORLD, &numbernodes);
	MPI_Get_processor_name(hostname, &len);
	printf ("Number of tasks= %d My rank= %d Running on %s\n", numbernodes,mynode,hostname);
	string s_mynode;
	ostringstream ss;
	ss<<mynode+1;
	s_mynode=ss.str();

	cout<<" string node is"<<s_mynode<<endl;

	double f_max=frequency_max;
	double f_min=frequency_min;
	double pts=points;
	double interval;
	interval=(frequency_max-frequency_min)/(points-1);
	vector<complex<double> > u_MATRIX;
	//TAKING OUT EACH STRING AND CONVERTING IT TO RESPECTIVE VALUES
	string line,s,com1,com2;
	char new_word,old_word='R';
	vector<double>values;
	vector<string>name;	
	vector<string>read; //WATCH FOR SPACE IN INPUT
	vector<vector<double> >G_MATRIX;
	vector<vector<double> >C_MATRIX;
	vector<vector<double> >B_MATRIX;
	int column=0,previous_section=1,new_section;
	vector<vector<double> > L_columns;
	vector<double>L_temp;	
	G_MATRIX.clear();
	C_MATRIX.clear();
	B_MATRIX.clear();


	//reading the input and output node from files created from pul.c code        
	vector<double> input_node,output_node;
	double temp_node;
	fstream read1("input_terminal.txt",fstream::in);
	while(read1>>temp_node)
		input_node.push_back(temp_node);
	read1.close();
	fstream read2("output_terminal.txt",fstream::in);
	while(read2>>temp_node)
		output_node.push_back(temp_node);
	read2.close();
	//printing the nodes
	for(int i=0;i!=input_node.size();++i)
		cout<<endl<<"PROCESS NUMBER : "<<mynode<<setw(mynode)<<"---"<<input_node[i]<<" ";
	cout<<endl;
	for(int i=0;i!=output_node.size();++i)
		cout<<endl<<"PROCESS NUMBER : "<<mynode<<setw(mynode)<<"---"<<output_node[i]<<" ";
	cout<<endl;
	MPI::COMM_WORLD.Barrier();
	cout<<endl<<" PROCESS NUMBER : "<< mynode <<"---Printing the element passing from main is "<< argv[mynode+1]<<endl;
	string filename=argv[mynode+1];
	ifstream myfile(filename.c_str());	
	cout<<endl<<"PROCESS NUMBER : "<<mynode<<setw(mynode)<<"---Generating the individual Matrices"<<endl;
	while(getline(myfile,line))/* && !line.empty())*/{

		for(int i=0;i!=(line.size());++i){
			if(isspace(line[i])&& !s.empty()){
				read.push_back(s);
				s.clear();
			}
			else if(!isspace(line[i])){
				s+=line[i];
			}
		}
		read.push_back(s);
		s.clear();//till here stored the strings from line in vector called read.
		//now reading that vector and taking out which are numbers and which are string.		
		for(int noofelements=0;noofelements!=read.size();++noofelements)
		{

			string &element=read[noofelements];
			//determining the which is string and which is value
			int count=0;
			for(int i=0;i!=element.size();++i)
			{
				if(isdigit(element[i])||ispunct(element[i])||(isdigit(element[i-1])&&element[i]=='e')||(isdigit(element[i-1])&&element[i]=='E'))
					++count;
			}
			if(count==element.size())
			{
				//taking each values and determining where is comma and word e or E
				int comma_place;
				int comma_counter=0;
				int e_counter=0;
				vector<int>e_place;
				for(int i=0;i!=element.size();++i)
				{
					if(element[i]==','){
						comma_place=i;
						++comma_counter;
					}
					else if(element[i]=='e'||element[i]=='E')
					{
						++e_counter;
						e_place.push_back(i);
					}
				}
				if(comma_counter==0&&e_counter==0)
				{
					double f;
					istringstream(element)>>f;
					values.push_back(f);
					values.push_back(0);
				}
				else if(comma_counter==1&&e_counter==0)
				{
					for(int z=0;z<comma_place;++z)
						s+=element[z];                                
					double f;
					istringstream(s)>>f;
					values.push_back(f);
					s.clear();
					for(int z=(comma_place+1);z!=element.size();++z)
						s+=element[z];
					double g;
					istringstream(s)>>g;
					values.push_back(g);
					s.clear();
				}

				else if(comma_counter==0&& e_counter==1)
				{
					double number,power,scientific;
					for(int z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number;
					s.clear();
					for(int z=(e_place[0]+1);z!=element.size();++z)
						s+=element[z];
					istringstream(s)>>power;
					s.clear();
					scientific=number*pow(10,power);
					values.push_back(scientific);
					values.push_back(0);
				}

				else if(comma_counter==1&& e_counter==2)
				{
					double number1,power1,scientific1,number2,power2,scientific2;
					for(int z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number1;
					s.clear();
					for(int z=(e_place[0]+1);z!=comma_place;++z)
						s+=element[z];
					istringstream(s)>>power1;
					s.clear();
					scientific1=number1*pow(10,power1);
					values.push_back(scientific1);
					for(int z=comma_place+1;z!=e_place[1];++z)
						s+=element[z];
					istringstream(s)>>number2;
					s.clear();
					for(int z=(e_place[1]+1);z!=element.size();++z)
						s+=element[z];
					istringstream(s)>>power2;
					s.clear();
					scientific2=number2*pow(10,power2);
					values.push_back(scientific2);



				}
				else if(comma_counter==1&& e_counter==1)
				{
					if(comma_place<e_place[0])
					{
						double f,number,power,scientific;
						for(int z=0;z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>f;
						s.clear();
						values.push_back(f);
						for(int z=comma_place+1;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(int z=(e_place[0]+1);z!=element.size();++z)
							s+=element[z];
						istringstream(s)>>power;
						scientific=number*pow(10,power);
						values.push_back(scientific);
					}
					else if(comma_place>e_place[0])
					{
						double f,number,power,scientific;
						for(int z=0;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(int z=(e_place[0]+1);z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>power;
						s.clear();						
						scientific=number*pow(10,power);
						values.push_back(scientific);
						for(int z=comma_place+1;z!=element.size();++z)
							s+=element[z];
						istringstream(s)>>f;
						s.clear();
						values.push_back(f);
					}

				}
			}
			else
				name.push_back(element);
		}
		addzeroes(G_MATRIX,C_MATRIX,B_MATRIX,name[0][0],values[0],values[2],column);
		if(name[0][0]=='G'|| name[0][0]=='g')	      
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		else if(name[0][0]=='R' || name[0][0]=='r'){
			double newvalue=(1/values[4]);
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],newvalue);
		}
		else if(name[0][0]=='C' || name[0][0]=='c'){
			insertvaluesC(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}		
		else if(name[0][0]=='L' || name[0][0]=='l'){
			insertvaluesL(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}
		else if(name[0][0]=='J' || name[0][0]=='j')
			insertvaluesJ(B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		else if(name[0][0]=='V' || name[0][0]=='v')
		{
			insertvaluesV(G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		}		
		else if(name[0][0]=='Z' || name[0][0]=='z')
			insertvaluesZ(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='E' || name[0][0]=='e')
			insertvaluesE(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='H' || name[0][0]=='h')
			insertvaluesH(G_MATRIX,C_MATRIX,L_columns,values[0],values[2],values[4],values[6],values[8],values[10],values[12]);	
		values.clear();
		name.clear();
		read.clear(); // to clear the read vector*/

	}
	;
	cout<<endl<<"PROCESS NUMBER : "<<mynode<<setw(mynode)<<"---Done generating individual matrices"<<endl;
	MPI::COMM_WORLD.Barrier();
	cout<<endl<<"PROCESS NUMBER : "<<mynode<<"--THE NUMBER OF VARIABLE ARE :"<<G_MATRIX.size()<<endl;
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(int i=0;i!=(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	//converting the B VECTOR MATRIX to B ARRAY MATRIX
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(int i=0;i!=(B_MATRIX).size();++i)
		for(int j=0;j!=(u_MATRIX).size();++j)
			B_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_MATRIX[i][j],0);

	//multiplying the B_COMPLEXARRAY AND uarraymatrix
	complex<double> Bucomplexmultiple[(B_MATRIX).size()][1];
	for(int i=0;i!=(B_MATRIX).size();++i)
		for(int j=0;j!=(u_MATRIX).size();++j)
			Bucomplexmultiple[i][0]+=B_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];
	MPI::COMM_WORLD.Barrier();
	string out_file1="./result/ORIGINAL_s"+s_mynode+".txt";
	cout<<endl<<"PROCESS NUMBER : "<<mynode<<"---output filename for orignial solution= "<<out_file1<<endl;
	MPI::COMM_WORLD.Barrier();
	double t1=MPI_Wtime();
	//comment for original solving starts here	
	fstream result11(out_file1.c_str(),fstream::out);
	if(result11.is_open())
	{
		double frequency;
		int timer1=0;
		while(timer1<points)
		{
			cout<<endl<<"MYNODE :"<<mynode<<setw(mynode)<<"original"<<timer1+1<<endl;	
			frequency=frequency_min+(interval*timer1);
			result11<<frequency<<'\t';
			//combining G_MATRIX and C_MATRIX into GplussC_MATRIX
			complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
			for(int i=0;i!=(G_MATRIX).size();++i){
				for(int j=0;j!=(G_MATRIX).size();++j)
				{
					GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
				}
			}
			//copying elements of GplusC_MATRIX in a one dimension array a
			complex<double> a[(G_MATRIX.size())*(G_MATRIX.size())];
			int k=0;
			for(int i=0;i!=(G_MATRIX).size();++i)
				for(int j=0;j!=(G_MATRIX).size();++j)
					a[k++]=GplussC_MATRIX[j][i];

			//copying elements of Bucomplexmultiple into b one dimension matrix
			complex<double> b[B_MATRIX.size()];
			for(int i=0;i!=(B_MATRIX).size();++i)
				b[i]=Bucomplexmultiple[i][0];

			int n=G_MATRIX.size();
			int nrhs=1;
			int lda=n;
			int ldb=n;
			int info;
			int ipiv[n];
			zgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
			if( info > 0 ) {
				printf( "The diagonal element of the triangular factor of A,\n" );
				printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
				printf( "the solution could not be computed.\n" );
				return( 1 );
			}
			complex<double> original_solution[n][1];
			for(int i=0;i!=n;++i)
				for(int j=0;j!=nrhs;++j)
					original_solution[i][j]=b[i+j*ldb];
			for(int i=0;i!=G_MATRIX.size();++i)
			{
				for(int k=0;k!=output_node.size();++k)
				{
					if(i==(output_node[k]-1))
					{
						if((original_solution[i][0]).imag()>=0)
							result11<<(original_solution[i][0]).real()<<'+'<<(original_solution[i][0]).imag()<<'i'<<'\t';
						else
							result11<<(original_solution[i][0]).real()<<(original_solution[i][0]).imag()<<'i'<<'\t';

					}

				}
			}
			result11<<';'<<endl;

			++timer1;
		}//while loop of frequency ends here


	}
	else cout<<"unable to open file"<<endl;

	cout<<endl<<"PROCESS NUMBER : "<<mynode<<"--- done computing Ax=b of original solution for different frequencies and writing for MATLAB READING in= "<<out_file1<<endl;
	//comment for original model solving ends here

	double t2=MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);




	//function to generate the X matrices containing the genrerated perpendicular columns for MOR
	cout<<endl<<"STARTING THE MOR FROM HERE"<<endl;	
	//char answer='y';//while loop starts here
	//while(answer=='y')
	//{	

	int expansion=1;
	cout<<"enter the number of expansion points :"<<expansion<<endl;
	//	cin>>expansion;
	double abcd[expansion],efgh[expansion],expan_freq;
	int q_val;
	for(int i=0;i<expansion;++i)

	{

		expan_freq=0;
		cout<<"enter the "<<i+1<<" expansion frequency: "<<expan_freq<<endl;


		//	cout<<"enter the "<<i+1<<" expansion frequency: ";
		//	cin>>expan_freq;
		abcd[i]=expan_freq;
		q_val=414;
		cout<<"enter the value of q : "<<q_val<<endl;

		//		cout<<"enter the value of q : ";
		//		cin>>q_val;
		efgh[i]=q_val;
	}




	vector<vector<double> > X_MATRIX;
	vector<double>range;	
	cout<<"computing the orthonormal matrices using block Arnoldi algorithm"<<endl;

	for(int count_Q=0;count_Q<expansion;++count_Q)
	{

		double expansion_frequency;
		expansion_frequency=abcd[count_Q];
		int Q_value;
		Q_value=efgh[count_Q];
		if(count_Q==0)
			blockArnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX,Q_value,range);
		else{
			cout<<"entering...................................................."<<endl;
			//			cout<<X_MATRIX.size()<<"................................."<<(X_MATRIX.at(0)).size()<<endl;
			multi_Arnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX,Q_value,expansion_frequency);
			cout<<"out..........................................................."<<endl;
		}


		cout<<"MYNODE :"<<mynode<<" the number of rows in "<< count_Q+1<<" are: "<<X_MATRIX.size()<<" the number of columns are :"<<(X_MATRIX.at(0)).size()<<endl;

	}
	cout<<"the number of rows in X_MATRIX are: "<<X_MATRIX.size()<<" the number of columns are :"<<(X_MATRIX[0]).size()<<endl;
	//checking each coulumvn is orthnormal to each other		
	//generating the reduced matrices
	cout<<"generating the reduced Matrices of orginal matrices"<<endl;
	vector<vector<double> > X_transpose,multiple,G_cap,C_cap,B_cap;
	vector_transpose(X_MATRIX,X_transpose);
	vectormultiplication(X_transpose,X_MATRIX,multiple);
	multiple.clear();
	vectormultiplication(X_transpose,G_MATRIX,multiple);
	vectormultiplication(multiple,X_MATRIX,G_cap);
	multiple.clear();
	vectormultiplication(X_transpose,C_MATRIX,multiple);
	vectormultiplication(multiple,X_MATRIX,C_cap);
	multiple.clear();
	vectormultiplication(X_transpose,B_MATRIX,B_cap);

	cout<<"done generating the reduced Matrices of orginal matrices"<<endl;	
	X_transpose.clear();
	vector_transpose(X_MATRIX,X_transpose);
	cout<<"the number of X reduced are :"<<G_cap.size()<<endl; 
	string X_str_col;
	ostringstream ss1;
	ss1<<(X_MATRIX.at(0)).size();
	X_str_col=ss1.str();
	string out_file2="./result/MOR_s"+s_mynode+"("+X_str_col+").txt";
	cout<<endl<<"PROCESS NUMBER : "<<mynode<<"---output filename for MOR solution= "<<out_file2<<endl;
	MPI::COMM_WORLD.Barrier();
	double t3=MPI_Wtime();


	fstream result22(out_file2.c_str(),fstream::out);
	if(result22.is_open() )
	{
		double frequency;
		int timer2=0;	
		while(timer2<points)
		{
			cout<<"computing MOR at "<<timer2+1<<" frequency point"<<endl;			
			frequency=frequency_min+(interval*timer2);

			result22<<frequency<<'\t';
			// Combining G_cap and C_cap into G_capplussC_cap
			complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
			for(int i=0;i!=(G_cap).size();++i){
				for(int j=0;j!=(G_cap).size();++j)
				{
					G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
				}
			}
			//copying elements of GcapplusCcap in a one dimension array acap
			complex<double> a_cap[(G_cap.size())*(G_cap.size())];
			int k_cap=0;
			for(int i=0;i!=(G_cap).size();++i)
				for(int j=0;j!=(G_cap).size();++j)
					a_cap[k_cap++]=G_capplussC_cap[j][i];

			//converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX
			complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
			for(int i=0;i!=(B_cap).size();++i)
				for(int j=0;j!=(u_MATRIX).size();++j)
					B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);

			//multiplying the B_COMPLEXARRAY AND uarraymatrix
			complex<double> B_capucomplexmultiple[(B_cap).size()][1];
			for(int i=0;i!=(B_cap).size();++i)
				for(int j=0;j!=(u_MATRIX).size();++j)
					B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];

			//copying elements of B_capucomplexmultiple into b_cap one dimension matrix
			complex<double> b_cap[(B_cap).size()];
			for(int i=0;i!=(B_cap).size();++i)
				b_cap[i]=B_capucomplexmultiple[i][0];

			int n_cap=G_cap.size();
			int nrhs_cap=1;

			int lda_cap=n_cap;
			int ldb_cap=n_cap;
			int info_cap;
			int ipiv_cap[n_cap];
			zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
			if( info_cap > 0 ) {
				printf( "The diagonal element of the triangular factor of A,\n" );
				printf( "U(%i,%i) is Jero, so that A is singular;\n", info_cap, info_cap );
				printf( "the solution could not be computed.\n" );
				return( 1 );
			}


			complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
			for(int i=0;i!=X_MATRIX.size();++i)
				for(int j=0;j!=(X_MATRIX[0]).size();++j)
					X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

			//converting b_cap from inverse into x_cap 2 dimension array
			complex<double> x_cap[G_cap.size()][1];

			for(int i=0;i!=n_cap;++i)
				for(int j=0;j!=nrhs_cap;++j)
					x_cap[i][j]=b_cap[i+j*ldb_cap];
			complex<double>MORsolutioncompare[X_MATRIX.size()][1];
			for(int i=0;i!=X_MATRIX.size();++i)
				for(int j=0;j!=(X_MATRIX[0]).size();++j)
					MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];
			for(int i=0;i!=G_MATRIX.size();++i)
			{
				for(int k=0;k!=output_node.size();++k)
				{
					if(i==(output_node[k]-1)){
						if((MORsolutioncompare[i][0]).imag()>=0)              
							result22<<(MORsolutioncompare[i][0]).real()<<'+'<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';
						else
							result22<<(MORsolutioncompare[i][0]).real()<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';

					}
				}
			}
			result22<<';'<<endl;
			++timer2;			
		}
		result22.close();

	}
	else cout<<"unable top open file"<<endl;

	cout<<"done computing Ax=b of MOR solution for different frequencies and writing for MATLAB READING in "<<out_file2<<endl;	
	MPI::COMM_WORLD.Barrier();
	double t4=MPI_Wtime();


	cout<<"Starting combining from here-------------------------"<<endl;
	int row_X=X_MATRIX.size();
	int column_X=(X_MATRIX.at(0)).size();
	double send_X[row_X][column_X];
	//converting X_MATRIX vector into array
	cout<<"Converting X_MATRIX vector into array============================="<<endl;
	for(int i=0;i<row_X;++i)
	{
		for(int j=0;j<column_X;++j)
		{
			send_X[i][j]=X_MATRIX[i][j];
		}
	}

	cout<<"Done converting X_Matrix vector into array.............................."<<endl;

	MPI_Barrier(MPI_COMM_WORLD);
	vector<vector<double> > X_big;

	if(mynode!=0)
	{
		cout<<"Sending X matrix from "<<mynode <<" node to 0 node"<<endl;
		MPI_Send(&row_X,1,MPI_INT,0,1,MPI_COMM_WORLD);
		cout<<"row send"<<endl;				
		MPI_Send(&column_X,1,MPI_INT,0,1,MPI_COMM_WORLD);
		cout<<"column send"<<endl;
		int dimension=row_X*column_X;
		cout<<"dimension :"<<dimension<<endl;
		MPI_Send(&send_X,dimension,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
		cout<<" 2D array send rank "<<mynode<<" OK"<<endl;
	}
	if(mynode==0)
	{
		vector<double> temp_X;
		if(X_big.empty())
		{
			for(int i=0;i<X_MATRIX.size();++i)
			{
				for(int j=0;j<(X_MATRIX.at(0)).size();++j)
				{
					temp_X.push_back(X_MATRIX[i][j]);
				}
				X_big.push_back(temp_X);
				temp_X.clear();
			}
		}
		for(int i=1;i<numbernodes;++i)
		{
			cout<<"receiving X_MATRIX from "<<i<<" node"<<endl;
			int row_rec,column_rec;
			MPI_Recv(&row_rec,1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
			MPI_Recv(&column_rec,1,MPI_INT,i,1,MPI_COMM_WORLD,&status);
			double X_rec[row_rec][column_rec];
			MPI_Recv(&X_rec,row_rec*column_rec,MPI_DOUBLE,i,1,MPI_COMM_WORLD,&status);

			vector<double> temp_X;
			if(X_big.empty())
			{
				for(int i=0;i<row_rec;++i)
				{
					for(int j=0;j<column_rec;++j)
					{
						temp_X.push_back(X_rec[i][j]);
					}
					X_big.push_back(temp_X);
					temp_X.clear();
				}
			}
			else
			{
				for(int i=0;i<row_rec;++i)
				{
					for(int j=0;j<column_rec;++j)
					{
						(X_big[i]).push_back(X_rec[i][j]);
					}
				}
			}

			cout<<"done sending X_matrix from each node to 0 node"<<endl;
		}
	}

	vector<vector<double> > X_svd;
	vector<double> sigma;
	MPI_Barrier(MPI_COMM_WORLD);
	if(mynode==0)
	{

		cout<<"0: the number of rows in X_big are : "<<X_big.size()<<" the number of columns in X_big are : "<<(X_big.at(0)).size()<<endl;
		cout<<"0: performing SVD................"<<endl;
		double t5=MPI_Wtime();

		SVD(X_svd,sigma,X_big);
		double t6=MPI_Wtime();
		cout<<"TIME TAKEN TO PERFORM SVD :"<<t6-t5<<endl;
		cout<<"0: the number of rows in X_svd are : "<<X_svd.size()<<" the number of columns in X_svd are : "<<(X_svd.at(0)).size()<<endl;
		cout<<"0: the number of elements in SIGMA are : "<<sigma.size()<<endl;
		fstream write("./result/sigma.txt",fstream::out);
		if(write.is_open())
		{
			for(int i=0;i<sigma.size();++i)
			{
				write<<i<<" "<<sigma[i]<<" "<<sigma[i]/sigma[0]<<endl;
			}
			write.close();
		}
		else cout<<"enable to open file"<<endl;
	}

	int row_S,column_S,sigma_size;
	if(mynode==0)
	{
		row_S=X_svd.size();
		column_S=(X_svd.at(0)).size();
		sigma_size=sigma.size();
		cout<<"0: row printt ---------------------------"<<row_S<<" column print ==================="<<column_S<<" sigma size :"<<sigma_size<<endl;

	}
	MPI_Bcast(&row_S,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&column_S,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&sigma_size,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	cout<<"MYNODE :"<<mynode<<", ROW :"<<row_S<<" "<<", COLUMN :"<<column_S<<", sigma size "<<sigma_size<<endl;
	double X_SVD_array[row_S][column_S];
	double sigma_array[sigma_size];
	if(mynode==0)
	{
		//converting X_MATRIX vector into array
		cout<<"0: Converting X_svd vector into array"<<endl;
		for(int i=0;i<row_S;++i)
		{
			for(int j=0;j<column_S;++j)
			{
				X_SVD_array[i][j]=X_svd[i][j];
			}
		}
		cout<<"0: converting sigma vector into array"<<endl;
		for(int i=0;i<sigma.size();++i)
		{
			sigma_array[i]=sigma[i];
		}
		cout<<"0: done converting X_svd into array"<<endl;


	}
	MPI_Bcast(&X_SVD_array,row_S*column_S,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sigma_array,sigma_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if(mynode!=0)
	{
		cout<<"MYNODE :"<<mynode<<" receiving in node "<<endl;

		vector<double> temp;
		for(int i=0;i<row_S;++i)
		{
			for(int j=0;j<column_S;++j)
			{
				temp.push_back(X_SVD_array[i][j]);
			}
			X_svd.push_back(temp);
			temp.clear();
		}
		for(int i=0;i<sigma_size;++i)
		{
			sigma.push_back(sigma_array[i]);    
		}
		cout<<"MYNODE :"<<mynode<<"  ,the number of rows in X_svd are : "<<X_svd.size()<<" ,the number of columns in X_svd are : "<<(X_svd.at(0)).size()<<" ,sigm size : "<<sigma.size()<<endl;
		cout<<"MYNODE :"<<mynode<<" done receiving in node "<<endl;

		cout<<"MYNODE :"<<mynode<<" DONE sending X_SVD_array and sigma_ array and storing in corresponding vectors---------------------------------"<<endl;
	}
	MPI_Barrier(MPI_COMM_WORLD);


	double thr=0;
	char ans_thr='y';
	while(ans_thr=='y')
	{

		if(mynode==0)
		{
			cout<<"ENTER THE THRESHOLD FOR COLUMNS: ";
			cin>>thr;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&thr,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		cout<<"MYNODE : "<<mynode<<" threshold received is : "<<thr<<endl;
		MPI_Barrier(MPI_COMM_WORLD);
		int count_a=0;
		for(int i=0;i<sigma.size();++i)
		{
			if(((sigma[i]/sigma[0])>thr) ||((sigma[i]/sigma[0])==thr) )
			{
				++count_a;
			}
		}
		vector<vector<double> > X_MATRIX;
		retrun_X_MATRIX(X_MATRIX,X_svd,count_a);
		cout<<"MYNODE :"<<mynode<<" the number of rows in X_MATRIX are : "<<X_MATRIX.size()<<", the number of columns in X_MATRIX are : "<<(X_MATRIX.at(0)).size()<<endl;
		original_plus_MOR_sending_X_MATRIX(X_MATRIX,f_max,f_min,points,filename,mynode,s_mynode);
		X_MATRIX.clear();
		MPI_Barrier(MPI_COMM_WORLD);


		if(mynode==0)
		{
			cout<<"do you want to continue for the new value of threshold(y/n) :";
			cin>>ans_thr;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&ans_thr,1,MPI_CHAR,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}



	MPI_Barrier(MPI_COMM_WORLD);
	cout<<"time taken to solve original model :"<<t2-t1<<endl;
	cout<<"time taken to solve reduced model :"<<t4-t3<<endl;
	MPI_Finalize();
	cout<<"end of program"<<endl;

	return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, complex<double> *a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%f,%f)", a[i+j*lda].real(), a[i+j*lda].imag() );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}




/*Function for printing the 2 dimension vector*/
void print( vector<vector<double> > &p)
{
	for(int i=0; i!=(p.size());++i){
		for(int j=0; j!=((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;}
}

/* Function for printing one dimension complex vector*/
void printcomplexvector(vector<complex<double> >&p)
{
	for(int i=0; i!=(p.size());++i){
		// for(int j=0; j!=((p[i]).size());++j)
		cout<<p[i]<<" ";
		cout<<";"<<endl;}
}

/* INSERTING ZEROES IN G,C,B MATRIC*/
void addzeroes(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,char S, double nodeI,double nodeJ,int &count)
{
	int size=(nodeI>nodeJ)?nodeI:nodeJ;
	int size1=G.size();//or size1=C.size,had to add this syntax because when the vectros were not empty,it was not able to take the size of the vectors.	
	int size2=B.size();	
	if(S=='R'||S=='r'||S=='G'||S=='g'||S=='C'||S=='c'){
		if(G.empty() && C.empty() && B.empty()  ){
			G=vector<vector<double> >(size,vector<double>(size,0));
			C=vector<vector<double> >(size,vector<double>(size,0));
			B=vector<vector<double> >(size,vector<double>(1,0));}
		else if(!G.empty() && !C.empty() && !B.empty() && size>G.size() && size>C.size()){
			for(int i=0; i!=size1;++i){
				for(int j=0; j!=(size-(G.size())) ;++j){
					G[i].push_back(0);
					C[i].push_back(0);
				}}
			for(int k=0; k!=(size-size1) ; ++k){
				G.push_back(vector<double>(size,0));
				C.push_back(vector<double>(size,0));
				B.push_back(vector<double>(1,0));
			}
		}

	}


	else if(S=='L'||S=='l'){
		for(int i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);}
			G.push_back(vector<double>((G.size()+1),0));
			C.push_back(vector<double>((C.size()+1),0));
			B.push_back(vector<double>(1,0));
	}

	else if(S=='J'||S=='j'){
		int s=0;
		for(int i=0;i!=B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(int i=0;i!=B.size();++i)
				B[i].push_back(0);
		}
	}
	else if(S=='V'||S=='v'){
		for(int i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
		int s=0;
		for(int i=0;i!=B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(int i=0;i!=B.size();++i)
				B[i].push_back(0);

		}
	}

	else if(S=='E'||S=='e'){
		for(int i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
	}
	else if(S=='H'||S=='h'){
		for(int i=0;i!=size1;++i){
			//	for(int j=0;j!=2;++j){
			G[i].push_back(0);
			C[i].push_back(0);
			//	}
		}
		G.push_back(vector<double>((size1+1),0));
		//	G.push_back(vector<double>((size1+2),0));
		C.push_back(vector<double>((size1+1),0));
		//	C.push_back(vector<double>((size1+2),0));
		B.push_back(vector<double>((count+1),0));
		//	B.push_back(vector<double>((count+1),0));
	}



}
void insertvaluesRG(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if(((i==I)&&(j==J))||((i==J)&&(j==I))&& I>=0 && J>=0)
				G[i][j]=G[i][j]-value;
			else if(((i==I)&&(j==I)&& I>=0)||((i==J)&&(j==J)&& J>=0))
				G[i][j]=G[i][j]+value;
		}
	}
}

void insertvaluesC(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i!=(C.size());++i){
		for(int j=0;j!=((C[i]).size());++j){
			if(((i==I)&&(j==J))||((i==J)&&(j==I))&& I>=0 && J>=0)
				C[i][j]=C[i][j]-value;
			else if(((i==I)&&(j==I)&& I>=0)||((i==J)&&(j==J)&& J>=0))
				C[i][j]=C[i][j]+value;
		}
	}
}

void insertvaluesL(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
{
	C[(C.size())-1][(C.size())-1]=C[(C.size())-1][(C.size())-1]+value;
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if(i==(G.size()-1) && j==I && I>=0)
				G[i][j]=G[i][j]-1;
			else if(i==(G.size()-1) && j==J && J>=0)
				G[i][j]=G[i][j]+1;
			else if(j==(G.size()-1) && i==I && I>=0)
				G[i][j]=G[i][j]+1;
			else if(j==(G.size()-1) && i==J && J>=0)
				G[i][j]=G[i][j]-1;
		}
	}
}

void insertvaluesJ(vector<vector<double> > &B,vector<complex<double> > &U,double nodeI,double nodeJ,double value1,double value2,int &count)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i!=B.size();++i)
	{
		if(i==I && I>=0)
			B[i][count]=B[i][count]+1;
		else if(i==J && J>=0)
			B[i][count]=B[i][count]-1;
	}

	complex<double> b= complex<double>(value1,value2);
	U.push_back(b);


}





void insertvaluesV(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<complex<double> >&U,double nodeI,double nodeJ,double value1,double value2,int &count)
{//cout<<"in the insert function"<<endl;
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if(i==(G.size()-1) && j==I && I>=0)
			{
				//	cout<<"inserted at row "<<i<<" and column "<<j<<endl;
				G[i][j]=G[i][j]-1;}
			else if(i==(G.size()-1) && j==J && J>=0)
			{//	cout<<"inserted at row "<<i<<" and column "<<j<<endl;

				G[i][j]=G[i][j]+1;
			}				
			else if(j==(G.size()-1) && i==I && I>=0)
			{// cout<<"inserted at row "<<i<<" and column "<<j<<endl;

				G[i][j]=G[i][j]+1;}
			else if(j==(G.size()-1) && i==J && J>=0)
			{G[i][j]=G[i][j]-1;
				//	cout<<"inserted at row "<<i<<" and column "<<j<<endl;

			}
		}
	}
	//cout<<"done with for loop"<<endl;
	B[(B.size()-1)][count]=B[(B.size()-1)][count]-1;
	complex<double> b= complex<double>(value1,value2);
	U.push_back(b);
	//cout<<"done with inserting value in B matrices"<<endl;
}

void insertvaluesZ(vector<vector<double> > &G,vector<vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
{
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double NCplus=nodeK-1;
	double NCminus=nodeL-1;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if((i==Nplus && j==NCplus && Nplus>=0 && NCplus>=0)||(i==Nminus && j==NCminus && Nminus>=0 && NCminus>=0))
			{					
				C[i][j]=C[i][j]+value;
				//	cout<<"row :"<<i<<" column :"<<j<<endl;
			}
			else if((i==Nplus && j==NCminus && Nplus>=0 && NCminus>=0)||(i==Nminus && j==NCplus && Nminus>=0 && NCplus>=0))
			{	C[i][j]=C[i][j]-value;
				//	cout<<"row :"<<i<<" column :"<<j<<endl;

			}

		}
	}
}

void insertvaluesE(vector<vector<double> > &G,vector<vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
{
	int last=(G.size())-1;
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double NCplus=nodeK-1;
	double NCminus=nodeL-1;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if(i==Nplus && j==last && Nplus>=0)
				G[i][j]=G[i][j]+1;
			else if(i==Nminus && j==last && Nminus>=0)
				G[i][j]=G[i][j]-1;
			else if(i==last && j==Nplus && Nplus>=0)
				G[i][j]=G[i][j]+1;
			else if(i==last && j==Nminus && Nminus>=0)
				G[i][j]=G[i][j]-1;
			else if(i==last && j==NCplus && NCplus>=0)
				G[i][j]=G[i][j]-value;
			else if(i==last && j==NCminus && NCminus>=0)
				G[i][j]=G[i][j]+value;
		}
	}
}



void insertvaluesH(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > & L_current_position,double nodeI,double nodeJ,double nodeK,double nodeL,double value,double row,double column)
{
	int last=(G.size())-1;
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double section=row-1;
	double dependent_line=column-1;
	//cout<<"section :"<<section<<" line :"<<dependent_line<<endl;
	for(int i=0;i!=(G.size());++i){
		for(int j=0;j!=((G[i]).size());++j){
			if(i==Nplus && j==last && Nplus>=0)
				G[i][j]=G[i][j]+1;
			else if(i==Nminus && j==last && Nminus>=0)
				G[i][j]=G[i][j]-1;
			else if(i==last && j==Nplus &&Nplus>=0)
				G[i][j]=G[i][j]+1;
			else if(i==last && j==Nminus &&Nminus>=0)
				G[i][j]=G[i][j]-1;
			else if(i==last && j==L_current_position[section][dependent_line])
			{
				//cout<<"the current position of H is "<<j<<endl;
				C[i][j]=C[i][j]-value;
			}
		}
	}
}




void AinverseB(vector<vector<double> > &Ad,vector<vector<double> > &Bd,vector<vector<double> > &Cd)
{
	int sizeA=Ad.size();
	int Nd=sizeA;
	int colB=(Bd[0]).size();
	int NRHSd=colB;
	int LDAd =Nd;
	int LDBd=Nd;
	int IPIVd[Nd];
	int INFOd;
	double ad[LDAd*Nd],bd[LDBd*NRHSd];
	int kd=0;
	int ld=0;
	for(int i=0;i!=(Ad[0]).size();++i)
		for(int j=0;j!=Ad.size();++j)
			ad[kd++]=Ad[j][i];
	for(int i=0;i!=(Bd[0]).size();++i)
		for(int j=0;j!=Bd.size();++j)
			bd[ld++]=Bd[j][i];
	dgesv_( &Nd, &NRHSd, ad, &LDAd, IPIVd, bd, &LDBd, &INFOd );
	if( INFOd > 0 ) {
		cout<<"The diagonal element of the triangular factor of A,\n";
		cout<< "U("<<INFOd<<","<<INFOd<<" is jero, so that A is singular;\n";
		cout<<( "the solution could not be computed.\n" );
		//	return(1);
	}

	//	cout << "===============================IN AINVERSEB 1" << endl;
	vector<double>tempd;
	for(int i=0;i!=Nd;++i)
	{
		for(int j=0;j!=NRHSd;++j)
			tempd.push_back(bd[i+j*LDBd]);
		Cd.push_back(tempd);
		tempd.clear();
	}
	//	cout << "===============================IN AINVERSEB 2" << endl;
}

void QR_factorization(vector<vector<double> > &a,vector<vector<double> > &E)
{
	vector<double> u;
	for(int i=0;i!=(a[0]).size();++i){
		for(int j=0;j!=a.size();++j)
			u.push_back(a[j][i]);
		//		cout<<"printing"<< i+1 <<"u vector"<<endl;
		//		for(int r=0;r!=u.size();++r)
		//			cout<<u[r]<<" ";
		//	cout<<endl;

		if(i==0)
		{
			double sum=0;
			for(int k=0;k!=u.size();++k)
				sum+=u[k]*u[k];
			for(int k=0;k!=u.size();++k){
				double root=sqrt(sum);
				double value=u[k]/root;
				E.push_back(vector<double>(1,value));
			}
			u.clear();
		}
		else{
			double sum=0;
			for(int k=0;k!=i;++k)
			{
				double prod=0;
				for(int z=0;z!=u.size();++z)
					prod+=a[z][i]*E[z][k];
				for(int z=0;z!=u.size();++z)
					u[z]=u[z]-(prod*E[z][k]);
			}
			for(int k=0;k!=u.size();++k)
				sum+=u[k]*u[k];
			for(int k=0;k!=u.size();++k){
				double root=sqrt(sum);
				double value=u[k]/root;
				E[k].push_back(value);
			}
			u.clear();
		}
	}
}

void vectormultiplication(vector<vector<double> > &a,vector<vector<double> > &b,vector<vector<double> > &c)
{
	if(((a.at(0)).size())!=b.size()){
		cout<<"columns of 1st matrix and row of 2nd matrix does not match"<<endl;
		cout<<((a.at(0)).size())<<";;;;;;;;;;;;;;;;;;;"<<b.size()<<endl;

	}


	else{
		vector<double>temp;
		int row=a.size();
		int column=(b[0]).size();
		int inner=b.size();
		for(int i=0;i!=row;++i){
			for(int j=0;j!=column;++j){
				double sum=0;
				for(int k=0;k!=inner;++k){
					sum+=a[i][k]*b[k][j];
				}
				temp.push_back(sum);
			}
			c.push_back(temp);
			temp.clear();
		}
	}
}

void vector_transpose(vector<vector<double> > &a,vector<vector<double> > &b)
{
	vector<double>temp;
	for(int i=0;i!=(a[0]).size();++i)
	{
		for(int j=0;j!=a.size();++j)
		{
			temp.push_back(a[j][i]);
		}
		b.push_back(temp);
		temp.clear();
	}
}

void blockArnoldi(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<vector<double> > &X,int &q_value, vector<double> &range)
{
	vector<double> insert;
	vector<vector<double> > R,Q,V,X_temp,Z,X_transpose,H,XH;
	AinverseB(G,B,R);
	//	cout << "HERE------------------------------------------->>>>>>>>>>>>>>>>>" << endl;
	QR_factorization(R,Q);
	cout<<"the number of columns in one moment are :"<<(Q.at(0)).size()<<" and row are :"<<Q.size()<<endl;	
	int q,n,input1,input2;
	//	cout << "===============Q=" << Q.size() << endl;
	X=Q;
	//	cout << "===============X=" << X.size() << endl;
	//	cout<<"printing the X 0"<<endl;
	cout<<"done computing X0"<<endl;
	//	print(X);	
	Q.clear();
	R.clear();
	range.push_back((X[0]).size());
	cout<<"Enter the value of q"<<endl;
	//	cin>>q;
	q=q_value;
	cout<<"the value of q is "<<q<<endl;		
	if((q%((B[0]).size()))==0)
		n=q/((B[0]).size());
	else
		n=floor(q/((B[0]).size()))+1;
	cout<<"the value of n is :"<<n<<endl;	
	for(int kA=1;kA<=n;++kA)
	{
		//	cout<<"getting previous columsn"<<endl;
		input1=kA-1;
		if(input1==0){
			for(int l=0;l!=X.size();++l)
			{
				for(int m=0;m!=range[input1];++m)
					insert.push_back(X[l][m]);
				X_temp.push_back(insert);
				insert.clear();
			}
		}
		else{
			for(int l=0;l!=X.size();++l)
			{
				for(int m=range[input1-1];m!=range[input1];++m)
					insert.push_back(X[l][m]);
				X_temp.push_back(insert);
				insert.clear();
			}
		}
		//		cout<<"done getting previous column"<<endl;
		//		cout<<"performing vector multiplication"<<endl;		
		vectormultiplication(C,X_temp,V);
		X_temp.clear();
		//		cout<<"doing inverse"<<endl;
		AinverseB(G,V,Z);
		V.clear();
		for(int pass=1;pass<=2;++pass)		
		{
			for(int jA=1;jA<=kA;++jA)
			{
				input2=kA-jA;
				if(input2==0){
					for(int l=0;l!=X.size();++l)
					{
						for(int m=0;m!=range[input2];++m)
							insert.push_back(X[l][m]);
						X_temp.push_back(insert);
						insert.clear();
					}
				}
				else{
					for(int l=0;l!=X.size();++l)
					{
						for(int m=range[input2-1];m!=range[input2];++m)
							insert.push_back(X[l][m]);
						X_temp.push_back(insert);
						insert.clear();
					}
				}
				vector_transpose(X_temp,X_transpose);
				vectormultiplication(X_transpose,Z,H);
				vectormultiplication(X_temp,H,XH);
				for(int i_A=0;i_A!=Z.size();++i_A)
					for(int j_A=0;j_A!=(Z[0]).size();++j_A)
						Z[i_A][j_A]=Z[i_A][j_A]-XH[i_A][j_A];
				X_temp.clear();
				X_transpose.clear();
				XH.clear();
				H.clear();
			}
		}
		//		cout<<"doing QrR"<<endl;		
		QR_factorization(Z,Q);
		for(int iAB=0;iAB!=Q.size();++iAB)
			for(int jAB=0;jAB!=(Q[0]).size();++jAB)
				(X[iAB]).push_back(Q[iAB][jAB]);
		range.push_back((X[0]).size());
		cout<<"printing the X"<<kA<<endl;
		//		print(Q);		
		Q.clear();
		Z.clear();
	}
}

void print_eigenvalues( char* desc, int n,double* wr,double* wi ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) {
		//      if( wi[j] == (double)0.0 ) {
		//       printf( " %6.2f", wr[j] );
		//  } else {
		printf( " (%f,%f)", wr[j], wi[j] );
		//  }
	}
	printf( "\n" );
}

int block_inverse(vector<vector<double> > A_real,vector<vector<double> > A_imag,vector<vector<double> > B,vector<vector<double> > &vec_soln,double expansion_frequency)
{
	int N=A_real.size();
	complex<double> a[N*N];
	int k1=0;
	//cout<<"combining:-----------------------------------------------------"<<endl;
	for(int i=0;i<N;++i){
		for(int j=0;j<N;++j)
		{

			a[k1++]=complex<double>(A_real[j][i],2*pi*expansion_frequency*A_imag[j][i]);
			//cout<<complex<double>(A_real[j][i],2*pi*expansion_frequency*A_imag[j][i])<<" ";
		}
		//cout<<endl;
	}
	//cout<<"B ===================================================="<<endl;
	complex<double> b[N*(B.at(0)).size()];
	int k2=0;
	for(int i=0;i<(B.at(0)).size();++i)
	{
		for(int j=0;j<N;++j)
		{
			b[k2++]=complex<double>(B[j][i],0);
			//cout<<complex<double>(B[j][i],0)<<" ";	
		}
		//cout<<endl;
	}

	int n=N;
	int nrhs=(B.at(0)).size();
	int lda=n,ldb=n,info;
	int ipiv[n];
	zgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );

	complex<double> sol[n][nrhs];
	//cout<<"array==================================="<<endl;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<nrhs;++j)
		{
			sol[i][j]=b[i+j*ldb];
			//cout<<sol[i][j]<<" ";		
		}
		//cout<<endl;
	}


	vector<double> temp;
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<nrhs;++j)
		{
			temp.push_back(real(sol[i][j]));
			temp.push_back(imag(sol[i][j]));
		}
		vec_soln.push_back(temp);
		temp.clear();
	}

	return 0;
}

int complex_vector_multiplication(vector<vector<complex<double> > > a,vector<vector<complex<double> > > b,vector<vector<complex<double> > > &c)
{
	if((a.at(0)).size()!=b.size()){
		cout<<"columns of 1st matrix and row of 2nd matrix does not match"<<endl;
		cout<<(a.at(0)).size()<<"''''''''''''''''''''''"<<b.size()<<endl;

	}

	else{
		vector<complex<double> > temp;
		int row=a.size();
		int column=(b.at(0)).size();
		int inner=b.size();
		for(int i=0;i<row;++i)
		{
			for(int j=0;j<column;++j){
				complex<double> sum=0;
				for(int k=0;k<inner;++k)
				{
					sum+=a[i][k]*b[k][j];
				}
				temp.push_back(sum);
			}
			c.push_back(temp);
			temp.clear();
		}
	}
	return 2;
}




int print_complex_vector(vector<vector<complex<double> > > a)
{
	for(int i=0;i<a.size();++i)
	{
		for(int j=0;j<(a.at(0)).size();++j)
		{
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
	return 3;
}



int multi_Arnoldi(vector<vector<double> > G,vector<vector<double> > C,vector<vector<double> > B,vector<vector<double> > &X,int q_value, double expansion_frequency)
{
	vector<double> insert,temp;
	vector<vector<double> > V_sep,GinverseC,GinverseB,Z,X_sep,X_transpose,complex_sep,H,Q,XH;
	vector<complex<double> >c1,c2;
	vector<vector<complex<double> > > GinverseC_complex,V,Z_complex;
	//doing Ginverc
	//cout<<X.size()<<"-------------------------------"<<(X.at(0)).size()<<endl;

	//cout<<"performing Ginversc-----------------------------------------------"<<endl;
	block_inverse(G,C,C,GinverseC,expansion_frequency);
	//converting to complex
	//cout<<"size of GinversC "<<(GinverseC.at(0)).size()<<"---------------------------"<<GinverseC.size()<<"----------------------------------------"<<endl;
	//cout<<"-------------------------------------------------------------"<<endl;
	for(int i=0;i<GinverseC.size();++i)
	{
		for(int j=0;j<((GinverseC.at(0)).size())/2;++j)
		{
			c1.push_back(complex<double> (GinverseC[i][2*j],GinverseC[i][(2*j)+1]));
		}
		GinverseC_complex.push_back(c1);
		c1.clear();
	}
	//cout<<"size of GinversCcomlex "<<(GinverseC_complex.at(0)).size()<<"------------------------"<<GinverseC_complex.size()<<"---------------------------------------------"<<endl;
	//Ginverse B
	block_inverse(G,C,B,GinverseB,expansion_frequency);
	cout<<"the number of columns in moments( GinversB) are : "<<(GinverseB.at(0)).size()<<" and the number of rows are :-"<<GinverseB.size()<<"-------------------"<<endl;
	for(int i=0;i<(GinverseB.at(0)).size();++i)
	{//#getting the each columns of GinverseB
		for(int j=0;j<GinverseB.size();++j)
		{
			insert.push_back(GinverseB[j][i]);
			Z.push_back(insert);
			insert.clear();
		}
		//cout<<i+1<<"---------------"<<Z.size()<<"-----------------------------"<<(Z.at(0)).size()<<endl;	
		for(int pass=1;pass<=2;++pass)
		{
			//cout<<pass<<"::::::::::::::::::::::::::::::::::::::::::"<<endl;
			//cout<<"X...."<<X.size()<<"-------------------------------"<<(X.at(0)).size()<<endl;
			for(int j=(((X.at(0)).size())-1);j>=0;--j)
			{
				for(int k=0;k<X.size();++k)
				{
					//cout<<k<<" ";				
					insert.push_back(X[k][j]);
					X_sep.push_back(insert);
					insert.clear();
				}
				//cout<<j<<"---"<<X_sep.size()<<"--------------------------------------------------"<<(X_sep.at(0)).size()<<endl;
				vector_transpose(X_sep,X_transpose);
				vectormultiplication(X_transpose,Z,H);
				vectormultiplication(X_sep,H,XH);
				for(int i_A=0;i_A<Z.size();++i_A)
				{
					for(int j_A=0;j_A<(Z.at(0)).size();++j_A)
					{
						Z[i_A][j_A]=Z[i_A][j_A]-XH[i_A][j_A];
					}
				}
				X_transpose.clear();
				XH.clear();
				H.clear();
				X_sep.clear();
			}
		}
		QR_factorization(Z,Q);
		Z.clear();
		for(int iAB=0;iAB<Q.size();++iAB)
		{
			for(int jAB=0;jAB<(Q.at(0)).size();++jAB)
			{
				(X[iAB]).push_back(Q[iAB][jAB]);
			}
		}

		//cout<<"X"<<(X.at(0)).size()<<endl;	
		if(complex_sep.empty())
		{
			complex_sep=Q;
			Q.clear();
		}
		else{
			for(int iAB=0;iAB<Q.size();++iAB)
			{
				for(int jAB=0;jAB<(Q.at(0)).size();++jAB)
				{
					(complex_sep[iAB]).push_back(Q[iAB][jAB]);
				}
			}
			Q.clear();
		}
	}

	cout<<"printing the X0"<<endl;
	//	cout<<"size of X "<<(X.at(0)).size()<<"---------------------------"<<X.size()<<"----------------------------------------"<<endl;
	//	cout<<"size of complex_sep "<<(complex_sep.at(0)).size()<<"---------------------------"<<complex_sep.size()<<"----------------------------------------"<<endl;	
	//cout<<"Enter the value of q"<<endl;
	cout<<"the value of q is "<<q_value<<endl;
	int n;
	if((q_value%((B[0]).size()))==0)
		n=q_value/((B[0]).size());
	else
		n=floor(q_value/((B[0]).size()))+1;
	cout<<"the value of n is :"<<n<<endl;
	//staring the GinversC part from here
	for(int k=1;k<=n;++k)
	{
		for(int i=0;i<complex_sep.size();++i)
		{
			for(int j=0;j<((complex_sep.at(0)).size()/2);++j)
			{
				c2.push_back(complex<double>(complex_sep[i][2*j],complex_sep[i][(2*j)+1]));
			}
			V.push_back(c2);
			c2.clear();
		}
		complex_sep.clear();
		//cout<<"size of V :"<<V.size()<<"----------------------------------"<<(V.at(1)).size()<<endl;
		complex_vector_multiplication(GinverseC_complex,V,Z_complex);
		//cout<<"size of Z_complex :"<<Z_complex.size()<<"----------------------------------"<<(Z_complex.at(0)).size()<<endl;

		V.clear();
		for(int i=0;i<Z_complex.size();++i)
		{
			for(int j=0;j<(Z_complex.at(0)).size();++j)
			{
				temp.push_back((Z_complex[i][j]).real());
				temp.push_back((Z_complex[i][j]).imag());
			}
			V_sep.push_back(temp);
			temp.clear();
		}
		//cout<<"size of V_sep :"<<V_sep.size()<<"----------------------------------"<<(V_sep.at(0)).size()<<endl;

		Z_complex.clear();
		for(int i=0;i<(V_sep.at(0)).size();++i)
		{//#getting the each columns of GinverseB
			for(int j=0;j<V_sep.size();++j)
			{
				insert.push_back(V_sep[j][i]);
				Z.push_back(insert);
				insert.clear();
			}
			for(int pass=1;pass<=2;++pass)
			{
				for(int j=((X.at(0)).size()-1);j>=0;--j)
				{
					for(int k=0;k<X.size();++k)
					{
						insert.push_back(X[k][j]);
						X_sep.push_back(insert);
						insert.clear();
					}
					vector_transpose(X_sep,X_transpose);
					vectormultiplication(X_transpose,Z,H);
					vectormultiplication(X_sep,H,XH);
					for(int i_A=0;i_A<Z.size();++i_A)
					{
						for(int j_A=0;j_A<(Z.at(0)).size();++j_A)
						{
							Z[i_A][j_A]=Z[i_A][j_A]-XH[i_A][j_A];
						}
					}
					X_transpose.clear();
					XH.clear();
					H.clear();
					X_sep.clear();
				}
			}
			QR_factorization(Z,Q);
			Z.clear();
			for(int iAB=0;iAB<Q.size();++iAB)
			{
				for(int jAB=0;jAB<(Q.at(0)).size();++jAB)
				{
					(X[iAB]).push_back(Q[iAB][jAB]);
				}
			}


			if(complex_sep.empty())
			{
				complex_sep=Q;
				Q.clear();
			}
			else{
				for(int iAB=0;iAB<Q.size();++iAB)
				{
					for(int jAB=0;jAB<(Q.at(0)).size();++jAB)
					{
						(complex_sep[iAB]).push_back(Q[iAB][jAB]);
					}
				}
				Q.clear();
			}
		}

		V_sep.clear();	

		cout<<"printing the X"<<k<<endl;
	}
	return 4;
}


int SVD(vector<vector<double> > &X_svd,vector<double> &sigma,vector<vector<double> >X_big)
{
	int m=X_big.size();
	int n=(X_big.at(0)).size();
	int lda =m,ldu=m,ldvt=n,info,lwork;
	double wkopt;
	double* work;
	int sn;

	if(m>n) sn=n;
	else  sn=m;

	double s[sn],u[ldu*m],vt[ldvt*n],a[lda*n];

	int count=0;
	for(int i=0;i<(X_big.at(0)).size();++i)
	{
		for(int j=0;j<X_big.size();++j)
		{
			a[count++]=X_big[j][i];
		}
	}

	lwork = -1;
	dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	/* Compute SVD */
	dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,&info );
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
	}
	vector<double> temp;

	for(int i=0;i<sn;++i)
	{
		sigma.push_back(s[i]);
	}

	for(int i=0;i<m;++i)
	{
		for(int j=0;j<sn;++j)
		{
			temp.push_back(u[i+j*ldu]);
		}
		X_svd.push_back(temp);
		temp.clear();

	}

	return 31;
}


int retrun_X_MATRIX(vector<vector<double> > &X_MATRIX,vector<vector<double> > X_svd,int count)
{
	vector<double> temp_svd;
	for(int i=0;i<X_svd.size();++i)
	{
		for(int j=0;j<count;++j)
		{
			temp_svd.push_back(X_svd[i][j]);
		}
		X_MATRIX.push_back(temp_svd);
		temp_svd.clear();
	}
	return 32;
}


int original_plus_MOR_sending_X_MATRIX(vector<vector<double> > X_MATRIX,double frequency_max1,double frequency_min1,double points1,string filename,int mynode,string s_mynode)
{	
	double interval;
	interval=(frequency_max1-frequency_min1)/(points1-1);
	string in_1,in_2,in_3;
	in_1=filename;
	cout<<"ENTER THE MAIN NETLIST FILE : ";
	cout<<"1) MAIN NETLIST FILE : "<<in_1;
	vector<complex<double> > u_MATRIX;
	vector<vector<double> >G_MATRIX;
	vector<vector<double> >C_MATRIX;
	vector<vector<double> >B_MATRIX;
	G_MATRIX.clear();
	C_MATRIX.clear();
	B_MATRIX.clear();


	//reading the input and output node from files created from pul.c code        
	vector<double> input_node,output_node;
	read_nodes("input_terminal.txt", input_node);
	read_nodes("output_terminal.txt", output_node);
	//printing the nodes
	for(int i=0;i<input_node.size();++i)
		cout<<input_node[i]<<" ";
	cout<<endl;
	for(int i=0;i<output_node.size();++i)
		cout<<output_node[i]<<" ";
	cout<<endl;

	populate_matrices(in_1, G_MATRIX, C_MATRIX, B_MATRIX,u_MATRIX);
	//	cout<<"G MATRIX"<<endl;	
	//	print(G_MATRIX);
	//	cout<<endl;
	//	cout<<"C MATRIX"<<endl;
	//	print(C_MATRIX);
	cout<<"THE NUMBER OF VARIABLE ARE :"<<G_MATRIX.size()<<endl;
	//	cout<<endl<<"B MATRIX"<<endl;
	//	print(B_MATRIX);
	//	cout<<endl<<"u complex vector"<<endl;
	//	printcomplexvector(u_MATRIX);
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(int i=0;i<(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	//printing tht ucomplexaray
	//	cout<<endl<<"ucomplexarray"<<endl;
	//	for(int i=0;i!=(u_MATRIX).size();++i)
	//		cout<<uarraymatrix[i][0]<<endl;
	//converting the B VECTOR MATRIX to B ARRAY MATRIX
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(int i=0;i<(B_MATRIX).size();++i)
		for(int j=0;j!=(u_MATRIX).size();++j)
			B_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_MATRIX[i][j],0);
	// print the the B ARRAY MATRIX
	//	cout<<endl<<"B_COMPLEXARRAYMATRIX"<<endl;
	//	for(int i=0;i!=(B_MATRIX).size();++i){
	//		for(int j=0;j!=(u_MATRIX).size();++j)
	//			cout<<B_COMPLEXARRAYMATRIX[i][j]<<"  ";
	//		cout<<endl;}

	//multiplying the B_COMPLEXARRAY AND uarraymatrix
	complex<double> Bucomplexmultiple[(B_MATRIX).size()][1];
	for(int i=0;i!=(B_MATRIX).size();++i)
		for(int j=0;j!=(u_MATRIX).size();++j)
			Bucomplexmultiple[i][0]+=B_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];

	/*	//comment for original solving starts here	
		char answe;
		cout<<"DO YOU WANT TO SOLVE FOR THE ORIGINAL VARIABLES(y/n) :";
		cin>>answe;
		if(answe=='y')
		{
		double t1=clock();

		fstream result11("SVD_MATLAB_result_at_output_node_at_multiple_frequencies.txt",fstream::out);
		if(result11.is_open() )
		{
		double frequency;
		int timer1=0;
		while(timer1<points1)
		{
		cout<<"original"<<timer1+1<<endl;	
		frequency=frequency_min1+(interval*timer1);
		result11<<frequency<<'\t';

	//combining G_MATRIX and C_MATRIX into GplussC_MATRIX
	complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
	for(int i=0;i<(G_MATRIX).size();++i)
	{
	for(int j=0;j<(G_MATRIX).size();++j)
	{
	GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
	}
	}
	//copying elements of GplusC_MATRIX in a one dimension array a
	complex<double> a[(G_MATRIX.size())*(G_MATRIX.size())];
	int k=0;
	for(int i=0;i!=(G_MATRIX).size();++i)
	for(int j=0;j!=(G_MATRIX).size();++j)
	a[k++]=GplussC_MATRIX[j][i];

	//printing the one dimension a matrix
	//		cout<<endl<<" one dimension 'a' matrix"<<endl;
	//		for(int i=0;i!=((G_MATRIX.size())*(G_MATRIX.size()));++i)
	//			cout<<a[i]<<" ";
	//		cout<<endl;
	//
	//copying elements of Bucomplexmultiple into b one dimension matrix

	complex<double> b[B_MATRIX.size()];
	for(int i=0;i!=(B_MATRIX).size();++i)
	b[i]=Bucomplexmultiple[i][0];

	// printing the b one dimension matrix
	//		cout<<"one dimension b matrix"<<endl;
	//		for(int i=0;i!=(B_MATRIX).size();++i)
	//			cout<<b[i];
	//		cout<<endl;
	//computing ax=b using zgesv routine
	//			cout<<"computing Ax=b of orginal solution"<<endl;
	int n=G_MATRIX.size();
	int nrhs=1;
	int lda=n;
	int ldb=n;
	int info;
	int ipiv[n];
	//	printf( " ZGESV Program Results\n" );
	zgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	// Check for the exact singularity 
	if( info > 0 ) {
	printf( "The diagonal element of the triangular factor of A,\n" );
	printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
	printf( "the solution could not be computed.\n" );
	return( 1 );
	}
	//	cout<<"done computing Ax=B of orginal solution"<<endl;
	complex<double> original_solution[n][1];
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<nrhs;++j)
		{
			original_solution[i][j]=b[i+j*ldb];
		}
	}
	for(int i=0;i<G_MATRIX.size();++i)
	{
		for(int k=0;k<output_node.size();++k)
		{
			if(i==(output_node[k]-1)){

				if((original_solution[i][0]).imag()>=0)
					result11<<(original_solution[i][0]).real()<<'+'<<(original_solution[i][0]).imag()<<'i'<<'\t';
				else
					result11<<(original_solution[i][0]).real()<<(original_solution[i][0]).imag()<<'i'<<'\t';

			}

		}
	}
	result11<<';'<<endl;

	++timer1;
}//for loop of frequency ends here

}
else cout<<"unable to open file"<<endl;
double t2=clock();
cout<<" TIME taken to compute the original model is : "<<(t2-t1)/double(CLOCKS_PER_SEC) << "seconds "<<endl;
cout<<"done computing Ax=b of original solution for different frequencies and writing for MATLAB READING in SVD_MATLAB_result_at_output_node_at_multiple_frequencies.txt"<<endl;
}
*/	//comment for original model solving ends here
/*
   cout<<"checking the matrices generated is orthonormal or not"<<endl;
   vector<vector<double> > X_MATRIX_TRANSPOSE,IDENTITY;
   vector_transpose(X_MATRIX,X_MATRIX_TRANSPOSE);
   vectormultiplication(X_MATRIX_TRANSPOSE,X_MATRIX,IDENTITY);
   cout<<"checking the orthonormality"<<endl;
   fstream wr_o("SVD_orthonormality_check.txt",fstream::out);

   if(wr_o.is_open())
   {
   for(int i=0;i<IDENTITY.size();++i)
   {
   for(int j=0;j<(IDENTITY.at(0)).size();++j)
   {
   wr_o<<IDENTITY[i][j]<<'\t';
   }
   wr_o<<endl;
   }
   wr_o.close();
   }
   else cout<<"unable to open file"<<endl;
   cout<<" wrting the results of orthonormality check in SVD_orthonormality_check.txt"<<endl;
 */

//generating the reduced matrices
cout<<"generating the reduced Matrices of orginal matrices"<<endl;
vector<vector<double> > X_transpose,multiple,G_cap,C_cap,B_cap;
vector_transpose(X_MATRIX,X_transpose);
vectormultiplication(X_transpose,X_MATRIX,multiple);
//		cout<<"printing the multiplication of Z matrices and its transpose"<<endl;
//		print(multiple);
multiple.clear();
vectormultiplication(X_transpose,G_MATRIX,multiple);
vectormultiplication(multiple,X_MATRIX,G_cap);
multiple.clear();
vectormultiplication(X_transpose,C_MATRIX,multiple);
vectormultiplication(multiple,X_MATRIX,C_cap);
multiple.clear();
vectormultiplication(X_transpose,B_MATRIX,B_cap);
//		cout<<"printing the Gcap"<<endl;
//		print(G_cap);
//		cout<<"printing the Ccap"<<endl;
//		print(C_cap);
//		cout<<"printing the Bcap"<<endl;
//		print(B_cap);

cout<<"done generating the reduced Matrices of orginal matrices"<<endl;	
X_transpose.clear();
vector_transpose(X_MATRIX,X_transpose);
cout<<"the number of X reduced are :"<<G_cap.size()<<endl; 
//double t3=clock();
string X_str_col;
ostringstream ss1;
ss1<<(X_MATRIX.at(0)).size();
X_str_col=ss1.str();
string out_SVD_file="./result/SVD_MOR_s"+s_mynode+"("+X_str_col+").txt";
cout<<"MYNODE :"<<mynode<<" writing in "<<out_SVD_file<<endl;
double t1=MPI_Wtime();	
fstream result22(out_SVD_file.c_str(),fstream::out);
if(result22.is_open() )
{
	double frequency;
	int timer2=0;	
	while(timer2<points1)
	{
		cout<<"MYNODE :"<<mynode<<" computing SVD MOR at "<<timer2+1<<" frequency point"<<endl;			
		frequency=frequency_min1+(interval*timer2);
		result22<<frequency<<'\t';
		// Combining G_cap and C_cap into G_capplussC_cap
		complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
		for(int i=0;i<(G_cap).size();++i){
			for(int j=0;j<(G_cap).size();++j)
			{
				G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
			}
		}

		//copying elements of GcapplusCcap in a one dimension array acap
		complex<double> a_cap[(G_cap.size())*(G_cap.size())];
		int k_cap=0;
		for(int i=0;i<(G_cap).size();++i)
			for(int j=0;j<(G_cap).size();++j)
				a_cap[k_cap++]=G_capplussC_cap[j][i];

		//converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX
		complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
		for(int i=0;i<(B_cap).size();++i)
			for(int j=0;j<(u_MATRIX).size();++j)
				B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);

		//multiplying the B_COMPLEXARRAY AND uarraymatrix
		complex<double> B_capucomplexmultiple[(B_cap).size()][1];
		for(int i=0;i<(B_cap).size();++i)
			for(int j=0;j<(u_MATRIX).size();++j)
				B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];


		//copying elements of B_capucomplexmultiple into b_cap one dimension matrix
		complex<double> b_cap[(B_cap).size()];
		for(int i=0;i<(B_cap).size();++i)
			b_cap[i]=B_capucomplexmultiple[i][0];

		//computing a_capx=b_cap using cgesv routine

		int n_cap=G_cap.size();
		int nrhs_cap=1;

		int lda_cap=n_cap;
		int ldb_cap=n_cap;
		int info_cap;
		int ipiv_cap[n_cap];
		zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
		if( info_cap > 0 ) {
			printf( "The diagonal element of the triangular factor of A,\n" );
			printf( "U(%i,%i) is Jero, so that A is singular;\n", info_cap, info_cap );
			printf( "the solution could not be computed.\n" );
			return( 1 );
		}

		//converting X_matrix into X complex arary matrices

		complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
		for(int i=0;i<X_MATRIX.size();++i)
			for(int j=0;j<(X_MATRIX[0]).size();++j)
				X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

		//converting b_cap from inverse into x_cap 2 dimension array
		complex<double> x_cap[G_cap.size()][1];

		for(int i=0;i<n_cap;++i)
			for(int j=0;j<nrhs_cap;++j)
				x_cap[i][j]=b_cap[i+j*ldb_cap];

		complex<double>MORsolutioncompare[X_MATRIX.size()][1];
		for(int i=0;i<X_MATRIX.size();++i)
			for(int j=0;j<(X_MATRIX[0]).size();++j)
				MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];
		for(int i=0;i<G_MATRIX.size();++i)
		{
			for(int k=0;k<output_node.size();++k)
			{
				if(i==(output_node[k]-1)){
					if((MORsolutioncompare[i][0]).imag()>=0)              
						result22<<(MORsolutioncompare[i][0]).real()<<'+'<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';
					else
						result22<<(MORsolutioncompare[i][0]).real()<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';

				}
			}
		}
		result22<<';'<<endl;
		++timer2;			
	}
	result22.close();

}
else cout<<"unable top open file"<<endl;
//double t4=clock();
//      cout<<" TIME taken to compute the MOR from SVD MATRIX solution is :"<<(t4-t3)/double(CLOCKS_PER_SEC) << "seconds"<<endl; 
cout<<"done computing Ax=b of MOR solution for different frequencies and writing for MATLAB READING in "<<out_SVD_file<<endl;		
double t2=MPI_Wtime();
cout<<"TIME taken to solve reduced model with after SVD :"<<t2-t1<<endl;
return 33;
}




int read_nodes(string filename, vector<double>& node_list)
{
	double temp_node;
	fstream fp(filename.c_str(),fstream::in);
	while(fp>>temp_node)
		node_list.push_back(temp_node);
	fp.close();

	return 0;
}


int populate_matrices(string filename, vector<vector<double> >& G_MATRIX, vector<vector<double> >& C_MATRIX, vector<vector<double> >& B_MATRIX, vector<complex<double> >& u_MATRIX)
{

	string line,s,com1,com2;
	vector<double>values;
	vector<string>name;	
	vector<string>read; //WATCH FOR SPACE IN INPUT    
	int column=0;

	ifstream myfile(filename.c_str());	
	cout<<"Generating the individual Matrices"<<endl;
	while(getline(myfile,line))/* && !line.empty())*/{

		for(int i=0;i<(line.size());++i){
			if(isspace(line[i])&& !s.empty()){
				read.push_back(s);
				s.clear();
			}
			else if(!isspace(line[i])){
				s+=line[i];
			}
		}
		read.push_back(s);
		s.clear();//till here stored the strings from line in vector called read.
		//now reading that vector and taking out which are numbers and which are string.		
		for(int noofelements=0;noofelements!=read.size();++noofelements)
		{

			string &element=read[noofelements];
			//determining the which is string and which is value
			int count=0;
			for(int i=0;i!=element.size();++i)
			{
				if(isdigit(element[i])||ispunct(element[i])||(isdigit(element[i-1])&&element[i]=='e')||(isdigit(element[i-1])&&element[i]=='E'))
					++count;
			}
			if(count==element.size())
			{
				//taking each values and determining where is comma and word e or E
				int comma_place;
				int comma_counter=0;
				int e_counter=0;
				vector<int>e_place;
				for(int i=0;i!=element.size();++i)
				{
					if(element[i]==','){
						comma_place=i;
						++comma_counter;
					}
					else if(element[i]=='e'||element[i]=='E')
					{
						++e_counter;
						e_place.push_back(i);
					}
				}
				if(comma_counter==0&&e_counter==0)
				{
					double f;
					istringstream(element)>>f;
					values.push_back(f);
					values.push_back(0);
				}
				else if(comma_counter==1&&e_counter==0)
				{
					for(int z=0;z<comma_place;++z)
						s+=element[z];                                
					double f;
					istringstream(s)>>f;
					values.push_back(f);
					s.clear();
					for(int z=(comma_place+1);z!=element.size();++z)
						s+=element[z];
					double g;
					istringstream(s)>>g;
					values.push_back(g);
					s.clear();
				}

				else if(comma_counter==0&& e_counter==1)
				{
					double number,power,scientific;
					for(int z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number;
					s.clear();
					for(int z=(e_place[0]+1);z!=element.size();++z)
						s+=element[z];
					istringstream(s)>>power;
					s.clear();
					scientific=number*pow(10,power);
					values.push_back(scientific);
					values.push_back(0);
				}

				else if(comma_counter==1&& e_counter==2)
				{
					double number1,power1,scientific1,number2,power2,scientific2;
					for(int z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number1;
					s.clear();
					for(int z=(e_place[0]+1);z!=comma_place;++z)
						s+=element[z];
					istringstream(s)>>power1;
					s.clear();
					scientific1=number1*pow(10,power1);
					values.push_back(scientific1);
					for(int z=comma_place+1;z!=e_place[1];++z)
						s+=element[z];
					istringstream(s)>>number2;
					s.clear();
					for(int z=(e_place[1]+1);z!=element.size();++z)
						s+=element[z];
					istringstream(s)>>power2;
					s.clear();
					scientific2=number2*pow(10,power2);
					values.push_back(scientific2);



				}
				else if(comma_counter==1&& e_counter==1)
				{
					if(comma_place<e_place[0])
					{
						double f,number,power,scientific;
						for(int z=0;z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>f;
						s.clear();
						values.push_back(f);
						for(int z=comma_place+1;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(int z=(e_place[0]+1);z!=element.size();++z)
							s+=element[z];
						istringstream(s)>>power;
						scientific=number*pow(10,power);
						values.push_back(scientific);
					}
					else if(comma_place>e_place[0])
					{
						double f,number,power,scientific;
						for(int z=0;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(int z=(e_place[0]+1);z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>power;
						s.clear();						
						scientific=number*pow(10,power);
						values.push_back(scientific);
						for(int z=comma_place+1;z!=element.size();++z)
							s+=element[z];
						istringstream(s)>>f;
						s.clear();
						values.push_back(f);
					}

				}
			}
			else
				name.push_back(element);
		}
		addzeroes(G_MATRIX,C_MATRIX,B_MATRIX,name[0][0],values[0],values[2],column);
		if(name[0][0]=='G'|| name[0][0]=='g')	
		{    
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}	
		else if(name[0][0]=='R' || name[0][0]=='r')
		{
			double newvalue=(1/values[4]);
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],newvalue);
		}
		else if(name[0][0]=='C' || name[0][0]=='c')
		{
			insertvaluesC(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}		
		else if(name[0][0]=='L' || name[0][0]=='l')
		{
			insertvaluesL(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}
		else if(name[0][0]=='J' || name[0][0]=='j')
		{
			insertvaluesJ(B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		}
		else if(name[0][0]=='V' || name[0][0]=='v')
		{
			insertvaluesV(G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		}		
		else if(name[0][0]=='Z' || name[0][0]=='z')
		{
			insertvaluesZ(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		}
		else if(name[0][0]=='E' || name[0][0]=='e')
		{
			insertvaluesE(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		}
		values.clear();
		name.clear();
		read.clear(); // to clear the read vector*/

	}
	cout<<"Done generating individual matrices"<<endl;

	return 27;
}





