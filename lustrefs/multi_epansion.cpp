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
using namespace std;
int print(vector<vector<double> >&);
int printcomplexvector(vector<complex<double> >&);
int addzeroes(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,char, double ,double ,int &);
int insertvaluesRG(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesC(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesL(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
int insertvaluesJ(vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
int insertvaluesV(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,int &);
int insertvaluesZ(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
int insertvaluesE(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
int insertvaluesH(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double,double,double);

#define pi 3.14
#define frequency_max 10e9
#define frequency_min 1e3
#define points 1000 
extern "C" int zgesv_( int *n , int *nrhs , complex<double> *a , int *lda , int *ipiv ,complex<double> *b , int *ldb , int *info  );
int print_matrix( char* desc, int m, int n, complex<double> *a, int lda );
int print_int_vector( char* desc, int n, int* a );
extern "C" int dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
int AinverseB(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
int QR_factorization(vector<vector<double> > &,vector<vector<double> > &);
int vectormultiplication(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
int vector_transpose(vector<vector<double> > &,vector<vector<double> > &);
int blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int & ,vector<double> &);
extern  "C" int dgeev_( char* jobvl, char* jobvr, int* n, double* a,int* lda,  double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info );
int multi_blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,int &,vector<double> &);
int print_eigenvalues( char* desc, int n, double* wr, double* wi );
int block_inverse(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,double );
int complex_vector_multiplication(vector<vector<complex<double> > > ,vector<vector<complex<double> > > ,vector<vector<complex<double> > > &);
int print_complex_vector(vector<vector<complex<double> > >);
int multi_Arnoldi(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &,int, double );




int main(int argc, char *argv[])
{       
	double interval;
	interval=(frequency_max-frequency_min)/(points-1);



	cout<<" DID U CHANGE THE FREQUENCY"<<endl/*<<"REMEMBER ARRAY IS A BITCH"<<endl*/<<"First enter RGC irrespective of the order,then enter L,then independent sources and then dependent sources"<<endl<<"R or r for resistance"<<endl<<"G or g for admittance"<<endl<<"C or c for capacitance"<<endl<<"L or l for inductance"<<endl<<"J or j for current source"<<endl<<"E or e for VCVS"<<endl<<"Z or z for VCCS"<<endl<<"H or h for CCVS"<<endl<<"V or v for independent voltage source"<<endl<<"PRESS ENTER TWICE AFTER DONE ENTERING THE CIRCUIT INFORAMTION IN THE FORM OF SPICE INPUT"<<endl;
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
		cout<<input_node[i]<<" ";
	cout<<endl;
	for(int i=0;i!=output_node.size();++i)
		cout<<output_node[i]<<" ";
	cout<<endl;

	ifstream myfile("stroud_1.txt");	
	cout<<"Generating the individual Matrices"<<endl;
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
			//	cout<<"the count is "<<count<<endl;			
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
		//cout<<name[0][0]<<endl;
		//	cout<<"the numbers in values are"<< endl;
		//		  for(int i=0;i!= values.size();++i)
		//		  cout<<values[i]<<endl;

		//	new_word=name[0][0];
		//cout<<"new_word is"<<new_word<<endl;
		addzeroes(G_MATRIX,C_MATRIX,B_MATRIX,name[0][0],values[0],values[2],column);
		/*	if((new_word=='H'||new_word=='h'||new_word=='v'||new_word=='V'||new_word=='Z'||new_word=='z') && (old_word=='L'||old_word=='l'))
			{
		//cout<<"HERE I AM"<<endl;
		L_columns.push_back(L_temp);
		if(name[0][0]=='V' || name[0][0]=='v')
		insertvaluesV(G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		else if(name[0][0]=='Z' || name[0][0]=='z')
		insertvaluesZ(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='H' || name[0][0]=='h')
		insertvaluesH(G_MATRIX,C_MATRIX,L_columns,values[0],values[2],values[4],values[6],values[8],values[10],values[12]);

		}*/
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
			/*	new_section=values[6];
				if(new_section==previous_section)
				L_temp.push_back(G_MATRIX.size());
				else
				{
				previous_section=new_section;
				L_columns.push_back(L_temp);
				L_temp.clear();
				L_temp.push_back(G_MATRIX.size());
				}
			 */
			insertvaluesL(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
		}
		else if(name[0][0]=='J' || name[0][0]=='j')
			insertvaluesJ(B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		else if(name[0][0]=='V' || name[0][0]=='v')
		{//	cout<<"inserting values of V here"<<endl;
			insertvaluesV(G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,values[0],values[2],values[4],values[5],column);
		}		
		else if(name[0][0]=='Z' || name[0][0]=='z')
			insertvaluesZ(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='E' || name[0][0]=='e')
			insertvaluesE(G_MATRIX,C_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='H' || name[0][0]=='h')
			insertvaluesH(G_MATRIX,C_MATRIX,L_columns,values[0],values[2],values[4],values[6],values[8],values[10],values[12]);	


		/*For debugging the complex values are stores in values
		  for(int r=0; r!=values.size();++r)
		  cout<<values[r]<<" ";
		  cout<<endl;*/
		values.clear();
		name.clear();
		read.clear(); // to clear the read vector*/

		//		old_word=new_word;
		//cout<<"old_word is "<<old_word<<endl;
	}
	/*	cout<<"columns of L of respective of sections are: "<<endl;
		cout<<"size of L_columns matrices is :"<<L_columns.size()<<endl;	
		for(int i=0;i!=L_columns.size();++i)
		{
		cout<<"the section is: "<<i<<endl;
		for(int j=0;j!=(L_columns[0]).size();++j)
		cout<<L_columns[i][j]<<" ";
		cout<<endl;
		}*/	
//	cout<<"G MATRIX"<<endl;	
//	print(G_MATRIX);
//	cout<<endl;
//	cout<<"C MATRIX"<<endl;
//	print(C_MATRIX);
	cout<<"Done generating individual matrices"<<endl;
cout<<"THE NUMBER OF VARIABLE ARE :"<<G_MATRIX.size()<<endl;
	//	cout<<endl<<"B MATRIX"<<endl;
	//	print(B_MATRIX);
	//	cout<<endl<<"u complex vector"<<endl;
	//	printcomplexvector(u_MATRIX);
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(int i=0;i!=(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	//printing tht ucomplexaray
	//	cout<<endl<<"ucomplexarray"<<endl;
	//	for(int i=0;i!=(u_MATRIX).size();++i)
	//		cout<<uarraymatrix[i][0]<<endl;
	//converting the B VECTOR MATRIX to B ARRAY MATRIX
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(int i=0;i!=(B_MATRIX).size();++i)
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

//comment for original solving starts here	
	  fstream result11("MATLAB_result_at_output_node_at_multiple_frequencies.txt",fstream::out);
	  fstream difference11("result_at_output_node_at_multiple_frequencies.txt",fstream::out);
	  if(result11.is_open() && difference11.is_open())
	  {
	  double frequency;
	  int timer1=0;
	  while(timer1<points)
	  {
	  cout<<"original"<<timer1+1<<endl;	
	  frequency=frequency_min+(interval*timer1);
	  result11<<frequency<<'\t';
	  difference11<<frequency<<'\t';
	//combining G_MATRIX and C_MATRIX into GplussC_MATRIX
	complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
	for(int i=0;i!=(G_MATRIX).size();++i){
	for(int j=0;j!=(G_MATRIX).size();++j)
	{
	GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
	}
	}
	//printing GplussC_MATRIX
	//	cout<<endl<<"GplussC_MATRIX"<<endl;
	//	for(int i=0;i!=(G_MATRIX).size();++i){
	//		//		<<cout<<i<<" th line: ";
	//		for(int j=0;j!=(G_MATRIX).size();++j)
	//			//	{if(i==2&&j==2)cout<<"this is the point"<<endl;
	//			cout<<GplussC_MATRIX[i][j]<<" ";
	//		cout<<";"<<endl;			
	//	}
	//	cout<<endl;

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
for(int i=0;i!=n;++i)
for(int j=0;j!=nrhs;++j)
original_solution[i][j]=b[i+j*ldb];
for(int i=0;i!=G_MATRIX.size();++i)
{
	for(int k=0;k!=output_node.size();++k)
	{
		if(i==(output_node[k]-1)){
			difference11<<original_solution[i][0]<<'\t';
			if((original_solution[i][0]).imag()>=0)
				result11<<(original_solution[i][0]).real()<<'+'<<(original_solution[i][0]).imag()<<'i'<<'\t';
			else
				result11<<(original_solution[i][0]).real()<<(original_solution[i][0]).imag()<<'i'<<'\t';

		}

	}
}
difference11<<endl;
result11<<';'<<endl;

++timer1;
}//for loop of frequency ends here

}
else cout<<"unable to open file"<<endl;

cout<<"done computing Ax=b of original solution for different frequencies and writing for MATLAB READING in MATLAB_result_at_output_node_at_multiple_frequencies.txt and for computing difference in result_at_output_node_at_multiple_frequencies.txt"<<endl;
 //comment for original model solving ends here






//function to generate the X matrices containing the genrerated perpendicular columns for MOR
cout<<endl<<"STARTING THE MOR FROM HERE"<<endl;	
char answer='y';//while loop starts here
while(answer=='y'){	

	int expansion;
	cout<<"enter the number of expansion points :";
	cin>>expansion;
	double abcd[expansion],efgh[expansion],expan_freq;
	int q_val;
	for(int i=0;i<expansion;++i)

	{
		cout<<"enter the "<<i+1<<" expansion frequency: ";
		cin>>expan_freq;
		abcd[i]=expan_freq;
		cout<<"enter the value of q : ";
		cin>>q_val;
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


		cout<<"the number of rows in "<< count_Q+1<<" are: "<<X_MATRIX.size()<<" the number of columns are :"<<(X_MATRIX.at(0)).size()<<endl;

	}
	cout<<"the number of rows in X_MATRIX are: "<<X_MATRIX.size()<<" the number of columns are :"<<(X_MATRIX[0]).size()<<endl;
	//blockArnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX);
	//	cout<<"printing the X_MATIRIX"<<endl;
	//	print(X_MATRIX);
	//checking each coulumvnis orthnormal to each other		
	cout<<"done computing the orthonormal matrices using block Arnoldi algorithm"<<endl;
	cout<<"checking the mtrices generated is orthonormal or not"<<endl;
	vector<vector<double> > X_MATRIX_TRANSPOSE,IDENTITY;
	vector_transpose(X_MATRIX,X_MATRIX_TRANSPOSE);
	vectormultiplication(X_MATRIX_TRANSPOSE,X_MATRIX,IDENTITY);
	cout<<"checking the orthonormality"<<endl;
	fstream wr_o("orthonormality_check.txt",fstream::out);

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
	cout<<" wrting the results of orthonormality check in orthonormality_check.txt"<<endl;

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

	   fstream difference22("MORresult_at_output_node_at_multiple_frequencies.txt",fstream::out);
	   fstream result22("MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt",fstream::out);
	   if(result22.is_open() && difference22.is_open())
	   {
	   double frequency;
	   int timer2=0;	
	   while(timer2<points)
	   {
	   cout<<"computing MOR at "<<timer2+1<<" frequency point"<<endl;			
	   frequency=frequency_min+(interval*timer2);
	   difference22<<frequency<<'\t';
	   result22<<frequency<<'\t';
	// Combining G_cap and C_cap into G_capplussC_cap
	complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
	for(int i=0;i!=(G_cap).size();++i){
	for(int j=0;j!=(G_cap).size();++j)
	{
	G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
	}
	}
	//printing GcapplussCcap
	//		cout<<endl<<"G_capplussC_cap"<<endl;
	//		for(int i=0;i!=(G_cap).size();++i){
	//			for(int j=0;j!=(G_cap).size();++j)
	//			{
	//				cout<<G_capplussC_cap[i][j]<<" ";
	//			}
	//			cout<<endl;
	//		}

	//copying elements of GcapplusCcap in a one dimension array acap
	complex<double> a_cap[(G_cap.size())*(G_cap.size())];
	int k_cap=0;
	for(int i=0;i!=(G_cap).size();++i)
	for(int j=0;j!=(G_cap).size();++j)
	a_cap[k_cap++]=G_capplussC_cap[j][i];
	//		cout<<endl;

	//printing the one dimension a_cap matrix
	//		cout<<endl<<" one dimension 'a_cap' matrix"<<endl;
	//		for(int i=0;i!=((G_cap.size())*(G_cap.size()));++i)
	//			cout<<a_cap[i]<<" ";
	//		cout<<endl;
	//
	//converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX
	complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
	for(int i=0;i!=(B_cap).size();++i)
	for(int j=0;j!=(u_MATRIX).size();++j)
	B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);
	// print the the B_cap ARRAY MATRIX
	//		cout<<endl<<"B_capCOMPLEXARRAYMATRIX"<<endl;
	//		for(int i=0;i!=(B_cap).size();++i){
	//			for(int j=0;j!=(u_MATRIX).size();++j)
	//				cout<<B_cap_COMPLEXARRAYMATRIX[i][j]<<"  ";
	//			cout<<endl;}

	//multiplying the B_COMPLEXARRAY AND uarraymatrix
	complex<double> B_capucomplexmultiple[(B_cap).size()][1];
	for(int i=0;i!=(B_cap).size();++i)
	for(int j=0;j!=(u_MATRIX).size();++j)
	B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];

	//printing the B_capucomplexmultiple
	//			cout<<endl<<"B_capucomplexmultiple"<<endl;
	//			for(int i=0;i!=(B_cap).size();++i)
	//				cout<<B_capucomplexmultiple[i][0]<<endl;

	//copying elements of B_capucomplexmultiple into b_cap one dimension matrix
	complex<double> b_cap[(B_cap).size()];
	for(int i=0;i!=(B_cap).size();++i)
	b_cap[i]=B_capucomplexmultiple[i][0];

	// printing the b_cap one dimension matrix
	//			cout<<"one dimension b_cap matrix"<<endl;
	//			for(int i=0;i!=(B_cap).size();++i)
	//				cout<<b_cap[i];
	//			cout<<endl;

	//computing a_capx=b_cap using cgesv routine
	//		cout<<"computing the Ax=B of the reduced equation of MOR method"<<endl;

	int n_cap=G_cap.size();
	int nrhs_cap=1;

	int lda_cap=n_cap;
	int ldb_cap=n_cap;
	int info_cap;
	int ipiv_cap[n_cap];
	//		printf( " CGESV Example Program Results\n" );
	zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
	// Check for the exact singularity 
	if( info_cap > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is Jero, so that A is singular;\n", info_cap, info_cap );
		printf( "the solution could not be computed.\n" );
		return( 1 );
	}

	//		cout<<"done computing the Ax=B of the reduced equation of MOR method"<<endl;
	//
	//converting X_matrix into X complex arary matrices

	complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
	for(int i=0;i!=X_MATRIX.size();++i)
		for(int j=0;j!=(X_MATRIX[0]).size();++j)
			X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

	//			cout<<" X MATRIX:-"<<endl;
	//			print(X_MATRIX);

	//printing the complex X complex array matrices
	//			cout<<"printing the X complex array"<<endl;
	//			for(int i=0;i!=(X_MATRIX).size();++i){
	//				for(int j=0;j!=(X_MATRIX[0]).size();++j)
	//					cout<<X_complexarray[i][j]<<"  ";
	//				cout<<endl;}

	//converting b_cap from inverse into x_cap 2 dimension array
	complex<double> x_cap[G_cap.size()][1];

	for(int i=0;i!=n_cap;++i)
		for(int j=0;j!=nrhs_cap;++j)
			x_cap[i][j]=b_cap[i+j*ldb_cap];

	//			*printing the x_cap array matrices
	//				cout<<endl<<"printing the x_cap array"<<endl;
	//				for(int i=0;i!=n_cap;++i)
	//				{
	//					for(int j=0;j!=nrhs_cap;++j)
	//						cout<<x_cap[i][j]<<" ";
	//					cout<<endl;
	//				}
	//multiplying the X_matrix with x_cap
	//		cout<<"computing comparison matirces by mutiplting the MOR solution with orthonormal matrices"<<endl;
	complex<double>MORsolutioncompare[X_MATRIX.size()][1];
	for(int i=0;i!=X_MATRIX.size();++i)
		for(int j=0;j!=(X_MATRIX[0]).size();++j)
			MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];
	for(int i=0;i!=G_MATRIX.size();++i)
	{
		for(int k=0;k!=output_node.size();++k)
		{
			if(i==(output_node[k]-1)){
				difference22<<MORsolutioncompare[i][0]<<'\t';							
				if((MORsolutioncompare[i][0]).imag()>=0)              
					result22<<(MORsolutioncompare[i][0]).real()<<'+'<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';
				else
					result22<<(MORsolutioncompare[i][0]).real()<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';

			}
		}
	}
	difference22<<endl;				
	result22<<';'<<endl;
	++timer2;			
}
difference22.close();
result22.close();

}
else cout<<"unable top open file"<<endl;

cout<<"done computing Ax=b of MOR solution for different frequencies and writing for MATLAB READING in MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt and for computing differnce in MORresult_at_output_node_at_multiple_frequencies.txt"<<endl;		
//cout<<"the size of Original solution is :"<<original_solution.size()<<" and the size of MOR solution compare is :"<<MORsolutioncompare.size()<<endl;
//	out<<" the original and MOR solution written in result.out"<<endl;

cout<<"do you want to continue for new value of Q(y/n):";
cin>>answer;
}//whille loop

cout<<"end of program"<<endl;
return 0;
}

/* Auxiliary routine: printing a matrix */
int print_matrix( char* desc, int m, int n, complex<double> *a, int lda ) 
{
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%f,%f)", a[i+j*lda].real(), a[i+j*lda].imag() );
		printf( "\n" );
	}
return 1;
}


/* Auxiliary routine: printing a vector of integers */
int print_int_vector( char* desc, int n, int* a ) 
{
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
return 2;
}




/*Function for printing the 2 dimension vector*/
int print( vector<vector<double> > &p)
{
	for(int i=0; i<(p.size());++i)
{
		for(int j=0; j<((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;
}
return 3;
}

/* Function for printing one dimension complex vector*/
int printcomplexvector(vector<complex<double> >&p)
{
	for(int i=0; i<(p.size());++i)
{
		// for(int j=0; j!=((p[i]).size());++j)
		cout<<p[i]<<" ";
		cout<<";"<<endl;
}
return 4;
}

/* INSERTING ZEROES IN G,C,B MATRIC*/
int addzeroes(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,char S, double nodeI,double nodeJ,int &count)
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
			for(int i=0; i<size1;++i)
{
				for(int j=0; j<(size-(G.size())) ;++j){
					G[i].push_back(0);
					C[i].push_back(0);
				}
}

			for(int k=0; k<(size-size1) ; ++k)
                    {
				G.push_back(vector<double>(size,0));
				C.push_back(vector<double>(size,0));
				B.push_back(vector<double>(1,0));
			}
		}

	}

	else if(S=='L'||S=='l')
{
		for(int i=0;i<(G.size());++i)
 {
			G[i].push_back(0);
			C[i].push_back(0);
}
			G.push_back(vector<double>((G.size()+1),0));
			C.push_back(vector<double>((C.size()+1),0));
			B.push_back(vector<double>(1,0));
	}

	else if(S=='J'||S=='j')
{
		int s=0;
		for(int i=0;i<B.size();++i)
{
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size())
{
			++count;
			for(int i=0;i<B.size();++i)
				B[i].push_back(0);
		}
	}
	else if(S=='V'||S=='v'){
		for(int i=0;i<(G.size());++i)
{
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
		int s=0;
		for(int i=0;i<B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(int i=0;i<B.size();++i)
				B[i].push_back(0);

		}
	}

	else if(S=='E'||S=='e'){
		for(int i=0;i<(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
	}
	else if(S=='H'||S=='h'){
		for(int i=0;i<size1;++i){
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

return 5;
}

int insertvaluesRG(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
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
return 6;
}

int insertvaluesC(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
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
return 7;
}

int insertvaluesL(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,double nodeI,double nodeJ,double value)
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
return 8;
}

int insertvaluesJ(vector<vector<double> > &B,vector<complex<double> > &U,double nodeI,double nodeJ,double value1,double value2,int &count)
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

return 9;
}





int insertvaluesV(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<complex<double> >&U,double nodeI,double nodeJ,double value1,double value2,int &count)
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
return 10;
}

int insertvaluesZ(vector<vector<double> > &G,vector<vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
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
return 11;
}

int insertvaluesE(vector<vector<double> > &G,vector<vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
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
return 12;
}



int insertvaluesH(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > & L_current_position,double nodeI,double nodeJ,double nodeK,double nodeL,double value,double row,double column)
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
return 13;
}





int AinverseB(vector<vector<double> > &Ad,vector<vector<double> > &Bd,vector<vector<double> > &Cd)
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
return 15;
}

int QR_factorization(vector<vector<double> > &a,vector<vector<double> > &E)
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
return 16;
}

int vectormultiplication(vector<vector<double> > &a,vector<vector<double> > &b,vector<vector<double> > &c)
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
return 17;
}

int vector_transpose(vector<vector<double> > &a,vector<vector<double> > &b)
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
return 18;
}

int blockArnoldi(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<vector<double> > &X,int &q_value, vector<double> &range)
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
return 19;
}

int print_eigenvalues( char* desc, int n,double* wr,double* wi ) {
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
return 20;
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

	return 20;
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
	return 21;
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
	return 22;
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
	return 24;
}


















