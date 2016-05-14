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
void print(vector<vector<double> >&);
void printcomplexvector(vector<complex<double> >&);
void addzeroes(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,char, double ,double ,size_t &);
void insertvaluesRG(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesC(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesL(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double);
void insertvaluesJ(vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,size_t &);
void insertvaluesV(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<complex<double> >&,double,double,double,double,size_t &);
void insertvaluesZ(vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double);
void insertvaluesE(vector<vector<double> > &,double,double,double,double,double);
void insertvaluesH(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,double,double,double,double,double,double,double);

#define pi 3.14
#define frequency 1e3 
extern "C" void zgesv_( int *n , int *nrhs , complex<double> *a , int *lda , int *ipiv ,complex<double> *b , int *ldb , int *info  );
void print_matrix( char* desc, int m, int n, complex<double> *a, int lda );
void print_int_vector( char* desc, int n, int* a );
extern "C" void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info );
void AinverseB(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
void QR_factorization(vector<vector<double> > &,vector<vector<double> > &);
void vectormultiplication(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
void vector_transpose(vector<vector<double> > &,vector<vector<double> > &);
void blockArnoldi(vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &,vector<vector<double> > &);
extern  "C" void dgeev_( char* jobvl, char* jobvr, int* n, double* a,int* lda,  double* wr, double* wi, double* vl, int* ldvl, double* vr, int* ldvr, double* work, int* lwork, int* info );

void print_eigenvalues( char* desc, int n, double* wr, double* wi );
int main()
{       
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
	size_t column=0,previous_section=1,new_section;
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
	for(size_t i=0;i!=input_node.size();++i)
		cout<<input_node[i]<<" ";
	cout<<endl;
	for(size_t i=0;i!=output_node.size();++i)
		cout<<output_node[i]<<" ";
	cout<<endl;



	ifstream myfile("text.txt");	
	cout<<"Generating the individual Matrices"<<endl;
	while(getline(myfile,line))/* && !line.empty())*/{

		for(size_t i=0;i!=(line.size());++i){
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
		for(size_t noofelements=0;noofelements!=read.size();++noofelements)
		{

			string &element=read[noofelements];
			//determining the which is string and which is value
			size_t count=0;
			for(size_t i=0;i!=element.size();++i)
			{
				if(isdigit(element[i])||ispunct(element[i])||(isdigit(element[i-1])&&element[i]=='e')||(isdigit(element[i-1])&&element[i]=='E'))
					++count;
			}
			//	cout<<"the count is "<<count<<endl;			
			if(count==element.size())
			{
				//taking each values and determining where is comma and word e or E
				size_t comma_place;
				size_t comma_counter=0;
				size_t e_counter=0;
				vector<size_t>e_place;
				for(size_t i=0;i!=element.size();++i)
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
					for(size_t z=0;z<comma_place;++z)
						s+=element[z];                                
					double f;
					istringstream(s)>>f;
					values.push_back(f);
					s.clear();
					for(size_t z=(comma_place+1);z!=element.size();++z)
						s+=element[z];
					double g;
					istringstream(s)>>g;
					values.push_back(g);
					s.clear();
				}

				else if(comma_counter==0&& e_counter==1)
				{
					double number,power,scientific;
					for(size_t z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number;
					s.clear();
					for(size_t z=(e_place[0]+1);z!=element.size();++z)
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
					for(size_t z=0;z!=e_place[0];++z)
						s+=element[z];
					istringstream(s)>>number1;
					s.clear();
					for(size_t z=(e_place[0]+1);z!=comma_place;++z)
						s+=element[z];
					istringstream(s)>>power1;
					s.clear();
					scientific1=number1*pow(10,power1);
					values.push_back(scientific1);
					for(size_t z=comma_place+1;z!=e_place[1];++z)
						s+=element[z];
					istringstream(s)>>number2;
					s.clear();
					for(size_t z=(e_place[1]+1);z!=element.size();++z)
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
						for(size_t z=0;z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>f;
						s.clear();
						values.push_back(f);
						for(size_t z=comma_place+1;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(size_t z=(e_place[0]+1);z!=element.size();++z)
							s+=element[z];
						istringstream(s)>>power;
						scientific=number*pow(10,power);
						values.push_back(scientific);
					}
					else if(comma_place>e_place[0])
					{
						double f,number,power,scientific;
						for(size_t z=0;z!=e_place[0];++z)
							s+=element[z];
						istringstream(s)>>number;
						s.clear();
						for(size_t z=(e_place[0]+1);z!=comma_place;++z)
							s+=element[z];
						istringstream(s)>>power;
						s.clear();						
						scientific=number*pow(10,power);
						values.push_back(scientific);
						for(size_t z=comma_place+1;z!=element.size();++z)
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
		//		  for(size_t i=0;i!= values.size();++i)
		//		  cout<<values[i]<<endl;

//		new_word=name[0][0];
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
		else if(name[0][0]=='C' || name[0][0]=='c')
			insertvaluesC(G_MATRIX,C_MATRIX,B_MATRIX,values[0],values[2],values[4]);
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
			insertvaluesE(G_MATRIX,values[0],values[2],values[4],values[6],values[8]);
		else if(name[0][0]=='H' || name[0][0]=='h')
			insertvaluesH(G_MATRIX,C_MATRIX,L_columns,values[0],values[2],values[4],values[6],values[8],values[10],values[12]);	


		/*For debugging the complex values are stores in values
		  for(size_t r=0; r!=values.size();++r)
		  cout<<values[r]<<" ";
		  cout<<endl;*/
		values.clear();
		name.clear();
		read.clear(); // to clear the read vector*/

	//	old_word=new_word;
		//cout<<"old_word is "<<old_word<<endl;
	}
	/*	cout<<"columns of L of respective of sections are: "<<endl;
		cout<<"size of L_columns matrices is :"<<L_columns.size()<<endl;	
		for(size_t i=0;i!=L_columns.size();++i)
		{
		cout<<"the section is: "<<i<<endl;
		for(size_t j=0;j!=(L_columns[0]).size();++j)
		cout<<L_columns[i][j]<<" ";
		cout<<endl;
		}*/	
//		cout<<"G MATRIX"<<endl;	
//		print(G_MATRIX);
//		cout<<endl;
//		cout<<"C MATRIX"<<endl;
//		print(C_MATRIX);
	cout<<"Done generating individual matrices"<<endl;
//		cout<<endl<<"B MATRIX"<<endl;
//		print(B_MATRIX);
//		cout<<endl<<"u complex vector"<<endl;
//		printcomplexvector(u_MATRIX);
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(size_t i=0;i!=(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	/*printing tht ucomplexaray*/
	//	cout<<endl<<"ucomplexarray"<<endl;
	//	for(size_t i=0;i!=(u_MATRIX).size();++i)
	//		cout<<uarraymatrix[i][0]<<endl;
	/*converting the B VECTOR MATRIX to B ARRAY MATRIX*/
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(size_t i=0;i!=(B_MATRIX).size();++i)
		for(size_t j=0;j!=(u_MATRIX).size();++j)
			B_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_MATRIX[i][j],0);
	/* print the the B ARRAY MATRIX*/
	//	cout<<endl<<"B_COMPLEXARRAYMATRIX"<<endl;
	//	for(size_t i=0;i!=(B_MATRIX).size();++i){
	//		for(size_t j=0;j!=(u_MATRIX).size();++j)
	//			cout<<B_COMPLEXARRAYMATRIX[i][j]<<"  ";
	//		cout<<endl;}

	/*multiplying the B_COMPLEXARRAY AND uarraymatrix*/
	complex<double> Bucomplexmultiple[(B_MATRIX).size()][1];
	for(size_t i=0;i!=(B_MATRIX).size();++i)
		for(size_t j=0;j!=(u_MATRIX).size();++j)
			Bucomplexmultiple[i][0]+=B_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];
	//		/*printing the Bucomplexmultiple*/
	//		cout<<endl<<"Bucomplexmultiple"<<endl;
	//		for(size_t i=0;i!=(B_MATRIX).size();++i)
	//			cout<<Bucomplexmultiple[i][0]<<endl;
	/*combining G_MATRIX and C_MATRIX into GplussC_MATRIX*/
	complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
	for(size_t i=0;i!=(G_MATRIX).size();++i){
		for(size_t j=0;j!=(G_MATRIX).size();++j)
		{
			GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
		}
	}
	/*printing GplussC_MATRIX*/
	//	cout<<endl<<"GplussC_MATRIX"<<endl;
	//	for(size_t i=0;i!=(G_MATRIX).size();++i){
	//		//		<<cout<<i<<" th line: ";
	//		for(size_t j=0;j!=(G_MATRIX).size();++j)
	//			//	{if(i==2&&j==2)cout<<"this is the point"<<endl;
	//			cout<<GplussC_MATRIX[i][j]<<" ";
	//		cout<<";"<<endl;			
	//	}
	//	cout<<endl;

	/*copying elements of GplusC_MATRIX in a one dimension array a*/
	complex<double> a[(G_MATRIX.size())*(G_MATRIX.size())];
	size_t k=0;
	for(size_t i=0;i!=(G_MATRIX).size();++i)
		for(size_t j=0;j!=(G_MATRIX).size();++j)
			a[k++]=GplussC_MATRIX[j][i];

	/*printing the one dimension a matrix*/
	//		cout<<endl<<" one dimension 'a' matrix"<<endl;
	//		for(size_t i=0;i!=((G_MATRIX.size())*(G_MATRIX.size()));++i)
	//			cout<<a[i]<<" ";
	//		cout<<endl;
	//
	/*copying elements of Bucomplexmultiple into b one dimension matrix*/

	complex<double> b[B_MATRIX.size()];
	for(size_t i=0;i!=(B_MATRIX).size();++i)
		b[i]=Bucomplexmultiple[i][0];

	/* printing the b one dimension matrix*/
	//		cout<<"one dimension b matrix"<<endl;
	//		for(size_t i=0;i!=(B_MATRIX).size();++i)
	//			cout<<b[i];
	//		cout<<endl;
	/*computing ax=b using zgesv routine*/
	cout<<"computing Ax=b of orginal solution"<<endl;
	int n=G_MATRIX.size();
	int nrhs=1;
	int lda=n;
	int ldb=n;
	int info;
	int ipiv[n];
	//	printf( " ZGESV Program Results\n" );
	zgesv_( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
	/* Check for the exact singularity */
	if( info > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", info, info );
		printf( "the solution could not be computed.\n" );
		return( 1 );
	}
	/* Print solution */
	//	print_matrix( "	ORIGINAL Solution", n, nrhs, b, ldb );
	/* Print details of LU factorization */
	//		print_matrix( "Details of LU factorization", n, n, a, lda );
	/* Print pivot indices */
	//		print_int_vector( "Pivot indices", n, ipiv );

	// converting the b matrices into 2 Dimension orginal_solution matrices in order to compare with MOR result in the end.
	cout<<"done computing Ax=B of orginal solution"<<endl;
	complex<double> original_solution[n][1];
	for(size_t i=0;i!=n;++i)
		for(size_t j=0;j!=nrhs;++j)
			original_solution[i][j]=b[i+j*ldb];
//cout<<"printing the original solution"<<endl;
//for(size_t i=0;i!=n;++i)
//cout<<original_solution[i][0]<<endl;
//
//cout<<endl;

/*	//Computing poles of original equaitons
	cout<<"computing the poles of original solution"<<endl;
	vector<vector<double> >A_original_poles;
	AinverseB(G_MATRIX,C_MATRIX,A_original_poles);
	//	cout<<" A matrix"<<endl;
	//	print(A_original_poles);
	//mutilpying -1 to A_orignial_poles matrices
	for(size_t i=0;i!=A_original_poles.size();++i)
		for(size_t j=0;j!=(A_original_poles[0]).size();++j)
			A_original_poles[i][j]=-1*A_original_poles[i][j];
	//storing in one dimesion to compute eigen value
	double one_dimension_A_original_poles[A_original_poles.size()*A_original_poles.size()];
	size_t count_original_poles=0;
	for(size_t i=0;i!=A_original_poles.size();++i)
		for(size_t j=0;j!=(A_original_poles[0]).size();++j)
			one_dimension_A_original_poles[count_original_poles++]=A_original_poles[j][i];

	int N_original_poles=A_original_poles.size();
	int info_original_poles,lwork_original_poles;
	int lda_original_poles=N_original_poles,ldvl_original_poles=N_original_poles,ldvr_original_poles=N_original_poles;
	double wkopt_original_poles;
	double* work_original_poles;
	double wr_original_poles[N_original_poles],wi_original_poles[N_original_poles],vl_original_poles[ldvl_original_poles*N_original_poles],vr_original_poles[ldvr_original_poles*N_original_poles];
	cout<<"computing eigenvalues of orginal solution"<<endl;
	lwork_original_poles=-1;
	dgeev_("N","N",&N_original_poles,one_dimension_A_original_poles,&lda_original_poles,wr_original_poles,wi_original_poles,vl_original_poles,&ldvl_original_poles,vr_original_poles,&ldvr_original_poles,&wkopt_original_poles,&lwork_original_poles,&info_original_poles);
	lwork_original_poles=(int)wkopt_original_poles;
	work_original_poles=(double*)malloc(lwork_original_poles*sizeof(double));
	dgeev_("N","N",&N_original_poles,one_dimension_A_original_poles,&lda_original_poles,wr_original_poles,wi_original_poles,vl_original_poles,&ldvl_original_poles,vr_original_poles,&ldvr_original_poles,work_original_poles,&lwork_original_poles,&info_original_poles);
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		return( 1 );
	}
	//print_eigenvalues( "Eigenvalues",N_original_poles , wr_original_poles, wi_original_poles );
	cout<<"done generating the eigen values of original solution"<<endl;
	//computing the poles of orginal solution
	complex<long double> ai_original_poles[N_original_poles],number1,number2;
	for(size_t i=0;i!=N_original_poles;++i)
		ai_original_poles[i]=complex<long double>(wr_original_poles[i],wi_original_poles[i]);
	cout<<"writing the poles of orginal solution in original_poles.out file"<<endl;
	fstream result1("original_poles.out",fstream::out);
	if(result1.is_open()){
		for(size_t i=0;i!=N_original_poles;++i)
		{
			number1=conj(ai_original_poles[i]);
			number2=number1/(number1*ai_original_poles[i]);
			result1<<number2<<endl;
		}
		result1.close();
	}
	else cout<<"unable to open file"<<endl;
	cout<<"done writing the poles of orginal solution in original_poles.out file"<<endl;*/

	//function to generate the X matrices containing the genrerated perpendicular columns for MOR
	
cout<<"the number of variable are"<<G_MATRIX.size()<<endl;
cout<<endl<<"STARTING THE MOR FROM HERE"<<endl;	
	vector<vector<double> > X_MATRIX;
	cout<<"computing the orthonormal matrices using block Arnoldi algorithm"<<endl;	
	blockArnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX);
cout<<"the rows of X_MATRIX are :"<<X_MATRIX.size()<<" and the columns of X_MATRIX are :"<<(X_MATRIX[0]).size()<<endl;
	//		cout<<"printing the X_MATIRIX"<<endl;
	//		print(X_MATRIX);
	//checking each coulumvnis orthnormal to each other		
	cout<<"done computing the orthonormal matrices using block Arnoldi algorithm"<<endl;
	cout<<"checking the mtrices generated is orthonormal or not"<<endl;
	for(size_t j=0;j!=(((X_MATRIX[0]).size())-1);++j)
	{
		double check=0;
		for(size_t i=0;i!=X_MATRIX.size();++i)
		{
			check+=X_MATRIX[i][j]*X_MATRIX[i][j+1];
		}
		int check1=check; // if check not check1 is checked than we find that it is not equal to zeros du to non zero fraction which shows the accuracy of grahm scnmidt method
		cout<<" value of product of column "<<j<<" and "<<j+1<<" of E = "<<check<<endl;
		if (check1 ==0)
			cout<<"ORTHONORMAL matrix"<<endl;
		else
			cout<<"NOT ORTHONORMAL matrix"<<endl;
	}
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
/*	//computing poles of MOR method
	cout<<"computing the poles of MOR method"<<endl;
	vector<vector<double> >A_MOR_poles;
	vectormultiplication(X_transpose,A_original_poles,multiple);
	vectormultiplication(multiple,X_MATRIX,A_MOR_poles);
	multiple.clear();
	double one_dimension_A_MOR_poles[A_MOR_poles.size()*A_MOR_poles.size()];
	size_t count_MOR_poles=0;
	for(size_t i=0;i!=A_MOR_poles.size();++i)
		for(size_t j=0;j!=(A_MOR_poles[0]).size();++j)
			one_dimension_A_MOR_poles[count_MOR_poles++]=A_MOR_poles[j][i];

	int N_MOR_poles=A_MOR_poles.size();
	int info_MOR_poles,lwork_MOR_poles;
	int lda_MOR_poles=N_MOR_poles,ldvl_MOR_poles=N_MOR_poles,ldvr_MOR_poles=N_MOR_poles;
	double wkopt_MOR_poles;
	double* work_MOR_poles;
	double wr_MOR_poles[N_MOR_poles],wi_MOR_poles[N_MOR_poles],vl_MOR_poles[ldvl_MOR_poles*N_MOR_poles],vr_MOR_poles[ldvr_MOR_poles*N_MOR_poles];
	cout<<"computing eigenvalues of MOR method"<<endl;
	lwork_MOR_poles=-1;
	dgeev_("N","N",&N_MOR_poles,one_dimension_A_MOR_poles,&lda_MOR_poles,wr_MOR_poles,wi_MOR_poles,vl_MOR_poles,&ldvl_MOR_poles,vr_MOR_poles,&ldvr_MOR_poles,&wkopt_MOR_poles,&lwork_MOR_poles,&info_MOR_poles);
	lwork_MOR_poles=(int)wkopt_MOR_poles;
	work_MOR_poles=(double*)malloc(lwork_MOR_poles*sizeof(double));
	dgeev_("N","N",&N_MOR_poles,one_dimension_A_MOR_poles,&lda_MOR_poles,wr_MOR_poles,wi_MOR_poles,vl_MOR_poles,&ldvl_MOR_poles,vr_MOR_poles,&ldvr_MOR_poles,work_MOR_poles,&lwork_MOR_poles,&info_MOR_poles);
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		return( 1 );
	}
	//print_eigenvalues( "Eigenvalues",N_MOR_poles , wr_MOR_poles, wi_MOR_poles );
	cout<<"done generating the eigen values of MOR solution"<<endl;
	//computing the poles of MOR solution
	complex<long double> ai_MOR_poles[N_MOR_poles],number3,number4;
	for(size_t i=0;i!=N_MOR_poles;++i)
		ai_MOR_poles[i]=complex<long double>(wr_MOR_poles[i],wi_MOR_poles[i]);
	cout<<"writing the poles of MOR solution in MOR_poles.out file"<<endl;
	fstream result2("MOR_poles.out",fstream::out);
	if(result2.is_open()){
		for(size_t i=0;i!=N_MOR_poles;++i)
		{
			number3=conj(ai_MOR_poles[i]);
			number4=number3/(number3*ai_MOR_poles[i]);
			result2<<number4<<endl;
		}
		result2.close();
	}
	else cout<<"unable to open file"<<endl;
	cout<<"done writing the poles of MOR method in MOR_poles.out file"<<endl;*/
	/* Combining G_cap and C_cap into G_capplussC_cap*/
	complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
	for(size_t i=0;i!=(G_cap).size();++i){
		for(size_t j=0;j!=(G_cap).size();++j)
		{
			G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
		}
	}
	/*printing GcapplussCcap*/
	//		cout<<endl<<"G_capplussC_cap"<<endl;
	//		for(size_t i=0;i!=(G_cap).size();++i){
	//			for(size_t j=0;j!=(G_cap).size();++j)
	//			{
	//				cout<<G_capplussC_cap[i][j]<<" ";
	//			}
	//			cout<<endl;
	//		}

	/*copying elements of GcapplusCcap in a one dimension array acap*/
	complex<double> a_cap[(G_cap.size())*(G_cap.size())];
	size_t k_cap=0;
	for(size_t i=0;i!=(G_cap).size();++i)
		for(size_t j=0;j!=(G_cap).size();++j)
			a_cap[k_cap++]=G_capplussC_cap[j][i];
	//		cout<<endl;

	/*printing the one dimension a_cap matrix*/
	//		cout<<endl<<" one dimension 'a_cap' matrix"<<endl;
	//		for(size_t i=0;i!=((G_cap.size())*(G_cap.size()));++i)
	//			cout<<a_cap[i]<<" ";
	//		cout<<endl;
	//
	/*converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX*/
	complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
	for(size_t i=0;i!=(B_cap).size();++i)
		for(size_t j=0;j!=(u_MATRIX).size();++j)
			B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);
	/* print the the B_cap ARRAY MATRIX*/
	//		cout<<endl<<"B_capCOMPLEXARRAYMATRIX"<<endl;
	//		for(size_t i=0;i!=(B_cap).size();++i){
	//			for(size_t j=0;j!=(u_MATRIX).size();++j)
	//				cout<<B_cap_COMPLEXARRAYMATRIX[i][j]<<"  ";
	//			cout<<endl;}

	/*multiplying the B_COMPLEXARRAY AND uarraymatrix*/
	complex<double> B_capucomplexmultiple[(B_cap).size()][1];
	for(size_t i=0;i!=(B_cap).size();++i)
		for(size_t j=0;j!=(u_MATRIX).size();++j)
			B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];

	/*printing the B_capucomplexmultiple*/
	//			cout<<endl<<"B_capucomplexmultiple"<<endl;
	//			for(size_t i=0;i!=(B_cap).size();++i)
	//				cout<<B_capucomplexmultiple[i][0]<<endl;

	/*copying elements of B_capucomplexmultiple into b_cap one dimension matrix*/
	complex<double> b_cap[(B_cap).size()];
	for(size_t i=0;i!=(B_cap).size();++i)
		b_cap[i]=B_capucomplexmultiple[i][0];

	/* printing the b_cap one dimension matrix*/
	//			cout<<"one dimension b_cap matrix"<<endl;
	//			for(size_t i=0;i!=(B_cap).size();++i)
	//				cout<<b_cap[i];
	//			cout<<endl;

	/*computing a_capx=b_cap using cgesv routine*/
	cout<<"computing the Ax=B of the reduced equation of MOR method"<<endl;

	int n_cap=G_cap.size();
	int nrhs_cap=1;

	int lda_cap=n_cap;
	int ldb_cap=n_cap;
	int info_cap;
	int ipiv_cap[n_cap];
	printf( " CGESV Example Program Results\n" );
	zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
	/* Check for the exact singularity */
	if( info_cap > 0 ) {
		printf( "The diagonal element of the triangular factor of A,\n" );
		printf( "U(%i,%i) is zero, so that A is singular;\n", info_cap, info_cap );
		printf( "the solution could not be computed.\n" );
		return( 1 );
	}
	/* Print solution */
	//	print_matrix( "Solution", n_cap, nrhs_cap, b_cap, ldb_cap );
	/* Print details of LU factorization */
	//	print_matrix( "Details of LU factorization", n_cap, n_cap, a_cap, lda_cap );
	/* Print pivot indices */
	//	print_int_vector( "Pivot indices", n_cap, ipiv_cap );

	cout<<"done computing the Ax=B of the reduced equation of MOR method"<<endl;

	/*converting X_matrix into X complex arary matrices*/

	complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
	for(size_t i=0;i!=X_MATRIX.size();++i)
		for(size_t j=0;j!=(X_MATRIX[0]).size();++j)
			X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

	//			cout<<" X MATRIX:-"<<endl;
	//			print(X_MATRIX);

	/*printing the complex X complex array matrices*/
	//			cout<<"printing the X complex array"<<endl;
	//			for(size_t i=0;i!=(X_MATRIX).size();++i){
	//				for(size_t j=0;j!=(X_MATRIX[0]).size();++j)
	//					cout<<X_complexarray[i][j]<<"  ";
	//				cout<<endl;}

	/*converting b_cap from inverse into x_cap 2 dimension array*/
cout<<"the rows of X reduced are :"<<G_cap.size()<<endl;	
complex<double> x_cap[G_cap.size()][1];
	for(size_t i=0;i!=n_cap;++i)
		for(size_t j=0;j!=nrhs_cap;++j)
			x_cap[i][j]=b_cap[i+j*ldb_cap];

	//				/*printing the x_cap array matrices*/
	//				cout<<endl<<"printing the x_cap array"<<endl;
	//				for(size_t i=0;i!=n_cap;++i)
	//				{
	//					for(size_t j=0;j!=nrhs_cap;++j)
	//						cout<<x_cap[i][j]<<" ";
	//					cout<<endl;
	//				}
	/*multiplying the X_matrix with x_cap*/
	cout<<"computing comparison matirces by mutiplting the MOR solution with orthonormal matrices"<<endl;
	complex<double>MORsolutioncompare[X_MATRIX.size()][1];
	for(size_t i=0;i!=X_MATRIX.size();++i)
		for(size_t j=0;j!=(X_MATRIX[0]).size();++j)
			MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];

	/*  writing the orginal and  MOR solution in result.out*/
	fstream result3("result_at_node.out",fstream::out);
	fstream result("result.out",fstream::out);
	cout<<endl<<"********************************************************************************************************************************************************************"<<endl;	
	if(result.is_open()&& result3.is_open()){
		cout<<"writing the original and MOR soluition in result.out and result_at_node.out"<<endl;
		result<<"ORIGINAL SOLUTION"<<'\t' <<"MORsolution"<<endl;
		//	cout<<"ORIGINAL SOLUTION"<<'\t' <<"MORsolution"<<endl;
		for(size_t i=0;i!=X_MATRIX.size();++i){
			//	cout<<original_solution[i][0]<<'\t'<<MORsolutioncompare[i][0]<<endl;		
			result<<original_solution[i][0]<<'\t'<<MORsolutioncompare[i][0]<<endl;

			for(size_t j=0;j!=input_node.size();++j){
				if(i==(input_node[j]-1))
					result3<<original_solution[i][0]<<'\t'<<MORsolutioncompare[i][0]<<"INPUT NODE : "<<input_node[j]<<endl;
			}
			for(size_t k=0;k!=output_node.size();++k){
				if(i==(output_node[k]-1))
					result3<<original_solution[i][0]<<'\t'<<MORsolutioncompare[i][0]<<"OUTPUT NODE : "<<output_node[k]<<endl;
			}

		}
		result.close();
result3.close();
	}
	else cout<<"unable to open file"<<endl;
	//cout<<"the size of Original solution is :"<<original_solution.size()<<" and the size of MOR solution compare is :"<<MORsolutioncompare.size()<<endl;
	cout<<" the original and MOR solution written in result.out"<<endl;


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
	for(size_t i=0; i!=(p.size());++i){
		for(size_t j=0; j!=((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;}
}

/* Function for printing one dimension complex vector*/
void printcomplexvector(vector<complex<double> >&p)
{
	for(size_t i=0; i!=(p.size());++i){
		// for(size_t j=0; j!=((p[i]).size());++j)
		cout<<p[i]<<" ";
		cout<<";"<<endl;}
}

/* INSERTING ZEROES IN G,C,B MATRIC*/
void addzeroes(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,char S, double nodeI,double nodeJ,size_t &count)
{
	size_t size=(nodeI>nodeJ)?nodeI:nodeJ;
	size_t size1=G.size();//or size1=C.size,had to add this syntax because when the vectros were not empty,it was not able to take the size of the vectors.	
	size_t size2=B.size();	
	if(S=='R'||S=='r'||S=='G'||S=='g'||S=='C'||S=='c'){
		if(G.empty() && C.empty() && B.empty()  ){
			G=vector<vector<double> >(size,vector<double>(size,0));
			C=vector<vector<double> >(size,vector<double>(size,0));
			B=vector<vector<double> >(size,vector<double>(1,0));}
		else if(!G.empty() && !C.empty() && !B.empty() && size>G.size() && size>C.size()){
			for(size_t i=0; i!=size1;++i){
				for(size_t j=0; j!=(size-(G.size())) ;++j){
					G[i].push_back(0);
					C[i].push_back(0);
				}}
			for(size_t k=0; k!=(size-size1) ; ++k){
				G.push_back(vector<double>(size,0));
				C.push_back(vector<double>(size,0));
				B.push_back(vector<double>(1,0));
			}
		}

	}


	else if(S=='L'||S=='l'){
		for(size_t i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);}
			G.push_back(vector<double>((G.size()+1),0));
			C.push_back(vector<double>((C.size()+1),0));
			B.push_back(vector<double>(1,0));
	}

	else if(S=='J'||S=='j'){
		size_t s=0;
		for(size_t i=0;i!=B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(size_t i=0;i!=B.size();++i)
				B[i].push_back(0);
		}
	}
	else if(S=='V'||S=='v'){
		for(size_t i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
		size_t s=0;
		for(size_t i=0;i!=B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(size_t i=0;i!=B.size();++i)
				B[i].push_back(0);

		}
	}

	else if(S=='E'||S=='e'){
		for(size_t i=0;i!=(G.size());++i){
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(vector<double>((G.size()+1),0));
		C.push_back(vector<double>((C.size()+1),0));
		B.push_back(vector<double>((count+1),0));
	}
	else if(S=='H'||S=='h'){
		for(size_t i=0;i!=size1;++i){
			//	for(size_t j=0;j!=2;++j){
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
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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
	for(size_t i=0;i!=(C.size());++i){
		for(size_t j=0;j!=((C[i]).size());++j){
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
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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

void insertvaluesJ(vector<vector<double> > &B,vector<complex<double> > &U,double nodeI,double nodeJ,double value1,double value2,size_t &count)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(size_t i=0;i!=B.size();++i)
	{
		if(i==I && I>=0)
			B[i][count]=B[i][count]+1;
		else if(i==J && J>=0)
			B[i][count]=B[i][count]-1;
	}

	complex<double> b= complex<double>(value1,value2);
	U.push_back(b);


}





void insertvaluesV(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<complex<double> >&U,double nodeI,double nodeJ,double value1,double value2,size_t &count)
{//cout<<"in the insert function"<<endl;
	double I=nodeI-1;
	double J=nodeJ-1;
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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

void insertvaluesE(vector<vector<double> > &G,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
{
	size_t last=(G.size())-1;
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double NCplus=nodeK-1;
	double NCminus=nodeL-1;
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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
	size_t last=(G.size())-1;
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double section=row-1;
	double dependent_line=column-1;
	//cout<<"section :"<<section<<" line :"<<dependent_line<<endl;
	for(size_t i=0;i!=(G.size());++i){
		for(size_t j=0;j!=((G[i]).size());++j){
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
	size_t sizeA=Ad.size();
	int Nd=sizeA;
	size_t colB=(Bd[0]).size();
	int NRHSd=colB;
	int LDAd =Nd;
	int LDBd=Nd;
	int IPIVd[Nd];
	int INFOd;
	double ad[LDAd*Nd],bd[LDBd*NRHSd];
	int kd=0;
	int ld=0;
	for(size_t i=0;i!=(Ad[0]).size();++i)
		for(size_t j=0;j!=Ad.size();++j)
			ad[kd++]=Ad[j][i];
	for(size_t i=0;i!=(Bd[0]).size();++i)
		for(size_t j=0;j!=Bd.size();++j)
			bd[ld++]=Bd[j][i];
	dgesv_( &Nd, &NRHSd, ad, &LDAd, IPIVd, bd, &LDBd, &INFOd );
	if( INFOd > 0 ) {
		cout<<"The diagonal element of the triangular factor of A,\n";
		cout<< "U("<<INFOd<<","<<INFOd<<" is zero, so that A is singular;\n";
		cout<<( "the solution could not be computed.\n" );
		//	return(1);
	}
	vector<double>tempd;
	for(size_t i=0;i!=Nd;++i)
	{
		for(size_t j=0;j!=NRHSd;++j)
			tempd.push_back(bd[i+j*LDBd]);
		Cd.push_back(tempd);
		tempd.clear();
	}
}

void QR_factorization(vector<vector<double> > &a,vector<vector<double> > &E)
{
	vector<double> u;
	for(size_t i=0;i!=(a[0]).size();++i){
		for(size_t j=0;j!=a.size();++j)
			u.push_back(a[j][i]);
		//		cout<<"printing"<< i+1 <<"u vector"<<endl;
		//		for(size_t r=0;r!=u.size();++r)
		//			cout<<u[r]<<" ";
		//	cout<<endl;

		if(i==0)
		{
			double sum=0;
			for(size_t k=0;k!=u.size();++k)
				sum+=u[k]*u[k];
			for(size_t k=0;k!=u.size();++k){
			//	double root=sqrt(sum);
			//	double value=u[k]/root;
				E.push_back(vector<double>(1,(u[k]/sqrt(sum))));
			}
			u.clear();
		}
		else{
			double sum=0;
			for(size_t k=0;k!=i;++k)
			{
				double prod=0;
				for(size_t z=0;z!=u.size();++z)
					prod+=a[z][i]*E[z][k];
				for(size_t z=0;z!=u.size();++z)
					u[z]=u[z]-(prod*E[z][k]);
			}
			for(size_t k=0;k!=u.size();++k)
				sum+=u[k]*u[k];
			for(size_t k=0;k!=u.size();++k){
			//	double root=sqrt(sum);
			//	double value=u[k]/root;
				E[k].push_back(u[k]/sqrt(sum));
			}
			u.clear();
		}
	}
	/*	for(size_t j=0;j!=((E[0].size())-1);++j)
		{
		double check=0;
		for(size_t i=0;i!=E.size();++i)
		{
		check+=E[i][j]*E[i][j+1];
		}
		int check1=check; // if check not check1 is checked than we find that it is not equal to zeros du to non zero fraction which shows the accuracy of grahm scnmidt method
		cout<<" value of product of column "<<j<<" and "<<j+1<<" of E = "<<check<<endl;
		if (check1 ==0)
		cout<<"ORTHONORMAL matrix"<<endl;
		else
		cout<<"NOT ORTHONORMAL matrix"<<endl;
		}

		cout<<endl;*/
}

void vectormultiplication(vector<vector<double> > &a,vector<vector<double> > &b,vector<vector<double> > &c)
{
	if((a[0]).size()!=b.size())
		cout<<"columns of 1st matrix and row of 2nd matrix does not match"<<endl;
	else{
		vector<double>temp;
		size_t row=a.size();
		size_t column=(b[1]).size();
		size_t inner=b.size();
		for(size_t i=0;i!=row;++i){
			for(size_t j=0;j!=column;++j){
				double sum=0;
				for(size_t k=0;k!=inner;++k){
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
	for(size_t i=0;i!=(a[0]).size();++i)
	{
		for(size_t j=0;j!=a.size();++j)
		{
			temp.push_back(a[j][i]);
		}
		b.push_back(temp);
		temp.clear();
	}
}

void blockArnoldi(vector<vector<double> > &G,vector<vector<double> > &C,vector<vector<double> > &B,vector<vector<double> > &X)
{
	vector<double> range,insert;
	vector<vector<double> > R,Q,V,X_temp,Z,X_transpose,H,XH;
	AinverseB(G,B,R);
//cout<<"the inverse matrix for the first time in block Arnoldi algorithm"<<endl;
//print(R);
	QR_factorization(R,Q);
	size_t q,n,input1,input2;
	X=Q;
//	cout<<"printing the X 0"<<endl;
//	print(X);	
	Q.clear();
	R.clear();
	range.push_back((X[0]).size());
	cout<<"Enter the value of q"<<endl;
	cin>>q;
	if((q%((B[0]).size()))==0){
		n=q/((B[0]).size());
cout<<"the value of n is :"<<n<<endl;
}

	else{
		n=floor(q/((B[0]).size()))+1;
cout<<"the value of n is :"<<n<<endl;
}
	for(size_t kA=1;kA<=n;++kA)
	{
		input1=kA-1;
		if(input1==0){
			for(size_t l=0;l!=X.size();++l)
			{
				for(size_t m=0;m!=range[input1];++m)
					insert.push_back(X[l][m]);
				X_temp.push_back(insert);
				insert.clear();
			}
		}
		else{
			for(size_t l=0;l!=X.size();++l)
			{
				for(size_t m=range[input1-1];m!=range[input1];++m)
					insert.push_back(X[l][m]);
				X_temp.push_back(insert);
				insert.clear();
			}
		}
		vectormultiplication(C,X_temp,V);
		X_temp.clear();
		AinverseB(G,V,Z);
		V.clear();
		for(size_t jA=1;jA<=kA;++jA)
		{
			input2=kA-jA;
			if(input2==0){
				for(size_t l=0;l!=X.size();++l)
				{
					for(size_t m=0;m!=range[input2];++m)
						insert.push_back(X[l][m]);
					X_temp.push_back(insert);
					insert.clear();
				}
			}
			else{
				for(size_t l=0;l!=X.size();++l)
				{
					for(size_t m=range[input2-1];m!=range[input2];++m)
						insert.push_back(X[l][m]);
					X_temp.push_back(insert);
					insert.clear();
				}
			}
			vector_transpose(X_temp,X_transpose);
			vectormultiplication(X_transpose,Z,H);
			vectormultiplication(X_temp,H,XH);
			for(size_t iA=0;iA!=Z.size();++iA)
				for(size_t jA=0;jA!=(Z[0]).size();++jA)
					Z[iA][jA]=Z[iA][jA]-XH[iA][jA];
			X_temp.clear();
			X_transpose.clear();
			XH.clear();
			H.clear();
		}
//printing the z for arnoldi process for SISO
for(size_t iA=0;iA!=Z.size();++iA)
                                for(size_t jA=0;jA!=(Z[0]).size();++jA)
cout<<Z[iA][jA]<<" ";
cout<<endl;
		QR_factorization(Z,Q);
		for(size_t iA=0;iA!=Q.size();++iA)
			for(size_t jA=0;jA!=(Q[0]).size();++jA)
				(X[iA]).push_back(Q[iA][jA]);
		range.push_back((X[0]).size());
			cout<<"printing the X"<<kA<<endl;
		//	print(Q);		
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

