#include<iostream>
#include<fstream>
#include<vector>
#include<algorithm>
#include"stdlib.h"
#include"stdio.h"
#include<cmath>
#define tr 0.1e-9
#define n 5 
#define l 10  // in cms
//#define m 2 // write code to automatically generate the m
using namespace std;
int print( vector<vector<double> > );
int vectormultiplication(vector<vector<double> > ,vector<vector<double> > ,vector<vector<double> > &);
extern "C" void dgeev_( char* jobvl, char* jobvr, int* N, double* a,int* lda,  double* wr, double* wi, double* vl, int* ldvl,double* vr, int* ldvr, double* work, int* lwork, int* info );
/* Auxiliary routines prototypes */
void print_eigenvalues( char* desc, int N,  double* wr, double* wi );

int main(int argc, char *argv[])
{
	cout<<" DID YOU CHANGE THE VALUE OF n"<<endl;

	fstream capacitance,resistance,inductance,admittance;
	capacitance.open("C_5.txt",fstream::in);
	resistance.open("R_5.txt",fstream::in);
	inductance.open("L_5.txt",fstream::in);
	admittance.open("G_5.txt",fstream::in);
	vector<vector<double> >C,R,L,G;
	vector<double>store1,store2,store3,store4;
	double cap_value,res_value,ind_value,adm_value;
	int s1=0;
	if(capacitance.is_open() && resistance.is_open() && inductance.is_open() && admittance.is_open())
	{
		while(capacitance>>cap_value && resistance>>res_value && inductance>>ind_value && admittance>>adm_value)
		{
			store1.push_back(cap_value/100);
			store2.push_back(res_value/100);
			store3.push_back(ind_value/100);
			store4.push_back(adm_value/100);
			++s1;
			if(s1==n){
				C.push_back(store1);
				store1.clear();
				R.push_back(store2);
				store2.clear();
				L.push_back(store3);
				store3.clear();
				G.push_back(store4);
				store4.clear();
				s1=0;
			}


		}
	}
	else cout<<"Unable to open file"<<endl;


	cout<<" C Matrix:"<<endl;
	print(C);
	cout<<" R Matrix: "<<endl;
	print(R);
	cout<<" L Matrix:"<<endl;
	print(L);
	cout<<" G Matrix: "<<endl;
	print(G);
	//computing the values between particular node and ground
	for(int i=0;i!=n;++i)
		for(int j=0;j!=n;++j)
			if(i!=j)
				C[i][i]+=C[i][j];

	cout<<" C modified Matrix: "<<endl;
	print(C);
	//finding the value of q by finding maximum eigen value
	vector<vector<double> > multiple;
	vectormultiplication(L,C,multiple);
	int N=multiple.size();
	int lda=N,ldvl=N,ldvr=N,info,lwork;
	double wkopt;
	double* work;
	double wr[N],wi[N];
	double vl[ldvl*N],vr[ldvr*N];
	double a_one_dimension[N*N];
	int count1=0;
	for(int i=0;i!=N;++i)
		for(int j=0;j!=N;++j)
			a_one_dimension[count1++]=multiple[j][i];

//	cout<<"one dimension a matrices"<<endl;
//	for(int i=0;i!=(N*N);++i)
//		cout<<a_one_dimension[i]<<" ";
	cout<<endl;
	lwork = -1;
	dgeev_( "N", "N", &N, a_one_dimension, &lda, wr, wi, vl, &ldvl, vr, &ldvr,&wkopt, &lwork, &info );
	lwork = (int)wkopt;
	work = (double*)malloc( lwork*sizeof(double) );
	dgeev_( "N", "N", &N, a_one_dimension, &lda, wr, wi, vl, &ldvl, vr, &ldvr,work, &lwork, &info );
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		return( 1 );
	}
//	cout<<"the LC marix :"<<endl;
//	print(multiple);
//	cout<<"printing directly:"<<endl;
//	for(int i=0;i!=N;++i)
//		cout<<"("<<wr[i]<<","<<wi[i]<<")"<<" ";
//	cout<<endl;
//	cout<<"the square root of eigen values"<<endl;	
//	for(int i=0;i!=N;++i)
//		cout<<sqrt(wr[i])<<" ";
//	cout<<endl;
//	cout<<" the number of section of each"<<endl;
//	for(int i=0;i!=N;++i)
//		cout<<20*sqrt(wr[i])*l/tr<<" ";
//	cout<<endl;


int R_count=0,C_count=0,L_count=0,G_count=0,E_count=0,Z_count=0;
	double eig[N],max,m;
	for(int i=0;i!=N;++i)
	{
		eig[i]=20*sqrt(wr[i])*l/tr;
	}
	max=eig[0];
	for(int i=1;i!=N;++i)
	{
		if(max<eig[i])
			max=eig[i];
	}
	m=ceil(max);
	cout<<"the max is"<<max<<endl;
//m=1;
	cout<<"the number of section are: "<<m<<endl;
	cout<<" enter the number of Multi conductor interconnects : ";
	int number;
	cin>>number;
	for(int no=0;no<number;++no)
	{
		cout<<"enter the staring node :";
		int s,s_1;
		cin>>s_1;
		s=s_1-1;


		vector<vector<vector<double> > > transmission_line;
		vector<vector<double> >twoD;
		vector<double> oneD,temp;

		for(int section=0;section!=m;++section)
		{
			if(section==0)
			{
				for(int line=0;line<n;++line)
				{

					for(int nodes=0;nodes<(n+2);++nodes)
					{	
						temp.push_back(++s);
					}
					twoD.push_back(temp);
					temp.clear();
				}
				transmission_line.push_back(twoD);
				twoD.clear();
			}
			else
			{
				for(int line=0;line<n;++line)
				{
					//		for(int nodes=0;nodes<n+1;++nodes)
					//			temp.push_back(++s);
					int section_size=(transmission_line[section-1][line]).size();
					temp.push_back(transmission_line[section-1][line][section_size-1]);
					for(int nodes=0;nodes<n+1;++nodes)
					{
						temp.push_back(++s);
					}	
					twoD.push_back(temp);
					temp.clear();
				}
				transmission_line.push_back(twoD);
				twoD.clear();
			}

		}
		cout<<"the end terminal is"<<s<<endl;
		cout<<"representing a TRANSMISSION line nodes: "<<endl;

		//printing as well as writing the input nodes
		fstream input("input_terminal.txt",fstream::out);
		cout<<"the input terminals are :"<<endl;
		for(int line=0;line<n;++line)
		{
			cout<<transmission_line[0][line][0]<<endl;
			input<<transmission_line[0][line][0]<<endl;
		}
		input.close();

		//printing as well as wrtining the output nodes
		fstream output("output_terminal.txt",fstream::out);
		cout<<"the output terminals are :"<<endl;
		for(int line=0;line<n;++line)
		{
			cout<<transmission_line[m-1][line][((transmission_line[m-1][line]).size())-1]<<endl;
			output<<transmission_line[m-1][line][((transmission_line[m-1][line]).size())-1]<<endl;
		}
		output.close();


/*		for(int section=0;section<m;++section)
		{
			cout<<" for section :"<<section+1<<endl;
			for(int line=0;line<n;++line)
			{
				for(int nodes=0;nodes<n+2;++nodes)
				{
					cout<<transmission_line[section][line][nodes]<<" ";
				}
				cout<<endl;
			}
	}*/
		cout<<"delat_x is: "<<l/m;

		fstream myfile;
		if(no==0)
			myfile.open("time_multi_10.txt",fstream::out);
		else
			myfile.open("time_multi_10.txt",fstream::out|fstream::app);

		if(myfile.is_open())
		{
			//writing R G C and L in input file
			//long double ans=R[0][0]*l/m;
			//cout<<ans<<endl;
			cout<<"file is being written"<<endl;
			for(int section=0;section<transmission_line.size();++section)
			{
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
					myfile<<"R"<<++R_count<<" "<<transmission_line[section][line][0]<<" "<<transmission_line[section][line][1]<<" "<<(R[line][line]*l/m)<<endl;
				}
			}
			for(int section=0;section<transmission_line.size();++section)
			{
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
//					myfile<<"G"<<++G_count<<" "<<transmission_line[section][line][((transmission_line[section][line]).size())-1]<<" "<<0<<" "<<(G[line][line]*l/m)<<endl;
				}
			}
			for(int section=0;section<transmission_line.size();++section)
			{
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
					myfile<<"C"<<++C_count<<" "<<transmission_line[section][line][((transmission_line[section][line]).size())-1]<<" "<<0<<" "<<(C[line][line]*l/m)<<endl;
				}
			}
			for(int section=0;section<transmission_line.size();++section)
			{
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
					myfile<<"L"<<++L_count<<" "<<transmission_line[section][line][1]<<" "<<transmission_line[section][line][2]<<" "<<(L[line][line]*l/m)<<endl;
				}
			}		// writng the VCCS 
			for(int section=0;section<transmission_line.size();++section)
			{
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
					vector<double> gain,k;
					for(int i=0;i<(C[line]).size();++i)
					{
						if(line!=i)
						{
							gain.push_back(C[line][i]);
						}
					}
					for(int i=0;i<(transmission_line[0]).size();++i)
					{
						if(line!=i)
						{  
							k.push_back(transmission_line[section][i][((transmission_line[section][line]).size())-1]);
						}            
					}


					for(int i=0;i<(n-1);++i)
					{
  //  myfile<<"C"<<++C_count<<" "<<transmission_line[section][line][((transmission_line[section][line]).size())-1]<<" "<<k[i]<<" "<<(-gain[i]*l/m)<<endl;
myfile<<"Z"<<++Z_count<<" "<<transmission_line[section][line][((transmission_line[section][line]).size())-1]<<" "<<0<<" "<<transmission_line[section][line][((transmission_line[section][line]).size())-1]<<" "<<k[i]<<" "<<(-gain[i]*l/m)<<endl;
	

				}			
				}
			}
			// wrting VCVS
			for(int section=0;section<transmission_line.size();++section)
				for(int line=0;line<(transmission_line[0]).size();++line)
				{
					vector<double> gain1,gain2,k1,k2;
					for(int i=0;i<(L[line]).size();++i)
					{
						if(line!=i){
							gain1.push_back(L[line][i]);
							gain2.push_back(L[i][i]);
						}
					}
					for(int i=0;i<(transmission_line[0]).size();++i)
					{
						if(line!=i)
						{
							k1.push_back(transmission_line[section][i][1]);
							k2.push_back(transmission_line[section][i][2]);
						}
					}
					/*cout<<" printing the off diagnol elements :"<<endl;
					  for(int i=0;i<gain1.size();++i)
					  cout<< gain1[i]<<" ";
					  cout<<endl;
					  for(int i=0;i<gain1.size();++i)
					  cout<< gain2[i]<<" ";
					  cout<<endl;*/
					for(int i=0;i<(n-1);++i)
					{
						myfile<<"E"<<++E_count<<" "<<transmission_line[section][line][2+i]<<" "<<transmission_line[section][line][3+i]<<" "<<k1[i]<<" "<<k2[i]<<" "<<((gain1[i]/gain2[i]))<<endl;

					}
				}
			cout<<"file is written"<<endl;
			myfile.close();
		}
		else cout<<"unable to open file"<<endl;
	}


	return 0;
}



int print( vector<vector<double> > p)
{
	for(int i=0; i!=(p.size());++i){
		for(int j=0; j!=((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;}
		return 1;
}







int vectormultiplication(vector<vector<double> > a,vector<vector<double> > b,vector<vector<double> > &c)
{
	if((a[0]).size()!=b.size())
		cout<<"columns of 1st matrix and row of 2nd matrix does not match"<<endl;
	else{
		vector<double>temp;
		int row=a.size();
		int column=(b[1]).size();
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
	return 2;
}
