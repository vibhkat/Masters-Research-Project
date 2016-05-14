#include <iostream>
#include <fstream>
#include<vector>
#include<cmath>
using namespace std;
#define tr 0.1e-9
#define n 1
#define l 15
//#define m 250
void print( vector<vector<double> > &);

int main () {
	cout<<"enter the number of transmission lines"<<endl;
	size_t number;
	cin>>number;
	for(size_t no=0;no!=number;++no)
	{
		if(n==1){
			fstream capacitance,resistance,inductance,admittance;
			capacitance.open("C_1.txt",fstream::in);
			resistance.open("R_1.txt",fstream::in);
			inductance.open("L_1.txt",fstream::in);
			admittance.open("G_1.txt",fstream::in);
			double C,R,L,G;
			double cap_value,res_value,ind_value,adm_value;
			if(capacitance.is_open() && resistance.is_open() && inductance.is_open() && admittance.is_open())
			{
				capacitance>>cap_value;
				C=cap_value/100;
				resistance>>res_value;
				R=res_value/100;
				inductance>>ind_value;
				L=ind_value/100;
				admittance>>adm_value;
				G=adm_value/100;		


			}
			else cout<<"Unable to open file"<<endl;

			cout<<"C value: "<<C<<endl;
			cout<<"R value: "<<R<<endl;
			cout<<"L value: "<<L<<endl;
			cout<<"G value: "<<G<<endl;


//			finding the number of sections
			double m=ceil(20*sqrt(L*C)*l/tr);
			cout<<"the sections are :"<<m<<endl;
cout<<" the Length is : "<<l;
			// storing the nodes in a vector
			vector<vector<double> > transmission_line;
			vector<double> oneD;
			size_t s,start_node;
			cout<<"enter the starting node for line number "<<no+1<<endl;
			cin>>start_node;
			s=start_node-1;
			for(size_t section=0;section!=m;++section)
			{
				if(section==0)
				{
					for(size_t i=0;i!=3;++i)
						oneD.push_back(++s);
					transmission_line.push_back(oneD);
					oneD.clear();
				}
				else
				{
					oneD.push_back(transmission_line[section-1][2]);
					for(size_t i=0;i!=2;++i)
						oneD.push_back(++s);
					transmission_line.push_back(oneD);
					oneD.clear();
				}

			}

			cout<<"the end terminal is :"<<s<<endl;


			//printing the nodes
			for(size_t i=0;i!=transmission_line.size();++i)
			{
				for(size_t j=0;j!=(transmission_line[i]).size();++j)
					cout<<transmission_line[i][j]<<" ";
				cout<<endl;
			}

			//writng the input and output nodes in a file
			fstream input("input_terminal.txt",fstream::out);
			cout<<"the input terminals are :"<<endl;
			cout<<transmission_line[0][0]<<endl;
			input<<transmission_line[0][0]<<endl;
			input.close();
			fstream output("output_terminal.txt",fstream::out);
			cout<<"the output terminals are :"<<endl;
			cout<<transmission_line[m-1][2]<<endl;
			output<<transmission_line[m-1][2]<<endl;
			output.close();
			//writing the spice input file that is to be read by solver
			fstream myfile;

			if(no==0)
				myfile.open("time_test_15.txt",fstream::out);
			else
				myfile.open("time_test_15.txt",fstream::out|fstream::app);

			if(myfile.is_open())
			{
				cout<<"file is being written"<<endl;
				for(size_t section=0;section!=transmission_line.size();++section)
					myfile<<"R"<<" "<<transmission_line[section][0]<<" "<<transmission_line[section][1]<<" "<<(R*l/m)<<endl;
				for(size_t section=0;section!=transmission_line.size();++section)
					myfile<<"G"<<" "<<transmission_line[section][2]<<" "<<0<<" "<<(G*l/m)<<endl;
				for(size_t section=0;section!=transmission_line.size();++section)
					myfile<<"C"<<" "<<transmission_line[section][2]<<" "<<0<<" "<<(C*l/m)<<endl;
				for(size_t section=0;section!=transmission_line.size();++section)
					myfile<<"L"<<" "<<transmission_line[section][1]<<" "<<transmission_line[section][2]<<" "<<(L*l/m)<<endl;
				cout<<"file is written"<<endl;
				myfile.close();
			}
			else cout<<"unable to open file"<<endl;
			transmission_line.clear();



		}
	}
	return 0;
}




void print( vector<vector<double> > &p)
{
	for(size_t i=0; i!=(p.size());++i){
		for(size_t j=0; j!=((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;}
}


