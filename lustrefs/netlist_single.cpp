#include <iostream>
#include <fstream>
#include<vector>
#include<cmath>
#include<fstream>
#include<cstring>
#include<sstream>
using namespace std;
#define tr 0.1e-9
#define n 1
#define l 15
#define m 250
void print( vector<vector<double> > &);

int main (int agrc,char *argv[] ) 
{
//string cap_file=argv[1],res_file=argv[2],induc_file=argv[3],admit_file=argv[4],netlist_file=argv[5];
int points;
cout<<"Enter the number of points:";
cin>>points;

for( int pts=1;pts<=points;++pts)
{
string s_pts;
ostringstream ss;
ss<<pts;
s_pts=ss.str();
string cap_file="C_"+s_pts+".txt";
string res_file="R_"+s_pts+".txt";
string induc_file="L_"+s_pts+".txt";
string admit_file="G_"+s_pts+".txt";
string netlist_file="time_test"+s_pts+".txt";
cout<<cap_file<<" "<<res_file<<" "<<induc_file<<" "<<admit_file<<" "<<netlist_file<<endl;


int number=6;
	cout<<"enter the number of transmission lines :"<<number<<endl;
//	int number;
//	cin>>number;
int in_node[6]={1,502,1002,1503,2004,2505};
	for(int no=0;no<number;++no)
	{
		if(n==1){
			fstream capacitance,resistance,inductance,admittance;
			capacitance.open(cap_file.c_str(),fstream::in);
			resistance.open(res_file.c_str(),fstream::in);
			inductance.open(induc_file.c_str(),fstream::in);
			admittance.open(admit_file.c_str(),fstream::in);
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
			else cout<<"Unable to open file.........................."<<endl;

			cout<<"C value: "<<C<<endl;
			cout<<"R value: "<<R<<endl;
			cout<<"L value: "<<L<<endl;
			cout<<"G value: "<<G<<endl;


			//finding the number of sections
			//			double m=ceil(20*sqrt(L*C)*l/tr);
			cout<<"the sections are :"<<m<<endl;
			cout<<" the Length is : "<<l<<endl;
			// storing the nodes in a vector
			vector<vector<double> > transmission_line;
			vector<double> oneD;
			int s,start_node;
			cout<<"enter the starting node for line number "<<no+1<<endl;
//			cin>>start_node;
start_node=in_node[no];
cout<<start_node<<endl;
			s=start_node-1;
			for(int section=0;section<m;++section)
			{
				if(section==0)
				{
					for(int i=0;i<3;++i)
						oneD.push_back(++s);
					transmission_line.push_back(oneD);
					oneD.clear();
				}
				else
				{
					oneD.push_back(transmission_line[section-1][2]);
					for(int i=0;i<2;++i)
						oneD.push_back(++s);
					transmission_line.push_back(oneD);
					oneD.clear();
				}

			}

			cout<<"the end terminal is :"<<s<<endl;

/*
			//printing the nodes
			for(int i=0;i<transmission_line.size();++i)
			{
				for(int j=0;j<(transmission_line[i]).size();++j)
					cout<<transmission_line[i][j]<<" ";
				cout<<endl;
			}

*/			//writng the input and output nodes in a file
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
				myfile.open(netlist_file.c_str(),fstream::out);
			else
				myfile.open(netlist_file.c_str(),fstream::out|fstream::app);

			if(myfile.is_open())
			{
				cout<<"file is being written"<<endl;
				for(int section=0;section<transmission_line.size();++section)
					myfile<<"R"<<" "<<transmission_line[section][0]<<" "<<transmission_line[section][1]<<" "<<(R*l/m)<<endl;
				for(int section=0;section<transmission_line.size();++section)
					myfile<<"G"<<" "<<transmission_line[section][2]<<" "<<0<<" "<<(G*l/m)<<endl;
				for(int section=0;section<transmission_line.size();++section)
					myfile<<"C"<<" "<<transmission_line[section][2]<<" "<<0<<" "<<(C*l/m)<<endl;
				for(int section=0;section<transmission_line.size();++section)
					myfile<<"L"<<" "<<transmission_line[section][1]<<" "<<transmission_line[section][2]<<" "<<(L*l/m)<<endl;
				cout<<"file is written"<<endl;
				myfile.close();
			}
			else cout<<"unable to open file*********************************************"<<endl;
			transmission_line.clear();



		}

	}


cout<<"inserting the extras in input file"<<endl;
 fstream read("extras.txt",fstream::in);// read other RLGC components in circuit
                if(read.is_open())
                {
                        fstream write(netlist_file.c_str(),fstream::out|fstream::app);
                        if(write.is_open())
                        {
                                string sentence;
                                while(getline(read,sentence))
                                {
                                        write<<sentence<<endl;
                                }
                                write.close();
                        }
                        else cout<<"unable to open file----------------------------------------"<<endl;

                        read.close();
                }

                else cout<<"unable to open file"<<endl;
cout<<"done inserting the extra in input file"<<endl;	
cout<<"rearranging the file"<<endl;
string line,filename=netlist_file;
cout<<" enter the file to be rearraged : ";
	fstream myfile1_1(filename.c_str(),fstream::in);
	fstream myfile2("temporary.txt",fstream::out|fstream::app);
	if(myfile1_1.is_open())
	{
		cout<<"reading and writng resistance"<<endl;
		while(getline(myfile1_1,line))
		{
			if(line[0]=='R'||line[0]=='r')
				myfile2<<line<<endl;
		}
		cout<<"done reading and writing the resistance"<<endl;
	}
	else cout<<"unable to open the file"<<endl;
	myfile1_1.close();


	fstream myfile1_2(filename.c_str(),fstream::in);
	if(myfile1_2.is_open())
	{
		cout<<"reading and writng admittance"<<endl;
		while(getline(myfile1_2,line))
		{
			if(line[0]=='G'||line[0]=='g')
				myfile2<<line<<endl;
		}
		myfile1_2.close();
		cout<<"done reading and writing the admittance"<<endl;
	}
	else cout<<"unable to open the file"<<endl;


	fstream myfile1_3(filename.c_str(),fstream::in);
	if(myfile1_3.is_open())
	{
		cout<<"reading and writng capacitance"<<endl;
		while(getline(myfile1_3,line))
		{
			if(line[0]=='C'||line[0]=='c')
				myfile2<<line<<endl;
		}
		myfile1_3.close();
		cout<<"done reading and writing the capacitance"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	fstream myfile1_4(filename.c_str(),fstream::in);
	if(myfile1_4.is_open())
	{
		cout<<"reading and writng current source"<<endl;
		while(getline(myfile1_4,line))
		{
			if(line[0]=='J'||line[0]=='j')
				myfile2<<line<<endl;
		}
		myfile1_4.close();
		cout<<"done reading and writing current source"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	fstream myfile1_5(filename.c_str(),fstream::in);
	if(myfile1_5.is_open())
	{
		cout<<"reading and writng inductance"<<endl;
		while(getline(myfile1_5,line))
		{
			if(line[0]=='L'||line[0]=='l')
				myfile2<<line<<endl;
		}
		myfile1_5.close();
		cout<<"done reading and writing inductance"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	fstream myfile1_6(filename.c_str(),fstream::in);
	if(myfile1_6.is_open())
	{
		cout<<"reading and writng VCCS"<<endl;
		while(getline(myfile1_6,line))
		{
			if(line[0]=='Z'||line[0]=='z')
				myfile2<<line<<endl;
		}
		myfile1_6.close();
		cout<<"done reading and VCCS"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	fstream myfile1_7(filename.c_str(),fstream::in);
	if(myfile1_7.is_open())
	{
		cout<<"reading and writng VCVS"<<endl;
		while(getline(myfile1_7,line))
		{
			if(line[0]=='E'||line[0]=='e')
				myfile2<<line<<endl;
		}
		myfile1_7.close();
		cout<<"done reading and writing CCVS"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	fstream myfile1_8(filename.c_str(),fstream::in);
	if(myfile1_8.is_open())
	{
		cout<<"reading and writng voltage source"<<endl;
		while(getline(myfile1_8,line))
		{
			if(line[0]=='V'||line[0]=='v')
				myfile2<<line<<endl;
		}
		myfile1_8.close();
		cout<<"done reading and writing inductance"<<endl;
	}
	else cout<<"unable to open the file"<<endl;

	myfile2.close();

	fstream write1(filename.c_str(),fstream::out);
	fstream read1("temporary.txt",fstream::in);
	while(getline(read1,line))
		write1<<line<<endl;

	write1.close();
	read1.close();

	if( remove( "temporary.txt" ) != 0 )
		perror( "Error deleting file" );
	else
		puts( "File successfully deleted" );


}
return 0;
}




void print( vector<vector<double> > &p)
{
	for(int i=0; i<(p.size());++i){
		for(int j=0; j!=((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;}
}


