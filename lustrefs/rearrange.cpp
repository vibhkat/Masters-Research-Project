#include<iostream>
#include<fstream>
#include<string>
#include"stdio.h"
using namespace std;
int main()
{
	string line,filename;
cout<<" enter the file to be rearraged : ";
cin>>filename;
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

	fstream write(filename.c_str(),fstream::out);
	fstream read("temporary.txt",fstream::in);
	while(getline(read,line))
		write<<line<<endl;

	write.close();
	read.close();

	if( remove( "temporary.txt" ) != 0 )
		perror( "Error deleting file" );
	else
		puts( "File successfully deleted" );

	return 0;
}
