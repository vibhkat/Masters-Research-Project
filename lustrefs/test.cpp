#include<iostream>
#include<vector>
using namespace std;
int main()
{
	char b='y';
	vector<double> x;
while(b=='y')
{
char a='y';
vector<double>y;
	while(a=='y')
	{
//		vector<double> y;
		double input;
		cout<<"enter the number: ";
		cin>>input;
		y.push_back(input);
		cout<<" do you want to continue(y/n):";
		cin>>a;
		if(a=='n')
		{
		cout<<"in --------------------"<<endl;
if(x.empty())
		x=y;
else
{
for(int i=0;i<y.size();++i)
x.push_back(y[i]);
}
                cout<<"out======================"<<endl;
		}
	}
cout<<"dou you want to start again(y/n):";
cin>>b;
}
cout<<"the size of X is :"<<x.size()<<endl;
	for(int i=0;i<x.size();++i)
	cout<<x[i]<<" ";
	cout<<endl;
	return 0;
}
