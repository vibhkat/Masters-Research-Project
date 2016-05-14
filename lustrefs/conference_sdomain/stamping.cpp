#include"stamping.h"
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
