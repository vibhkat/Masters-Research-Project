#include"stamping.h"
#include"functions.h"



int read_nodes(std::string filename, std::vector<double>& node_list)
{
	double temp_node;
	std::fstream fp(filename.c_str(),std::fstream::in);
	while(fp>>temp_node)
		node_list.push_back(temp_node);
	fp.close();

	return 0;
}

int populate_matrices(std::string filename,std::vector<std::vector<double> >& G_MATRIX, std::vector<std::vector<double> >& C_MATRIX, std::vector<std::vector<double> >& B_MATRIX,std::vector<std::vector<double> >& u_MATRIX,double points,double delta_t)
{
	std::cout<<"FILE NAME :"<<filename<<std::endl;
	std::string line,s;
	int column=0;
	std::ifstream myfile(filename.c_str());

	while(getline(myfile,line))
	{
		std::string word;
		netlist_t info;
		std::istringstream record(line);
		record>>word;
		info.name.push_back(word);
		while(record>>word)
		{
			double f;
			std::istringstream(word)>>f;
			info.values.push_back(f);   
		}


		/*


		//printing the name and values
		std::cout<<"string :"<<std::endl;
		for(int i=0;i<info.name.size();++i)
		{
		std::cout<<info.name[i]<<std::endl;
		}
		std::cout<<"values :"<<std::endl;
		for(int i=0;i<info.values.size();++i)
		{
		std::cout<<info.values[i]<<" ";
		}
		std::cout<<std::endl;

		 */
		addzeroes(G_MATRIX,C_MATRIX,B_MATRIX,info.name[0][0],info.values[0],info.values[1],column);

		if(info.name[0][0]=='G'|| info.name[0][0]=='g')
		{	      
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,info.values[0],info.values[1],info.values[2]);
		}
		else if(info.name[0][0]=='R' || info.name[0][0]=='r'){
			double newvalue=(1/info.values[2]);
			insertvaluesRG(G_MATRIX,C_MATRIX,B_MATRIX,info.values[0],info.values[1],newvalue);
		}
		else if(info.name[0][0]=='C' || info.name[0][0]=='c')
		{
			insertvaluesC(G_MATRIX,C_MATRIX,B_MATRIX,info.values[0],info.values[1],info.values[2]);
		}
		else if(info.name[0][0]=='L' || info.name[0][0]=='l')
		{
			insertvaluesL(G_MATRIX,C_MATRIX,B_MATRIX,info.values[0],info.values[1],info.values[2]);
		}
		else if(info.name[0][0]=='J' || info.name[0][0]=='j')
		{

			std::vector<double> time,magnitude;
			for(int i=2;i<info.values.size();++i)
			{
				if(i%2==0)
				{
					time.push_back(info.values[i]);
				}
				else
				{
					magnitude.push_back(info.values[i]);
				}
			}
			insertvaluesJ(B_MATRIX,u_MATRIX,info.values[0],info.values[1],time,magnitude,column,points,delta_t);


		}
		else if(info.name[0][0]=='V' || info.name[0][0]=='v')
		{
			std::vector<double> time,magnitude;
			for(int i=2;i<info.values.size();++i)
			{
				if(i%2==0)
				{
					time.push_back(info.values[i]);
				}
				else
				{
					magnitude.push_back(info.values[i]);
				}
			}

			insertvaluesV(G_MATRIX,C_MATRIX,B_MATRIX,u_MATRIX,info.values[0],info.values[1],time,magnitude,column,points,delta_t);
		}

		else if(info.name[0][0]=='Z' || info.name[0][0]=='z')
		{
			insertvaluesZ(G_MATRIX,C_MATRIX,info.values[0],info.values[1],info.values[2],info.values[3],info.values[4]);
		}
		else if(info.name[0][0]=='E' || info.name[0][0]=='e')
		{
			insertvaluesE(G_MATRIX,C_MATRIX,info.values[0],info.values[1],info.values[2],info.values[3],info.values[4]);
		}


	}

	return 2;
}








/* INSERTING ZEROES IN G,C,B MATRIC*/
int addzeroes(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,std::vector<std::vector<double> > &B,char S, double nodeI,double nodeJ,int &count)
{
	int size=(nodeI>nodeJ)?nodeI:nodeJ;
	int size1=G.size();//or size1=C.size,had to add this syntax because when the vectros were not empty,it was not able to take the size of the std::vectors.	
	int size2=B.size();	
	if(S=='R'||S=='r'||S=='G'||S=='g'||S=='C'||S=='c')
	{
		if(G.empty() && C.empty() && B.empty()  )
		{
			G=std::vector<std::vector<double> >(size,std::vector<double>(size,0));
			C=std::vector<std::vector<double> >(size,std::vector<double>(size,0));
			B=std::vector<std::vector<double> >(size,std::vector<double>(1,0));}
		else if(!G.empty() && !C.empty() && !B.empty() && size>G.size() && size>C.size())
		{
			for(int i=0; i<size1;++i)
			{
				for(int j=0; j<(size-(G.size())) ;++j)
				{
					G[i].push_back(0);
					C[i].push_back(0);
				}
			}
			for(int k=0; k<(size-size1) ; ++k){
				G.push_back(std::vector<double>(size,0));
				C.push_back(std::vector<double>(size,0));
				B.push_back(std::vector<double>(1,0));
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
		G.push_back(std::vector<double>((G.size()+1),0));
		C.push_back(std::vector<double>((C.size()+1),0));
		B.push_back(std::vector<double>(1,0));
	}

	else if(S=='J'||S=='j')
	{
		int s=0;
		for(int i=0;i<B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size()){
			++count;
			for(int i=0;i<B.size();++i)
			{
				B[i].push_back(0);
			}
		}
	}
	else if(S=='V'||S=='v')
	{
		for(int i=0;i<(G.size());++i)
		{
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(std::vector<double>((G.size()+1),0));
		C.push_back(std::vector<double>((C.size()+1),0));
		B.push_back(std::vector<double>((count+1),0));
		int s=0;
		for(int i=0;i<B.size();++i){
			if(B[i][count]==0)
				++s;
		}
		if(s<B.size())
		{
			++count;
			for(int i=0;i<B.size();++i)
			{
				B[i].push_back(0);
			}

		}
	}

	else if(S=='E'||S=='e')
	{
		for(int i=0;i<(G.size());++i)
		{
			G[i].push_back(0);
			C[i].push_back(0);
		}
		G.push_back(std::vector<double>((G.size()+1),0));
		C.push_back(std::vector<double>((C.size()+1),0));
		B.push_back(std::vector<double>((count+1),0));
	}


	return 3;
}
int insertvaluesRG(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,std::vector<std::vector<double> > &B,double nodeI,double nodeJ,double value)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i<(G.size());++i){
		for(int j=0;j<((G[i]).size());++j){
			if(((i==I)&&(j==J))||((i==J)&&(j==I))&& I>=0 && J>=0)
				G[i][j]=G[i][j]-value;
			else if(((i==I)&&(j==I)&& I>=0)||((i==J)&&(j==J)&& J>=0))
				G[i][j]=G[i][j]+value;
		}
	}
	return 4;
}

int insertvaluesC(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,std::vector<std::vector<double> > &B,double nodeI,double nodeJ,double value)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i<(C.size());++i){
		for(int j=0;j<((C[i]).size());++j){
			if(((i==I)&&(j==J))||((i==J)&&(j==I))&& I>=0 && J>=0)
				C[i][j]=C[i][j]-value;
			else if(((i==I)&&(j==I)&& I>=0)||((i==J)&&(j==J)&& J>=0))
				C[i][j]=C[i][j]+value;
		}
	}
	return 5;
}

int insertvaluesL(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,std::vector<std::vector<double> > &B,double nodeI,double nodeJ,double value)
{
	C[(C.size())-1][(C.size())-1]=C[(C.size())-1][(C.size())-1]+value;
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i<(G.size());++i){
		for(int j=0;j<((G[i]).size());++j){
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
	return 6;
}


int insertvaluesE(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
{
	int last=(G.size())-1;
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double NCplus=nodeK-1;
	double NCminus=nodeL-1;
	for(int i=0;i<(G.size());++i){
		for(int j=0;j<((G[i]).size());++j){
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
	return 8;
}

int insertvaluesZ(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,double nodeI,double nodeJ,double nodeK,double nodeL,double value)
{
	double Nplus=nodeI-1;
	double Nminus=nodeJ-1;
	double NCplus=nodeK-1;
	double NCminus=nodeL-1;
	for(int i=0;i<(G.size());++i){
		for(int j=0;j<((G[i]).size());++j){
			if((i==Nplus && j==NCplus && Nplus>=0 && NCplus>=0)||(i==Nminus && j==NCminus && Nminus>=0 && NCminus>=0))
			{
				C[i][j]=C[i][j]+value;
			}
			else if((i==Nplus && j==NCminus && Nplus>=0 && NCminus>=0)||(i==Nminus && j==NCplus && Nminus>=0 && NCplus>=0))
			{       C[i][j]=C[i][j]-value;

			}

		}
	}
	return 9;
}

int insertvaluesV(std::vector<std::vector<double> > &G,std::vector<std::vector<double> > &C,std::vector<std::vector<double> > &B,std::vector<std::vector<double> > &U,double nodeI,double nodeJ,std::vector<double> time,std::vector<double> magnitude,int &column, double points, double delta_t)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i<G.size();++i)
	{
		for(int j=0;j<((G[i]).size());++j)
		{
			if(i==(G.size()-1) && j==I && I>=0)
			{
				G[i][j]=G[i][j]-1;}
			else if(i==(G.size()-1) && j==J && J>=0)
			{
				G[i][j]=G[i][j]+1;
			}				
			else if(j==(G.size()-1) && i==I && I>=0)
			{
				G[i][j]=G[i][j]+1;
			}
			else if(j==(G.size()-1) && i==J && J>=0)
			{
				G[i][j]=G[i][j]-1;

			}
		}
	}
	B[(B.size()-1)][column]=B[(B.size()-1)][column]-1;

	std::vector<double> u_temp;
	int n=0;
	while(n<=points)
	{
		double t=n*delta_t;
		for(int i=1;i<time.size();++i)
		{
			if(t>=time[i-1] && t<=time[i])
			{
				u_temp.push_back( (((magnitude[i]-magnitude[i-1])/(time[i]-time[i-1]))*(t-time[i-1]))+magnitude[i-1]);

			}
		}
		++n;
	}
	U.push_back(u_temp);
	u_temp.clear();
	return 10;

}


int insertvaluesJ(std::vector<std::vector<double> > &B,std::vector<std::vector<double> > &U,double nodeI,double nodeJ,std::vector<double> time,std::vector<double> magnitude,int &count, double points, double delta_t)
{
	double I=nodeI-1;
	double J=nodeJ-1;
	for(int i=0;i<B.size();++i)
	{
		if(i==I && I>=0)
			B[i][count]=B[i][count]+1;
		else if(i==J && J>=0)
			B[i][count]=B[i][count]-1;
	}

	std::vector<double> u_temp;
	int n=0;
	while(n<=points)
	{
		double t=n*delta_t;
		for(int i=1;i<time.size();++i)
		{
			if(t>=time[i-1] && t<=time[i])
			{
				u_temp.push_back( (((magnitude[i]-magnitude[i-1])/(time[i]-time[i-1]))*(t-time[i-1]))+magnitude[i-1]);

			}
		}
		++n;
	}
	U.push_back(u_temp);
	u_temp.clear();

	return 7;
}

