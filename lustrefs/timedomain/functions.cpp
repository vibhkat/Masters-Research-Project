#include"functions.h"

int print_2D(std::vector<std::vector<double> > p)
{
	for(int i=0; i<(p.size());++i)
	{
		for(int j=0; j<((p[i]).size());++j)
			std::cout<<p[i][j]<<" ";
		std::cout<<";"<<std::endl;
	}
	return 3;
}

int print_1D(std::vector<double> p)
{
	for(int i=0;i<p.size();++i)
	{
		std::cout<<p[i]<<" ";
	}
	std::cout<<std::endl;

	return 4;
}

int vectoraddition(std::vector<std::vector<double> > A,std::vector<std::vector<double> > B,std::vector<std::vector<double> > &C)
{

	if(A.size()==B.size() && (A.at(0)).size()==(B.at(0)).size())
	{
		std::vector<double> temp;
		for(int i=0;i<A.size();++i)
		{
			for(int j=0;j<(A.at(0)).size();++j)
			{
				temp.push_back(A[i][j]+B[i][j]);
			}
			C.push_back(temp);
			temp.clear();
		}
	}
	else std::cout<<" The dimension of both Matrices for addtion do not match "<<std::endl;



	return 5;
}

int vectormultiplication(std::vector<std::vector<double> > a,std::vector<std::vector<double> > b,std::vector<std::vector<double> > &c)
{
	if(((a.at(0)).size())!=b.size()){
		std::cout<<"columns of 1st matrix and row of 2nd matrix does not match"<<std::endl;
		std::cout<<((a.at(0)).size())<<";;;;;;;;;;;;;;;;;;;"<<b.size()<<std::endl;

	}


	else{
		std::vector<double>temp;
		int row=a.size();
		int column=(b[0]).size();
		int inner=b.size();
		for(int i=0;i<row;++i){
			for(int j=0;j<column;++j){
				double sum=0;
				for(int k=0;k<inner;++k){
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

int AinverseB(std::vector<std::vector<double> > Ad,std::vector<std::vector<double> > Bd,std::vector<std::vector<double> > &Cd)
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
	for(int i=0;i<(Ad[0]).size();++i)
		for(int j=0;j<Ad.size();++j)
			ad[kd++]=Ad[j][i];
	for(int i=0;i<(Bd[0]).size();++i)
		for(int j=0;j<Bd.size();++j)
			bd[ld++]=Bd[j][i];
	dgesv_( &Nd, &NRHSd, ad, &LDAd, IPIVd, bd, &LDBd, &INFOd );
	if( INFOd > 0 ) {
		std::cout<<"The diagonal element of the triangular factor of A,\n";
		std::cout<< "U("<<INFOd<<","<<INFOd<<" is jero, so that A is singular;\n";
		std::cout<<( "the solution could not be computed.\n" );
		//	return(1);
	}

	//	cout << "===============================IN AINVERSEB 1" << endl;
	std::vector<double>tempd;
	for(int i=0;i<Nd;++i)
	{
		for(int j=0;j<NRHSd;++j)
			tempd.push_back(bd[i+j*LDBd]);
		Cd.push_back(tempd);
		tempd.clear();
	}
	//	cout << "===============================IN AINVERSEB 2" << endl;
	return 15;
}


int original_solver(std::string filename,std::vector<std::vector<double> > G_MATRIX,std::vector<std::vector<double> > C_MATRIX,std::vector<std::vector<double> > B_MATRIX,std::vector<std::vector<double> > u_MATRIX, double delta_t, int points,std::vector<double> input_node, std::vector<double> output_node)
{
	std::vector<double> temp;
	std::vector<std::vector<double> > A_plus,A_minus;

	//computing A_plus
	for(int i=0;i<C_MATRIX.size();++i)
	{
		for(int j=0;j<(C_MATRIX.at(0)).size();++j)
		{
			temp.push_back((2*C_MATRIX[i][j]/delta_t)+G_MATRIX[i][j]);
		}
		A_plus.push_back(temp);
		temp.clear();
	}


	//computing A_minus
	for(int i=0;i<C_MATRIX.size();++i)
	{
		for(int j=0;j<(C_MATRIX.at(0)).size();++j)
		{
			temp.push_back((2*C_MATRIX[i][j]/delta_t)-G_MATRIX[i][j]);
		}
		A_minus.push_back(temp);
		temp.clear();
	}
	std::vector<std::vector<double> > X_new,X_previous;
	std::fstream write(filename.c_str(),std::fstream::out);
	int n=0;
	while(n<=points)
	{
		std::cout<<" ORIGINAL POINT : "<<n<<std::endl;
		write<<n*delta_t<<'\t';

		if(n==0)
		{
			for(int i=0;i<G_MATRIX.size();++i)
			{
				X_new.push_back(std::vector<double>(1,0));
			}
		}
		else
		{
			std::vector<std::vector<double> > A_minus_X_previous,Bu,u_temp,B;
			for(int i=0;i<u_MATRIX.size();++i)
			{
				double value=u_MATRIX[i][n] + u_MATRIX[i][n-1];
				u_temp.push_back(std::vector<double>(1,value));
			}
			vectormultiplication(B_MATRIX,u_temp,Bu);
			vectormultiplication(A_minus,X_previous,A_minus_X_previous);
			vectoraddition(A_minus_X_previous,Bu,B);
			AinverseB(A_plus,B,X_new);
			u_temp.clear();
			Bu.clear();
			B.clear();
			A_minus_X_previous.clear();
		}
		for(int i=0;i<X_new.size();++i)
		{
			for(int k=0;k<output_node.size();++k)
			{
				if(i==(output_node[k]-1))
					write<<X_new[i][0]<<'\t';
			}
		}
		write<<';'<<std::endl;
		X_previous.clear();
		X_previous=X_new;
		X_new.clear();

		++n;
	}

	return 10;
}


int QR_factorization(std::vector<std::vector<double> > &a,std::vector<std::vector<double> > &E)
{
	std::vector<double> u;
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
				E.push_back(std::vector<double>(1,value));
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

int vector_transpose(std::vector<std::vector<double> > a,std::vector<std::vector<double> > &b)
{
	std::vector<double>temp;
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


int blockArnoldi(std::vector<std::vector<double> > G,std::vector<std::vector<double> > C,std::vector<std::vector<double> > B,std::vector<std::vector<double> > &X,int q_value)
{
	std::vector<double> insert,range;
	std::vector<std::vector<double> > R,Q,V,X_temp,Z,X_transpose,H,XH;
	AinverseB(G,B,R);
	QR_factorization(R,Q);
	std::cout<<"the number of columns in one moment are :"<<(Q.at(0)).size()<<" and row are :"<<Q.size()<<std::endl;	
	int q,n,input1,input2;
	X=Q;
	std::cout<<"done computing X0"<<std::endl;
	Q.clear();
	R.clear();
	range.push_back((X[0]).size());
	std::cout<<"Enter the value of q"<<std::endl;
	//	cin>>q;
	q=q_value;
	std::cout<<"the value of q is "<<q<<std::endl;		
	if((q%((B[0]).size()))==0)
		n=q/((B[0]).size());
	else
		n=floor(q/((B[0]).size()))+1;
	std::cout<<"the value of n is :"<<n<<std::endl;	
	for(int kA=1;kA<=n;++kA)
	{
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
		vectormultiplication(C,X_temp,V);
		X_temp.clear();
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
		QR_factorization(Z,Q);
		for(int iAB=0;iAB!=Q.size();++iAB)
			for(int jAB=0;jAB!=(Q[0]).size();++jAB)
				(X[iAB]).push_back(Q[iAB][jAB]);
		range.push_back((X[0]).size());
		std::cout<<"printing the X"<<kA<<std::endl;
		Q.clear();
		Z.clear();
	}
	return 19;
}



int reduced_MATRIX_generator(std::vector<std::vector<double> > G_MATRIX,std::vector<std::vector<double> > C_MATRIX,std::vector<std::vector<double> > B_MATRIX,std::vector<std::vector<double> > X_MATRIX,std::vector<std::vector<double> > &G_cap,std::vector<std::vector<double> > &C_cap,std::vector<std::vector<double> > &B_cap)
{
	std::vector<std::vector<double> > X_transpose,multiple;
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
	return 20;
}



int reduced_solver(std::string filename,std::vector<std::vector<double> > G_cap,std::vector<std::vector<double> > C_cap,std::vector<std::vector<double> > B_cap,std::vector<std::vector<double> > u_MATRIX, double delta_t, int points,std::vector<double> input_node, std::vector<double> output_node, std::vector<std::vector<double> > X_MATRIX)
{
	std::vector<double> temp;
	std::vector<std::vector<double> > A_plus,A_minus;

	//computing A_plus
	for(int i=0;i<C_cap.size();++i)
	{
		for(int j=0;j<(C_cap.at(0)).size();++j)
		{
			temp.push_back((2*C_cap[i][j]/delta_t)+G_cap[i][j]);
		}
		A_plus.push_back(temp);
		temp.clear();
	}


	//computing A_minus
	for(int i=0;i<C_cap.size();++i)
	{
		for(int j=0;j<(C_cap.at(0)).size();++j)
		{
			temp.push_back((2*C_cap[i][j]/delta_t)-G_cap[i][j]);
		}
		A_minus.push_back(temp);
		temp.clear();
	}
	std::vector<std::vector<double> > X_new,X_previous,X_new_original;
	std::fstream write(filename.c_str(),std::fstream::out);
	int n=0;
	while(n<=points)
	{
		std::cout<<" MOR POINT : "<<n<<std::endl;
		write<<n*delta_t<<'\t';

		if(n==0)
		{
			for(int i=0;i<G_cap.size();++i)
			{
				X_new.push_back(std::vector<double>(1,0));
			}
		}
		else
		{
			std::vector<std::vector<double> > A_minus_X_previous,Bu,u_temp,B;
			for(int i=0;i<u_MATRIX.size();++i)
			{
				double value=u_MATRIX[i][n] + u_MATRIX[i][n-1];
				u_temp.push_back(std::vector<double>(1,value));
			}
			vectormultiplication(B_cap,u_temp,Bu);
			vectormultiplication(A_minus,X_previous,A_minus_X_previous);
			vectoraddition(A_minus_X_previous,Bu,B);
			AinverseB(A_plus,B,X_new);
			u_temp.clear();
			Bu.clear();
			B.clear();
			A_minus_X_previous.clear();
		}
		vectormultiplication(X_MATRIX,X_new,X_new_original);		
		for(int i=0;i<X_new_original.size();++i)
		{
			for(int k=0;k<output_node.size();++k)
			{
				if(i==(output_node[k]-1))
					write<<X_new_original[i][0]<<'\t';
			}
		}
		write<<';'<<std::endl;
		X_previous.clear();
		X_previous=X_new;
		X_new.clear();
		X_new_original.clear();
		++n;
	}

	return 100;
}


