#include"functions.h"
#include"stamping.h"

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



int original_plus_MOR_giving_X_big(vector<vector<double> > &X_big,double frequency_max,double frequency_min, double points,string input_file)
{	
double interval;
	interval=(frequency_max-frequency_min)/(points-1);
	string in_1,in_2,in_3;
	cout<<"ENTER THE MAIN NETLIST FILE : ";
//	cin>>in_1;
in_1=input_file;
//	cout<<"enter the file containing  INPUT TERMINAL : " ;
//	cin>>in_2;
//	cout<<"enter the file containing OUTPUT TERMINALS : ";
//	cin>>in_3;
	cout<<"1) MAIN NETLIST FILE : "<<in_1<<endl;//<<"2) INPUT TERMINAL FILE : "<<in_2<<endl<<"3) OUTPUT TERMINAL FILE : "<<in_3<<endl;

	cout<<" DID U CHANGE THE FREQUENCY"<<endl/*<<"REMEMBER ARRAY IS A BITCH"<<endl*/<<"First enter RGC irrespective of the order,then enter L,then independent sources and then dependent sources"<<endl<<"R or r for resistance"<<endl<<"G or g for admittance"<<endl<<"C or c for capacitance"<<endl<<"L or l for inductance"<<endl<<"J or j for current source"<<endl<<"E or e for VCVS"<<endl<<"Z or z for VCCS"<<endl<<"H or h for CCVS"<<endl<<"V or v for independent voltage source"<<endl<<"PRESS ENTER TWICE AFTER DONE ENTERING THE CIRCUIT INFORAMTION IN THE FORM OF SPICE INPUT"<<endl;
	vector<complex<double> > u_MATRIX;
	vector<vector<double> >G_MATRIX;
	vector<vector<double> >C_MATRIX;
	vector<vector<double> >B_MATRIX;
	G_MATRIX.clear();
	C_MATRIX.clear();
	B_MATRIX.clear();


	//reading the input and output node from files created from pul.c code        
	vector<double> input_node,output_node;
	read_nodes("input_terminal.txt", input_node);
	read_nodes("output_terminal.txt", output_node);
	//printing the nodes
	for(int i=0;i<input_node.size();++i)
		cout<<input_node[i]<<" ";
	cout<<endl;
	for(int i=0;i<output_node.size();++i)
		cout<<output_node[i]<<" ";
	cout<<endl;
 double t1=clock();
	populate_matrices(in_1, G_MATRIX, C_MATRIX, B_MATRIX,u_MATRIX);
string C_matrix="SERIAL_C_"+in_1;
string G_matrix="SERIAL_G_"+in_1;
file_write(C_matrix,C_MATRIX);
file_write(G_matrix,G_MATRIX);
cout<<" C_MATRIX AND G_MATRIX written in "<<C_matrix<<"  and "<<G_matrix<<endl;

double t2=clock();
 cout<<" TIME taken to populate G,C and B matrices is : "<<(t2-t1)/double(CLOCKS_PER_SEC) << "seconds" <<endl;

	//	cout<<"G MATRIX"<<endl;	
	//	print(G_MATRIX);
	//	cout<<endl;
	//	cout<<"C MATRIX"<<endl;
	//	print(C_MATRIX);
	cout<<"THE NUMBER OF VARIABLE ARE :"<<G_MATRIX.size()<<endl;
	//	cout<<endl<<"B MATRIX"<<endl;
	//	print(B_MATRIX);
	//	cout<<endl<<"u complex vector"<<endl;
	//	printcomplexvector(u_MATRIX);
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(int i=0;i<(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	//printing tht ucomplexaray
	//	cout<<endl<<"ucomplexarray"<<endl;
	//	for(int i=0;i!=(u_MATRIX).size();++i)
	//		cout<<uarraymatrix[i][0]<<endl;
	//converting the B VECTOR MATRIX to B ARRAY MATRIX
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(int i=0;i<(B_MATRIX).size();++i)
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

/*	//comment for original solving starts here	
	fstream result11("MATLAB_result_at_output_node_at_multiple_frequencies.txt",fstream::out);
	if(result11.is_open() )
	{
		double frequency;
		int timer1=0;
		while(timer1<points)
		{
			cout<<"original"<<timer1+1<<endl;	
			frequency=frequency_min+(interval*timer1);
			result11<<frequency<<'\t';

			//combining G_MATRIX and C_MATRIX into GplussC_MATRIX
			complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
			for(int i=0;i<(G_MATRIX).size();++i)
			{
				for(int j=0;j<(G_MATRIX).size();++j)
				{
					GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
				}
			}
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
			for(int i=0;i<n;++i)
			{
				for(int j=0;j<nrhs;++j)
				{
					original_solution[i][j]=b[i+j*ldb];
				}
			}
			for(int i=0;i<G_MATRIX.size();++i)
			{
				for(int k=0;k<output_node.size();++k)
				{
					if(i==(output_node[k]-1)){

						if((original_solution[i][0]).imag()>=0)
							result11<<(original_solution[i][0]).real()<<'+'<<(original_solution[i][0]).imag()<<'i'<<'\t';
						else
							result11<<(original_solution[i][0]).real()<<(original_solution[i][0]).imag()<<'i'<<'\t';

					}

				}
			}
			result11<<';'<<endl;

			++timer1;
		}//for loop of frequency ends here

	}
	else cout<<"unable to open file"<<endl;

	cout<<"done computing Ax=b of original solution for different frequencies and writing for MATLAB READING in MATLAB_result_at_output_node_at_multiple_frequencies.txt"<<endl;
*/	//comment for original model solving ends here






	//function to generate the X matrices containing the genrerated perpendicular columns for MOR
	cout<<endl<<"STARTING THE MOR FROM HERE"<<endl;	
	char answer='y';//while loop starts here
	while(answer=='y')
	{	

		int expansion=1;
cout<<"enter the number of expansion points :"<<expansion<<endl;

//		cout<<"enter the number of expansion points :";
	//	cin>>expansion;
		double abcd[expansion],efgh[expansion],expan_freq;
		int q_val;
		for(int i=0;i<expansion;++i)

		{
                expan_freq=0;
 cout<<"enter the "<<i+1<<" expansion frequency: "<<expan_freq<<endl;
	
//	cout<< enter the "<<i+1<<" expansion frequency: ";
//		cin>>expan_freq;

abcd[i]=expan_freq;


q_val=414;			
cout<<"enter the value of q : "<<q_val<<endl;
			


//			cout<<"enter the value of q : ";

//cin>>q_val;
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
string s1="SERIAL_X_MATRIX_"+in_1;
file_write(s1,X_MATRIX);
cout<<"written X_MATRIX in "<<s1<<endl;
		//blockArnoldi(G_MATRIX,C_MATRIX,B_MATRIX,X_MATRIX);
		//	cout<<"printing the X_MATIRIX"<<endl;
		//	print(X_MATRIX);
		//checking each coulumvnis orthnormal to each other		
/*		cout<<"done computing the orthonormal matrices using block Arnoldi algorithm"<<endl;
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
double t5=clock();

		fstream result22("MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt",fstream::out);
		if(result22.is_open() )
		{
			double frequency;
			int timer2=0;	
			while(timer2<points)
			{
				cout<<"computing MOR at "<<timer2+1<<" frequency point"<<endl;			
				frequency=frequency_min+(interval*timer2);
				result22<<frequency<<'\t';
				// Combining G_cap and C_cap into G_capplussC_cap
				complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
				for(int i=0;i<(G_cap).size();++i){
					for(int j=0;j<(G_cap).size();++j)
					{
						G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
					}
				}

				//copying elements of GcapplusCcap in a one dimension array acap
				complex<double> a_cap[(G_cap.size())*(G_cap.size())];
				int k_cap=0;
				for(int i=0;i<(G_cap).size();++i)
					for(int j=0;j<(G_cap).size();++j)
						a_cap[k_cap++]=G_capplussC_cap[j][i];

				//converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX
				complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
				for(int i=0;i<(B_cap).size();++i)
					for(int j=0;j<(u_MATRIX).size();++j)
						B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);

				//multiplying the B_COMPLEXARRAY AND uarraymatrix
				complex<double> B_capucomplexmultiple[(B_cap).size()][1];
				for(int i=0;i<(B_cap).size();++i)
					for(int j=0;j<(u_MATRIX).size();++j)
						B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];


				//copying elements of B_capucomplexmultiple into b_cap one dimension matrix
				complex<double> b_cap[(B_cap).size()];
				for(int i=0;i<(B_cap).size();++i)
					b_cap[i]=B_capucomplexmultiple[i][0];

				//computing a_capx=b_cap using cgesv routine

				int n_cap=G_cap.size();
				int nrhs_cap=1;

				int lda_cap=n_cap;
				int ldb_cap=n_cap;
				int info_cap;
				int ipiv_cap[n_cap];
				zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
				if( info_cap > 0 ) {
					printf( "The diagonal element of the triangular factor of A,\n" );
					printf( "U(%i,%i) is Jero, so that A is singular;\n", info_cap, info_cap );
					printf( "the solution could not be computed.\n" );
					return( 1 );
				}

				//converting X_matrix into X complex arary matrices

				complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
				for(int i=0;i<X_MATRIX.size();++i)
					for(int j=0;j<(X_MATRIX[0]).size();++j)
						X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

				//converting b_cap from inverse into x_cap 2 dimension array
				complex<double> x_cap[G_cap.size()][1];

				for(int i=0;i<n_cap;++i)
					for(int j=0;j<nrhs_cap;++j)
						x_cap[i][j]=b_cap[i+j*ldb_cap];

				complex<double>MORsolutioncompare[X_MATRIX.size()][1];
				for(int i=0;i<X_MATRIX.size();++i)
					for(int j=0;j<(X_MATRIX[0]).size();++j)
						MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];
				for(int i=0;i<G_MATRIX.size();++i)
				{
					for(int k=0;k<output_node.size();++k)
					{
						if(i==(output_node[k]-1)){
							if((MORsolutioncompare[i][0]).imag()>=0)              
								result22<<(MORsolutioncompare[i][0]).real()<<'+'<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';
							else
								result22<<(MORsolutioncompare[i][0]).real()<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';

						}
					}
				}
				result22<<';'<<endl;
				++timer2;			
			}
			result22.close();

		}
		else cout<<"unable top open file"<<endl;
double t6=clock();
cout<<" TIME taken to compute MOR from individual orhtonormal matrices is : "<<(t6-t5)/double(CLOCKS_PER_SEC) << "seconds"<<endl;

		cout<<"done computing Ax=b of MOR solution for different frequencies and writing for MATLAB READING in MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt "<<endl;		
		//cout<<"the size of Original solution is :"<<original_solution.size()<<" and the size of MOR solution compare is :"<<MORsolutioncompare.size()<<endl;
		//	out<<" the original and MOR solution written in result.out"<<endl;
*/
answer='n';	
	cout<<"do you want to continue for new value of Q (y/n)---- -:"<<answer<<endl;
//	cin>>answer;
		if(answer=='n')
		{
			if(X_big.empty())
				X_big=X_MATRIX;
			else
			{
				for(int i=0;i<X_MATRIX.size();++i)
				{
					for(int j=0;j<(X_MATRIX.at(0)).size();++j)
					{
						(X_big[i]).push_back(X_MATRIX[i][j]);
					}
				}
			}
		}

	}//whille loop


	return 30;
}

int SVD(vector<vector<double> > &X_svd,vector<double> &sigma,vector<vector<double> >X_big)
{
        int m=X_big.size();
        int n=(X_big.at(0)).size();
        int lda =m,ldu=m,ldvt=n,info,lwork;
        double wkopt;
        double* work;
        int sn;

        if(m>n) sn=n;
        else  sn=m;

        double s[sn],u[ldu*m],vt[ldvt*n],a[lda*n];

        int count=0;
        for(int i=0;i<(X_big.at(0)).size();++i)
        {
                for(int j=0;j<X_big.size();++j)
                {
                        a[count++]=X_big[j][i];
                }
        }

        lwork = -1;
        dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        work = (double*)malloc( lwork*sizeof(double) );
        /* Compute SVD */
        dgesvd_( "All", "All", &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork,&info );
        if( info > 0 ) {
                printf( "The algorithm computing SVD failed to converge.\n" );
        }
        vector<double> temp;

        for(int i=0;i<sn;++i)
        {
                sigma.push_back(s[i]);
        }

        for(int i=0;i<m;++i)
        {
                for(int j=0;j<sn;++j)
                {
                        temp.push_back(u[i+j*ldu]);
                }
                X_svd.push_back(temp);
                temp.clear();

        }

return 31;
}





int retrun_X_MATRIX(vector<vector<double> > &X_MATRIX,vector<vector<double> > X_svd,int count)
{
vector<double> temp_svd;
for(int i=0;i<X_svd.size();++i)
        {
                for(int j=0;j<count;++j)
                {
                        temp_svd.push_back(X_svd[i][j]);
                }
                X_MATRIX.push_back(temp_svd);
                temp_svd.clear();
        }
return 32;
}




int original_plus_MOR_sending_X_MATRIX(vector<vector<double> > X_MATRIX,double frequency_max,double frequency_min, double points)
{	
double interval;
	interval=(frequency_max-frequency_min)/(points-1);
	string in_1,in_2,in_3;
	cout<<"ENTER THE MAIN NETLIST FILE : ";
	cin>>in_1;

//	cout<<"enter the file containing  INPUT TERMINAL : " ;
//	cin>>in_2;
//	cout<<"enter the file containing OUTPUT TERMINALS : ";
//	cin>>in_3;
	cout<<"1) MAIN NETLIST FILE : "<<in_1;//<<endl<<"2) INPUT TERMINAL FILE : "<<in_2<<endl<<"3) OUTPUT TERMINAL FILE : "<<in_3<<endl;

	cout<<" DID U CHANGE THE FREQUENCY"<<endl/*<<"REMEMBER ARRAY IS A BITCH"<<endl*/<<"First enter RGC irrespective of the order,then enter L,then independent sources and then dependent sources"<<endl<<"R or r for resistance"<<endl<<"G or g for admittance"<<endl<<"C or c for capacitance"<<endl<<"L or l for inductance"<<endl<<"J or j for current source"<<endl<<"E or e for VCVS"<<endl<<"Z or z for VCCS"<<endl<<"H or h for CCVS"<<endl<<"V or v for independent voltage source"<<endl<<"PRESS ENTER TWICE AFTER DONE ENTERING THE CIRCUIT INFORAMTION IN THE FORM OF SPICE INPUT"<<endl;
	vector<complex<double> > u_MATRIX;
	vector<vector<double> >G_MATRIX;
	vector<vector<double> >C_MATRIX;
	vector<vector<double> >B_MATRIX;
	G_MATRIX.clear();
	C_MATRIX.clear();
	B_MATRIX.clear();


	//reading the input and output node from files created from pul.c code        
	vector<double> input_node,output_node;
	read_nodes("input_terminal.txt", input_node);
	read_nodes("output_terminal.txt", output_node);
	//printing the nodes
	for(int i=0;i<input_node.size();++i)
		cout<<input_node[i]<<" ";
	cout<<endl;
	for(int i=0;i<output_node.size();++i)
		cout<<output_node[i]<<" ";
	cout<<endl;

	populate_matrices(in_1, G_MATRIX, C_MATRIX, B_MATRIX,u_MATRIX);
	//	cout<<"G MATRIX"<<endl;	
	//	print(G_MATRIX);
	//	cout<<endl;
	//	cout<<"C MATRIX"<<endl;
	//	print(C_MATRIX);
	cout<<"THE NUMBER OF VARIABLE ARE :"<<G_MATRIX.size()<<endl;
	//	cout<<endl<<"B MATRIX"<<endl;
	//	print(B_MATRIX);
	//	cout<<endl<<"u complex vector"<<endl;
	//	printcomplexvector(u_MATRIX);
	/*converting ucomplex vectro matrix into ucomplex array matrix*/
	complex<double> uarraymatrix[(u_MATRIX).size()][1];
	for(int i=0;i<(u_MATRIX).size();++i)
		uarraymatrix[i][0]=u_MATRIX[i];
	//printing tht ucomplexaray
	//	cout<<endl<<"ucomplexarray"<<endl;
	//	for(int i=0;i!=(u_MATRIX).size();++i)
	//		cout<<uarraymatrix[i][0]<<endl;
	//converting the B VECTOR MATRIX to B ARRAY MATRIX
	complex<double> B_COMPLEXARRAYMATRIX[(B_MATRIX).size()][(u_MATRIX).size()];
	for(int i=0;i<(B_MATRIX).size();++i)
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

/*	//comment for original solving starts here	
	char answe;
		cout<<"DO YOU WANT TO SOLVE FOR THE ORIGINAL VARIABLES(y/n) :";
	cin>>answe;
	if(answe=='y')
	{
	                double t1=clock();

	fstream result11("SVD_MATLAB_result_at_output_node_at_multiple_frequencies.txt",fstream::out);
	if(result11.is_open() )
	{
		double frequency;
		int timer1=0;
		while(timer1<points)
		{
			cout<<"original"<<timer1+1<<endl;	
			frequency=frequency_min+(interval*timer1);
			result11<<frequency<<'\t';

			//combining G_MATRIX and C_MATRIX into GplussC_MATRIX
			complex<double> GplussC_MATRIX[(G_MATRIX).size()][(G_MATRIX).size()];
			for(int i=0;i<(G_MATRIX).size();++i)
			{
				for(int j=0;j<(G_MATRIX).size();++j)
				{
					GplussC_MATRIX[i][j]=complex<double>(G_MATRIX[i][j],(2*pi*frequency*C_MATRIX[i][j]));
				}
			}
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
			for(int i=0;i<n;++i)
			{
				for(int j=0;j<nrhs;++j)
				{
					original_solution[i][j]=b[i+j*ldb];
				}
			}
			for(int i=0;i<G_MATRIX.size();++i)
			{
				for(int k=0;k<output_node.size();++k)
				{
					if(i==(output_node[k]-1)){

						if((original_solution[i][0]).imag()>=0)
							result11<<(original_solution[i][0]).real()<<'+'<<(original_solution[i][0]).imag()<<'i'<<'\t';
						else
							result11<<(original_solution[i][0]).real()<<(original_solution[i][0]).imag()<<'i'<<'\t';

					}

				}
			}
			result11<<';'<<endl;

			++timer1;
		}//for loop of frequency ends here

	}
	else cout<<"unable to open file"<<endl;
double t2=clock();
                cout<<" TIME taken to compute the original model is : "<<(t2-t1)/double(CLOCKS_PER_SEC) << "seconds "<<endl;
	cout<<"done computing Ax=b of original solution for different frequencies and writing for MATLAB READING in SVD_MATLAB_result_at_output_node_at_multiple_frequencies.txt"<<endl;
}
*/	//comment for original model solving ends here

cout<<"checking the mtrices generated is orthonormal or not"<<endl;
                vector<vector<double> > X_MATRIX_TRANSPOSE,IDENTITY;
                vector_transpose(X_MATRIX,X_MATRIX_TRANSPOSE);
                vectormultiplication(X_MATRIX_TRANSPOSE,X_MATRIX,IDENTITY);
                cout<<"checking the orthonormality"<<endl;
                fstream wr_o("SVD_orthonormality_check.txt",fstream::out);

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
                cout<<" wrting the results of orthonormality check in SVD_orthonormality_check.txt"<<endl;


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
double t3=clock();
		fstream result22("SVD_MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt",fstream::out);
		if(result22.is_open() )
		{
			double frequency;
			int timer2=0;	
			while(timer2<points)
			{
				cout<<"computing MOR at "<<timer2+1<<" frequency point"<<endl;			
				frequency=frequency_min+(interval*timer2);
				result22<<frequency<<'\t';
				// Combining G_cap and C_cap into G_capplussC_cap
				complex<double> G_capplussC_cap[(G_cap).size()][(G_cap).size()];
				for(int i=0;i<(G_cap).size();++i){
					for(int j=0;j<(G_cap).size();++j)
					{
						G_capplussC_cap[i][j]=complex<double>(G_cap[i][j],(2*pi*frequency*C_cap[i][j]));
					}
				}

				//copying elements of GcapplusCcap in a one dimension array acap
				complex<double> a_cap[(G_cap.size())*(G_cap.size())];
				int k_cap=0;
				for(int i=0;i<(G_cap).size();++i)
					for(int j=0;j<(G_cap).size();++j)
						a_cap[k_cap++]=G_capplussC_cap[j][i];

				//converting the B_cap VECTOR MATRIX to B_cap ARRAY MATRIX
				complex<double> B_cap_COMPLEXARRAYMATRIX[(B_cap).size()][(u_MATRIX).size()];
				for(int i=0;i<(B_cap).size();++i)
					for(int j=0;j<(u_MATRIX).size();++j)
						B_cap_COMPLEXARRAYMATRIX[i][j]=complex<double>(B_cap[i][j],0);

				//multiplying the B_COMPLEXARRAY AND uarraymatrix
				complex<double> B_capucomplexmultiple[(B_cap).size()][1];
				for(int i=0;i<(B_cap).size();++i)
					for(int j=0;j<(u_MATRIX).size();++j)
						B_capucomplexmultiple[i][0]+=B_cap_COMPLEXARRAYMATRIX[i][j]*uarraymatrix[j][0];


				//copying elements of B_capucomplexmultiple into b_cap one dimension matrix
				complex<double> b_cap[(B_cap).size()];
				for(int i=0;i<(B_cap).size();++i)
					b_cap[i]=B_capucomplexmultiple[i][0];

				//computing a_capx=b_cap using cgesv routine

				int n_cap=G_cap.size();
				int nrhs_cap=1;

				int lda_cap=n_cap;
				int ldb_cap=n_cap;
				int info_cap;
				int ipiv_cap[n_cap];
				zgesv_( &n_cap, &nrhs_cap, a_cap, &lda_cap, ipiv_cap, b_cap, &ldb_cap, &info_cap );
				if( info_cap > 0 ) {
					printf( "The diagonal element of the triangular factor of A,\n" );
					printf( "U(%i,%i) is Jero, so that A is singular;\n", info_cap, info_cap );
					printf( "the solution could not be computed.\n" );
					return( 1 );
				}

				//converting X_matrix into X complex arary matrices

				complex<double> X_complexarray[X_MATRIX.size()][(X_MATRIX[0]).size()];
				for(int i=0;i<X_MATRIX.size();++i)
					for(int j=0;j<(X_MATRIX[0]).size();++j)
						X_complexarray[i][j]=complex<double>(X_MATRIX[i][j],0);

				//converting b_cap from inverse into x_cap 2 dimension array
				complex<double> x_cap[G_cap.size()][1];

				for(int i=0;i<n_cap;++i)
					for(int j=0;j<nrhs_cap;++j)
						x_cap[i][j]=b_cap[i+j*ldb_cap];

				complex<double>MORsolutioncompare[X_MATRIX.size()][1];
				for(int i=0;i<X_MATRIX.size();++i)
					for(int j=0;j<(X_MATRIX[0]).size();++j)
						MORsolutioncompare[i][0]+=X_MATRIX[i][j]*x_cap[j][0];
				for(int i=0;i<G_MATRIX.size();++i)
				{
					for(int k=0;k<output_node.size();++k)
					{
						if(i==(output_node[k]-1)){
							if((MORsolutioncompare[i][0]).imag()>=0)              
								result22<<(MORsolutioncompare[i][0]).real()<<'+'<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';
							else
								result22<<(MORsolutioncompare[i][0]).real()<<(MORsolutioncompare[i][0]).imag()<<'i'<<'\t';

						}
					}
				}
				result22<<';'<<endl;
				++timer2;			
			}
			result22.close();

		}
		else cout<<"unable top open file"<<endl;
double t4=clock();
        cout<<" TIME taken to compute the MOR from SVD MATRIX solution is :"<<(t4-t3)/double(CLOCKS_PER_SEC) << "seconds"<<endl; 
		cout<<"done computing Ax=b of MOR solution for different frequencies and writing for MATLAB READING in SVD_MATLAB_MORresult_at_output_node_at_multiple_frequencies.txt "<<endl;		
		//cout<<"the size of Original solution is :"<<original_solution.size()<<" and the size of MOR solution compare is :"<<MORsolutioncompare.size()<<endl;
		//	out<<" the original and MOR solution written in result.out"<<endl;

		return 33;
}





















