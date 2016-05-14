#include"print.h"


/*Function for printing the 2 dimension vector*/
int print( vector<vector<double> > &p)
{
	for(int i=0; i<(p.size());++i)
	{
		for(int j=0; j<((p[i]).size());++j)
			cout<<p[i][j]<<" ";
		cout<<";"<<endl;
	}
	return 3;
}

/* Function for printing one dimension complex vector*/
int printcomplexvector(vector<complex<double> >&p)
{
	for(int i=0; i<(p.size());++i)
	{
		// for(int j=0; j!=((p[i]).size());++j)
		cout<<p[i]<<" ";
		cout<<";"<<endl;
	}
	return 4;
}

int print_eigenvalues( char* desc, int n,double* wr,double* wi ) {
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
	return 20;
}

/* Auxiliary routine: printing a matrix */
int print_matrix( char* desc, int m, int n, complex<double> *a, int lda ) 
{
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%f,%f)", a[i+j*lda].real(), a[i+j*lda].imag() );
		printf( "\n" );
	}
	return 1;
}
/* Auxiliary routine: printing a vector of integers */
int print_int_vector( char* desc, int n, int* a ) 
{
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
	return 2;
}


int print_complex_vector(vector<vector<complex<double> > > a)
{
	for(int i=0;i<a.size();++i)
	{
		for(int j=0;j<(a.at(0)).size();++j)
		{
			cout<<a[i][j]<<" ";
		}
		cout<<endl;
	}
	return 22;
}


int file_write(string filename, vector<vector<double> > V)
{
        fstream write(filename.c_str(),fstream::out);
        if (write.is_open())
        {
                for(int i=0;i<V.size();++i)
                {
                        for(int j=0;j<(V.at(0)).size();++j)
                        {
                                write<<V[i][j]<<" ";
                        }
                        write<<endl;
                }
        }
        else cout<<"unable to open file"<<endl;

        return 28;
}

