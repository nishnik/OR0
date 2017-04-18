#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

int main()
  {
  mat A = randu<mat>(4,5);
  mat B = randu<mat>(4,5);

  cout<< A.t() <<endl;
  cout<< A <<endl;
  
  cout << A*B.t() << endl;
  
  return 0;
  }

int main()
{
	int n, m, i, j, geq, leq, eq, M=1000, max_min;
	float **A, *b, **activeA, *z, temp;

	cout<<"Enter No of Equations : ";
	cin>>m;

	cout<<"Enter No of Unknwons : ";
	cin>>n;

    A= new float* [m]; 
    activeA = new float* [m+1+1+1];

    for (i=0; i<m; i++){
        A[i]= new float [n];
        activeA[i]= new float [n+geq+1+1];
    }
    activeA[m]= new float [n+geq+1];
    activeA[m+1]= new float [n+geq+1];
    activeA[m+2]= new float [n+geq+1];

	cout<<"Enter the matrix A in order (e.g first 1) >= 2) = 3) <= ,contraints) : \n";
	for (i=0; i<m; i++){
		for (j=0; j<n; j++){
			cin>>A[i][j];
			activeA[i][j]=A[i][j];
		}
	}

	for (i=0; i<geq; i++)
	{
		for (j=0; j<m; j++){
			if (i==j)
				activeA[j][n+i]=-1;
			else
				activeA[j][n+i]=0;
		}
	}

	cout<<"Enter the column vector 'b' : \n";
    b=new float [m];
	for (i=0; i<m; i++){
		cin>>b[i];
		activeA[i][n+geq]=b[i];
	}

	for (i=0; i<geq+eq; i++)
	{
		for (j=0; j<n+geq+1; j++)
			activeA[m][j]-=activeA[i][j];
	}

	for (i=0; i<n+geq; i++)
	{
		activeA[m+2][i]=i+1;
	}

	for (i=0; i<m; i++)
	{
		activeA[i][n+geq+1]=51+i;
	}

	cout<<"Enter the choice : 1 -- Maximize, 2 -- Minimize"<<endl;
	cin>>max_min;

	cout<<"Enter the objective function z :\n";
	z=new float [n];
	for (i=0; i<n; i++){
		cin>>z[i];
		if (max_min==1)
			activeA[m+1][i]=-z[i];
		if (max_min==2)
			activeA[m+1][i]=z[i];
	}