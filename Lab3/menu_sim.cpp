//				Lab3 Question

#include <iostream>
#include <iomanip> 
#include <vector>

using namespace std;

int min_index(float *lastrow, int length)
{
	int index=0;
	for (int i=1; i<length; i++)
	{
		if (lastrow[i]<lastrow[index])
			index=i;
	}
	if (lastrow[index]>=0)
		return -1;
	else
		return index;
}

int min_col_index(float *temp, int length)
{
	int index=0;
	for (int i=1; i<length; i++)
	{
		if (temp[i]<temp[index])
			index=i;
	}
	return index;
}

void print(int n, int m, float **A)
{
	int i, j;
	cout<<"----------------"<<endl;
	for (i=0; i<m; i++){
		for (j=0; j<(n); j++)
			cout<<setw(10)<<setprecision(5)<<A[i][j];
		cout<<endl;
	}
	// cout<<"----------------"<<endl;
}

void simplex(float **activeA, int n, int m)
{
	float lastrow[n+1], *temp;
	int i, j, pivot_col=0, pivot_row, count;
	bool stop=true, unbounded=false;

	for (i=0; i<n; i++)
	{
		if (activeA[m][i]<=0)
			stop=false;

	}

	if (stop==true){
		cout<<endl<<"Optimal Solution is "<<activeA[m][n]<<endl;
		
		for (i=0; i<m; i++)
		{
			if (activeA[i][n+1]<50)
				cout<<"x"<<activeA[i][n+1]<<" is "<<activeA[i][n]<<endl;
		}

		for (i=0; i<n; i++)
		{
			if (activeA[m+1][i]<50)
				cout<<"x"<<activeA[m+2][i]<<" is "<<0<<endl;
		}
		return;
	}

	pivot_col=min_index(activeA[m], n+1);

	//If the last element corresponding to pivot choosen pivot col is zero, Alternate solution exists.
	if (activeA[m][pivot_col]==0)
	{
		cout<<endl<<"Alternate Solution Exists. Optimal Solution is "<<activeA[m][n]<<endl;
		for (i=0; i<m; i++)
		{
			if (activeA[i][n+1]<50)
				cout<<"x"<<activeA[i][n+1]<<" is "<<activeA[i][n]<<endl;
		}

		for (i=0; i<n; i++)
		{
			if (activeA[m+1][i]<50)
				cout<<"x"<<activeA[m+2][i]<<" is "<<0<<endl;
		}
		
		return;
	}


	//Unboundedness Checking
	count=0;
	for (i=0; i<m; i++)
	{
		if (activeA[i][pivot_col]<=0)
			count++;
	}
	if (count==m){
		cout<<"Solution is Unbounded"<<endl;
		return;
	}

	//creating temp array for choosing pivot_row
	temp=new float [m];
	for (i=0; i<m; i++)
	{
		if (activeA[i][pivot_col]>0)
			temp[i]=activeA[i][n]/activeA[i][pivot_col];
		else
			temp[i]=1000000;
	}
	pivot_row=min_col_index(temp, m);
	cout<<"Row Column"<<pivot_row<<"  "<<pivot_col<<endl;

	float swap;
	swap=activeA[m+1][pivot_col];
	activeA[m+1][pivot_col]=activeA[pivot_row][n+1];
	activeA[pivot_row][n+1]=swap;

	//make new table
	float pivot=activeA[pivot_row][pivot_col];
	for (i=0; i<=m; i++)
	{
		for(j=0; j<=n; j++)
		{
			if (i!= pivot_row && j!= pivot_col)
			{
				activeA[i][j]=(activeA[i][j]*activeA[pivot_row][pivot_col]-activeA[pivot_row][j]*activeA[i][pivot_col])/activeA[pivot_row][pivot_col];
			}
		}
	}

	for (i=0; i<=n; i++)
		activeA[pivot_row][i]=activeA[pivot_row][i]/pivot;

	for (i=0; i<=m; i++)
		activeA[i][pivot_col]=-activeA[i][pivot_col]/pivot;
	activeA[pivot_row][pivot_col]=1/pivot;
	
	print(n+1, m+1, activeA);


	simplex(activeA, n, m);

}


int main()
{
	int n, m, i, j;
	float **A, *b, **activeA, *z;

	cout<<"Enter No of Equations : ";
	cin>>m;

	cout<<"Enter No of Unknwons : ";
	cin>>n;

    A= new float* [m]; 
    activeA = new float* [m+1+1];

    for (i=0; i<m; i++){
        A[i]= new float [n];
        activeA[i]= new float [n+1+1];
    }
    activeA[m]= new float [n+1];
    activeA[m+1]= new float [n+1];


	cout<<"Enter the matrix A : \n";
	for (i=0; i<m; i++){
		for (j=0; j<n; j++){
			cin>>A[i][j];
			activeA[i][j]=A[i][j];
		}
	}

	cout<<"Enter the column vector 'b' : \n";
    b=new float [m];
	for (i=0; i<m; i++){
		cin>>b[i];
		activeA[i][n]=b[i];
	}

	for (i=0; i<n; i++)
	{
		activeA[m+1][i]=i+1;
	}

	for (i=0; i<m; i++)
	{
		activeA[i][n+1]=51+i;
	}

	cout<<"Enter the objective function z (with the constant term) :\n";
	z=new float [n];
	for (i=0; i<n; i++){
		cin>>z[i];
		activeA[m][i]=-z[i];
	}
	cin>>activeA[m][n];


	print(n+2, m+2, activeA);


	simplex(activeA, n, m);

	return 0;


}
