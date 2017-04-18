//				Lab6 : Dual Simplex Question

#include <iostream>
#include <iomanip> 
#include <vector>
#include <cmath>

using namespace std;

int iteration, usr_iter;

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

void simplex(float **activeA, int n, int m, int max_min, int choice)
{
	float lastrow[n+1], *temp;
	int i, j, pivot_col=0, pivot_row, count;
	bool stop=true, unbounded=false;

	iteration++;

	for (i=0; i<m; i++)
	{
		if (activeA[i][n]<0)
			stop=false;

	}

	if (stop==true){
		if (choice==4)
		{
			if (max_min==1)
				cout<<endl<<"Optimal Solution is "<<activeA[m][n]<<endl;
			if (max_min==2)
				cout<<endl<<"Optimal Solution is "<<(-1)*activeA[m][n]<<endl;
			
			for (i=0; i<m; i++)
			{
				if (activeA[i][n+1]<50)
					cout<<"x"<<activeA[i][n+1]<<" is "<<activeA[i][n]<<endl;
			}

			for (i=0; i<n; i++)
			{
				if (activeA[m+1][i]<50)
					cout<<"x"<<activeA[m+1][i]<<" is "<<0<<endl;
			}
		}
		return;
	}

	float temp1[m];
	for (i=0; i<m; i++)
		temp1[i]=activeA[i][n];

	pivot_row=min_index(temp1, m);

	//Unboundedness Checking
	count=0;
	for (i=0; i<n; i++)
	{
		if (activeA[pivot_row][i]>=0)
			count++;
	}
	if (count==n){
		if (choice==4)
			cout<<"Solution is Unbounded"<<endl;
		return;
	}

	//creating temp array for choosing pivot_row
	temp=new float [m];
	for (i=0; i<n; i++)
	{
		if (activeA[pivot_row][i]<0)
			temp[i]=abs(activeA[m][i]/activeA[pivot_row][i]);
		else
			temp[i]=1000000;
	}
	pivot_col=min_col_index(temp, n);
	// cout<<"Row Column"<<pivot_row<<"  "<<pivot_col<<endl;

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
	
	// print(n+1, m+1, activeA);

	if (choice==3 && iteration==usr_iter)
		print(n+1, m+1, activeA);	

	if (choice==2 && iteration==usr_iter)
	{
		cout<<"Basic Variables are : "<<endl;
		for (i=0; i<m; i++)
		{
			if (activeA[i][n+1]<50)
				cout<<"x"<<activeA[i][n+1]<<endl;
			else
				cout<<"z"<<activeA[i][n+1]-50<<endl;
		}

		cout<<"Non Basic Variables are :"<<endl;
		for (i=0; i<n; i++)
		{
			if (activeA[m+1][i]<50)
				cout<<"x"<<activeA[m+1][i]<<endl;
			else
				cout<<"z"<<activeA[m+1][i]-50<<endl;
		}
	}	


	simplex(activeA, n, m, max_min, choice);

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
    activeA[m]= new float [n+1+1];
    activeA[m+1]= new float [n+1+1];


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

	int max_min;
	cout<<"Enter the choice : 1 -- Maximize, 2 -- Minimize"<<endl;
	cin>>max_min;

	cout<<"Enter the objective function z :\n";
	z=new float [n];
	for (i=0; i<n; i++){
		cin>>z[i];
		if (max_min==1)
			activeA[m][i]=-z[i];
		if (max_min==2)
			activeA[m][i]=z[i];
	}
	z=new float [n];
	cin>>activeA[m][n];


	print(n+2, m+2, activeA);

	float tempA[m+2][n+2];
	for (i=0; i<m+2; i++)
		for (j=0; j<n+2; j++)
			tempA[i][j]=activeA[i][j];

	int choice;
	do{
	cout<<"Choose the Option:"<<endl;
	cout<<"1) Initial Table"<<endl;
	cout<<"2) Basic and Non-Basic at ith iteration"<<endl;
	cout<<"3) Table at ith iteration"<<endl;
	cout<<"4) Optimal Solution"<<endl;
	cout<<"5) Exit"<<endl;
	cin>>choice;
		if (choice==1){
			cout<<"Initial Table is :"<<endl;
			print(n+1, m+1, activeA);
		}
		if (choice==3 || choice==2)
		{
			cout<<"Enter the iteration: ";
			cin>>usr_iter;
		}
		iteration=0;
		simplex(activeA, n, m, max_min, choice);
		for (i=0; i<m+2; i++)
			for (j=0; j<n+2; j++)
				activeA[i][j]=tempA[i][j];

	}
	while(choice<5);


	

	return 0;


}
