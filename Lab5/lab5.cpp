//				Lab3 Question

#include <iostream>
#include <iomanip> 
#include <vector>

using namespace std;

int checkSolution;

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

void simplex(float **activeA, int n, int m, int *index, int max_min)
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
		checkSolution=1; // Solution exists
		return;
	}

	pivot_col=min_index(activeA[m], n);

	//If the last element corresponding to pivot choosen pivot col is zero, Alternate solution exists.
	if (activeA[m][pivot_col]==0)
	{
		checkSolution=2; // Alternate Solution Exists
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
		checkSolution=3;   // Unboundedness
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
	// cout<<"Row Column"<<pivot_row<<"  "<<pivot_col<<endl;

	index[pivot_row]=pivot_col;

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


	simplex(activeA, n, m, index, max_min);

}


int main()
{
	int n, m, i, j, *sol_index, geq, leq, eq, M=1000, max_min, choice;
	float **A, *b, **activeA, *z, temp;

	cout<<"Enter No of Equations : ";
	cin>>m;

	cout<<"Enter No of Unknwons : ";
	cin>>n;

    // arr = new int [n];
    cout<<"Enter the number of >=, <=, = contraints : "<<endl;
    cin>>geq>>leq>>eq;

    if (geq==0 && eq==0)
    	M=0;

    A= new float* [m]; 
    activeA = new float* [m+1];

    for (i=0; i<m; i++){
        A[i]= new float [n];
        activeA[i]= new float [n+geq+1];
    }
    activeA[m]= new float [n+geq+1];

	cout<<"Enter the matrix A in order (e.g first 1) >= 2) <= 3) = ,contraints) : \n";
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
		activeA[m][n+i]=M;
	}

	cout<<"Enter the column vector 'b' : \n";
    b=new float [m];
    temp=0;
	for (i=0; i<m; i++){
		cin>>b[i];
		temp=temp+b[i];
		activeA[i][n+geq]=b[i];
	}
	activeA[m][n+geq]= temp*M*(-1);

	cout<<"Enter the choice : 1 -- Maximize, 2 -- Minimize"<<endl;
	cin>>max_min;

	cout<<"Enter the objective function z :\n";
	z=new float [n];
	for (i=0; i<n; i++){
		cin>>z[i];
		// activeA[m][i]=-z[i];
		temp=0;
		for (j=0; j<m; j++)
			temp=temp+A[j][i];
		if (max_min==1)
			activeA[m][i]=-z[i]-temp*M;
		if (max_min==2)
			activeA[m][i]=z[i]+temp*M;
	}



	// print(n+geq+1, m+1, activeA);

	// for ()


	// print(n+1, m+1, activeA);
	sol_index = new int[m];

	for (i=0; i<m; i++)
		sol_index[i]=-1;

	simplex(activeA, n+geq, m, sol_index, max_min);

	cout<<"Enter your choice :"<<endl;
	cout<<"(a) List of all BFS"<<endl;
	cout<<"(b) Number of Iterations to solve the problem"<<endl;
	cout<<"(c) List of all Non-basic variables along with net evaluations in ith iteration"<<endl;
	cout<<"(d) List of Basic variables along with min ratios in ith iteration"<<endl;
	cout<<"(e) simplex table of ith iteration"<<endl;
	cout<<"(f) optimal solution (if exists otherwise generate report for infeasibility, unboundedness, alternative optimum etc.)"<<endl;
	cin>>choice;

	if (choice==5)
	{
		if (checkSolution==1)
			cout<<"Optimal Solution is "<<endl;
		else if (checkSolution==2)
			cout<<"Alternate Solution Exists. One optimal Solution is :"<<endl;
		if (checkSolution==3){
			cout<<"Unboundedness"<<endl;
			return 0;
		}
		
		if (max_min==1)
			cout<<activeA[m][n]<<endl;
		if (max_min==2)
			cout<<(-1)*activeA[m][n]<<endl;
		
		float *sol_;
		sol_ = new float [m];
		for (i=0; i<m; i++)
		{
			if (sol_index[i]>=0)
				sol_[sol_index[i]]=activeA[i][n];
		}
		for (i=0; i<m; i++)
		{
			cout<<"x"<<i+1<<" is "<<sol_[i]<<endl;
		}
	}

	return 0;


}
