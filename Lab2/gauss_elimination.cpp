				//GAUSS_ELIMINATION.cpp

#include <iostream>
#include <iomanip> 
using namespace std;

void opr_divide(int n ,int row, float divide, float**A, float *b)
{
	int i;
	for (i=0; i<n; i++)
		A[row][i]=A[row][i]/divide;
	b[row]=b[row]/divide;
}

void opr_multi(int n ,int row, float multiply, float**A, float *b)
{
	int i;
	for (i=0; i<n; i++)
		A[row][i]=A[row][i]*multiply;
	b[row]=b[row]*multiply;
}

//R2 -> R2 + 2R1
void opr_row(int n, int row2, int row1, float multiply, float **A, float *b)
{
	int i;
	for (i=0; i<n; i++){
			A[row2][i]=A[row2][i]+A[row1][i]*multiply;
	}
	b[row2]=b[row2]+b[row1]*multiply;


}

void print(int n, float **A, float* b)
{
	int i, j;
	cout<<"----------------"<<endl;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++)
			cout<<setw(3)<<setprecision(3)<<A[i][j];
		cout<<" | "<<b[i];
		cout<<endl;
	}
	// cout<<"----------------"<<endl;
}

int findSwap(float **A, int row, int n, float *b)
{
	if (row==n-1){
		// cout<<"Last Row: Cannot find a row to swap"<<endl;
		return 0;
	}
	int i, swapRow, temp;
	for (i=row+1; i<n; i++)
		if (A[i][row]!=0){
			swapRow=i;
			break;
		}

	temp=b[row];
	b[row]=b[swapRow];
	b[swapRow]=temp;

	for (i=0; i<n; i++)
	{
		temp=A[row][i];
		A[row][i]=A[swapRow][i];
		A[swapRow][i]=temp;
	}

	return 1;

}

int gauss_elimination( int n,float **A, float *b, float *sol)
{
	int i, j, ret=-1;
	// float sol[n];
	print(n, A, b);

	// opr_divide(n, 0, .3, A);
	// print(n, A);

	// opr_multi(n, 1, .1, A);
	// print(n, A);

	// opr_row(n, 1, 0, 3, A);
	// print(n, A);

	// Gauss Elimination

	for (i=0; i<n; i++)
	{
		if (A[i][i]==0)
			ret=findSwap(A, i, n, b);

		if (ret==0)
		{
			print(n, A, b);
			cout<<"Unique Solution Does not Exist"<<endl;
			return 0;
		}

		opr_divide(n, i, A[i][i], A, b);
		for (j=i+1; j<n; j++)
		{
			opr_row(n, j, i, -A[j][i], A, b);
		}
		// print(n, A, b);

	}

	sol[n-1]=b[n-1];
	for (i=n-2; i>=0; i--)
	{
		sol[i]=b[i];
		for (j=i+1; j<n; j++)
			sol[i]=sol[i]-sol[j]*A[i][j];
	}

	// cout<<"Solution is : ";
	// for (i=0; i<n; i++)
	// 	cout<<setprecision(2)<<"  "<<sol[i];

	// cout<<endl;

	return 1;

	

}
