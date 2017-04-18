#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gauss_seidel.c"

int gauss_seidel(int n, float er, float **A, float *b);

float er=0.01;

void combinationUtil(float **A, int arr[], int data[], int start, int end,
                     int index, int r, float* b)
{
    int i, j;
    // Current combination is ready to be printed, print it
    if (index == r)
    {
    	float **temp_A;
    	temp_A = (float **)malloc(r * sizeof(float));
	    for (i=0; i<r; i++)
	         temp_A[i] = (float *)malloc(r * sizeof(float));
        for (j=0; j<r; j++)
        {
        	for(i=0; i<r; i++)
        	{
        		temp_A[j][i]=A[j][data[i]];
        		printf("%.2f ", temp_A[j][i]);
        	}
        	printf("\n");
        }
        gauss_seidel(r, er, temp_A, b);
        printf("\n-------------------\n");
        return;
    }
 
    for (i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(A, arr, data, i+1, end, index+1, r, b);
    }
}

int main()
{
	int n, m, i, j, *arr;
	float **A, *b;

	printf("Enter No of Equations : ");
	scanf("%d", &m);

	printf("Enter No of Unknwons : ");
	scanf("%d", &n);

	A = (float **)malloc(m * sizeof(float));
    for (i=0; i<m; i++)
         A[i] = (float *)malloc(n * sizeof(float));

    arr=(int *)malloc(n*sizeof(int));

	printf("Enter the matrix A : \n");
	for (i=0; i<m; i++){
		for (j=0; j<n; j++)
			scanf("%f", &A[i][j]);
	}

	printf("Enter the column vector 'b' : \n");
	b= (float *)malloc(m*sizeof(float));
	for (i=0; i<m; i++)
		scanf("%f", &b[i]);


	for (i=0; i<n; i++)
		arr[i]=i;


	int data[m]; 
	combinationUtil(A , arr, data, 0,n-1, 0, m, b);
	
	return 0;
}