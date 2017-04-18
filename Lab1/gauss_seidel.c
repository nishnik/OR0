#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int gauss_seidel(int n, float er, float **A, float *b)
{
	int i, j, key, iteration=0;
	float *x0, sum, *x;

	//Initialize solution
	x0=(float *)malloc(n*sizeof(float));
	x=(float *)malloc(n*sizeof(float));

	for (i=0; i<n; i++)
		x0[i]=0;

	while(1)
	{	
		iteration++;
		key=0;
		for (i=0; i<n; i++)
		{
			sum=b[i];
			for (j=0; j<n; j++)
			{
				if (j!=i)
					sum=sum-A[i][j]*x0[j];
			}
			
			x[i]=sum/A[i][i];
			if (((x[i]-x0[i])/x[i])>er || ((x[i]-x0[i])/x[i])<(-er))
				key=1;
			x0[i]=x[i];
		}
		if (key==0)
		{	
			printf("Solution is : \n");
			for (i=0; i<n; i++)
			{
				printf("%f ", x0[i]);
			}
			printf("\nIterations : %d", iteration);
			break;
		}
	}

	return 0;
}

