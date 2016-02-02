/* 
Description		:Make array with random numbers, and find
				highest number
Author			:Andreas Slyngstad
								*/

#include <stdio.h>

void main()
	{
		matrix2();
	}
/*
int *matrix()
{
	int n = 5;
	int i;
	double*A_storage=(double*)malloc(n*n*sizeof(double));
	double **A = (double**) malloc(n *sizeof(double*));
	for (i=0; i<n; i++)
		A[i] = &(A_storage[i*n]);
}
*/
int matrix2()
{
	int n = 3;
	int i, j;
	int a[n][n];
	a[2][2] = 1;
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < n; j++)
		{
			printf(" array[%d][%d] %d \n",i,j, &a[i][j]);
		}

	}
}