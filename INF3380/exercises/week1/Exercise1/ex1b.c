/* 
Description		:Make array with random numbers, and find
				highest number
Author			:Andreas Slyngstad
								*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int *array()
{
	int n = 10;
	static int numbers[10]; /*Static to change variable later */
	int i;
	time_t t;

	/*Initialize random number generator */
	srand((unsigned) time(&t));
	for(i = 0; i < n; i++ )
	{
		 numbers[i] = rand() % 20;
	}
	return numbers;
}

void minmax()
{
	int *ar;
	ar = array();
	int min, max, k;
	min = ar[0];
	max = ar[0];
	
	printf("Random generated array \n");
	for(k = 0; k < 10; k++)
	{

		if(min > ar[k])
		{
			min = ar[k];
		}

		if(max < ar[k])
		{
			max = ar[k];
		}
	printf("%d \n", ar[k]);
	}
	printf("The minimum number is %d \n", min);
	printf("The maximum number is %d \n", max);
}

int main()
{
	int length = 10;
	array();
	minmax();
	
}

