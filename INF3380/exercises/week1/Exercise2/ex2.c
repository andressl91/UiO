/* 
Description		:Reading a file, and extracting data
Author			:Andreas Slyngstad
								*/

#include <stdio.h>
#include <stdlib.h>

int *getdata()
{
	FILE *fp;
	char data[10];
	static int temp[6];
	int as_int, i;

	fp = fopen("temp.txt", "r");

	for (i = 0; i < 6; i++)
	{
		fscanf(fp, "%s", data);
		fgets(data, 50, (FILE*)fp);
		as_int = atoi(data);
		temp[i] = as_int;
		
	}

	fclose(fp);
	return temp;
}	

int main()
{
	int *da;
	int i;
	da = getdata();
		
	for (i = 0; i < 6; i++)
	{
	printf("%d \n", da[i]);
	}
	return 0;
}