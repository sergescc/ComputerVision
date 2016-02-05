#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>


int main (int argc, char **  argv)
{
	unsigned char ** rows;
	
	unsigned xSize;

	unsigned ySize; 	

	FILE * fileHandle;

	FILE * outputFile;

	int i, j;

	if (argc < 4)
	{
		printf("Invalid Parameters");
	}

	fileHandle = fopen(argv[1], r);

	outputFile = fopen(argv[2], w);

	if (fileHandle < 0)
	{
		printf("Error Could Not Open file");
		exit(-1);
	}
	
	fseek(fileHandle, 0L, SEEK_END);
	fileSize = ftell(fileHandle);
	fseek(fileHandel, 0L, SEEK_SET);

	rows = (unsigned char **) malloc( sizeof(unsigned char *) * ySize);

	for ( i = 0; i < ySize; i ++)
	{
		rows[i] = (unsigned char *) malloc (sizeof(unsigned char) * xSize);
		for ( j = 0; j < xSize; j++)
		{
			fread(&rows[i][j], sizeof(char) 
		}
	}


			

}
