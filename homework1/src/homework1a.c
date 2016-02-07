#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CursorCntl.h"
#include <sys/wait.h>


#define SINGLE_BYTE 1
#define HEADER_SIZE 20
#define NUM_ARGS 3
#define CMD_LENGTH 2


static char *CMD = "/usr/bin/xv";

int main (int argc, char **  argv)
{
	unsigned char ** rows;
	unsigned xSize;
	unsigned ySize;
	unsigned numRows;
	unsigned char * zeroBuff;
	pid_t xvFork;
	char headerBuffer[HEADER_SIZE];
	int bytesRead;
	char * xvArgs[NUM_ARGS];
	int status;


	FILE * fileHandle;
	FILE * outputFile;

	int i, j;

	if (argc < 5)
	{
		printError("Invalid Parameters\n");
		exit (-1);
	}

	fileHandle = fopen(argv[1], "r");

	outputFile = fopen(argv[2], "w");

	xSize = atoi(argv[3]);

	ySize = atoi(argv[4]);

	if (fileHandle == NULL)
	{
		printError("Error Could Not Open file\n");
		exit (-1);
	}

	rows = (unsigned char **) malloc( sizeof(unsigned char *) * ySize);

	zeroBuff = (unsigned char *) calloc (xSize/2, sizeof(unsigned char ));

	numRows = 0;

	printOK("Reading Image\n");

	for ( i = 0; i < ySize; i ++)
	{
		rows[i] = (unsigned char *) malloc (sizeof(unsigned char) * xSize);
		bytesRead = fread(rows[i], sizeof(unsigned char), xSize , fileHandle );
		if (bytesRead < xSize)
		{
			if (feof(fileHandle))
			{
				printWarning("Reached Unexpected EOF Attemptin padding to size\n");
					for (j = (bytesRead -1); j < xSize; j++)
					{
						rows[i][j] = 0;
					}
			}
			else
			{
				printError("Error Reading File");
			}

		}
		numRows++;

	}

	printOK("Building Image\n");

	sprintf(headerBuffer, "P5\n%d %d\n255\n", xSize, ySize);

	fwrite(headerBuffer, sizeof(unsigned char), strlen(headerBuffer), outputFile);

	for ( i = 0; i < numRows; i++)
	{
		fwrite( zeroBuff, sizeof(unsigned char), xSize/2, outputFile);
		for (j = xSize/2 ; j < xSize; j ++)
		{
			fwrite(&rows[i][j], sizeof(unsigned char), SINGLE_BYTE, outputFile);
		}
	}

	fclose(outputFile);
	fclose(fileHandle);

	printOK("Image Built Starting XV\n");

	xvArgs[0] = CMD;
	xvArgs[1] = argv[2];
	xvArgs[2] = NULL;

	xvFork = fork();

	if (xvFork == 0)
	{
		status = execv(CMD, xvArgs);
		if (status < 0)
		{
			printError("Could Not Initiate XV\n");
		}

	}
	else
	{
		waitpid(xvFork, 0, 0);
	}

return 0;

}
