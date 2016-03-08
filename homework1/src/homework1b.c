///////////////////////// homework1b.c /////////////////////////////////////////
/*
By:		Sergio Coronado
				16.484 Computer Vision
				Assignemnt #1
				Part 2

PURPOSE:
	Program reads in an image specfied by user and outputs a file specified by user
	that has reduced to a quarter of the original size by keeping even pixels and
	discarding odd ones.

USAGE:

	HW1b <input-file> <output-file> <xSize> <ySize>

*/

///////////////////////// Includes /////////////////////////////////////////////


#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CursorCntl.h"
#include <sys/wait.h>

///////////////////////// Constant Definitions /////////////////////////////////

#define SINGLE_BYTE 1
#define HEADER_SIZE 20
#define NUM_ARGS 3
#define CMD_LENGTH 2

/* Command String for launching XV
	 Note:
	 		For Runnign on School Lab this must be : /usr/X11R6/bin/xv
*/

//static char *CMD = "/usr/bin/xv";
static char *CMD = "/usr/X11R6/bin/xv";

///////////////////////// Function Prototypes //////////////////////////////////


unsigned char ** ReadImage(
 FILE * openedFile, unsigned xSize, unsigned ySize, unsigned * numRows);

///////////////////////// Main /////////////////////////////////////////////////


int main (int argc, char **  argv)
{

	/////////////////////// Variables ////////////////////////////////////////////


	unsigned char ** img;						//Matrix Holding Image Values
	unsigned xSize;									//NUmber of horizontal pixels
	unsigned ySize;									//NUmber of vertical pixels
	unsigned numRows;								//rows Read in
	pid_t xvFork;										//pid for xv launched
	char headerBuffer[HEADER_SIZE]; //array contianign header params for the image
	char * xvArgs[NUM_ARGS];				//Arguments for XV
	int status;											//Status variable


	FILE * openedFile;							// file to be read in
	FILE * outputFile;							// file to be wrtten out

	int i, j;				//counter varaibles

	/* Check Correct number of argumetns */

	if (argc < 5)
	{
		printError("Invalid Parameters\n");
		exit (-1);
	}

	/* open image file */

	openedFile = fopen(argv[1], "r");

	xSize = atoi(argv[3]);

	ySize = atoi(argv[4]);

	if (openedFile == NULL)
	{
		printError("Error Could Not Open image file\n");
		exit (-1);
	}

	numRows = 0;

	printOK("Reading Image\n");

	/* Read Image */

	img = ReadImage(  openedFile, xSize, ySize, &numRows);

	/* Close input file and open output file */
	/* Designed in this order in case user wishes to over-write image */

	fclose(openedFile);

	outputFile = fopen(argv[2], "w");

	if (outputFile == NULL)
	{
		printError("Could Not Open or create ouput file\n");
		exit (-1);
	}

	printOK("Building Image\n");

	sprintf(headerBuffer, "P5\n%d %d\n255\n", xSize/2, ySize/2);

	fwrite(headerBuffer, sizeof(unsigned char), strlen(headerBuffer), outputFile);

	for ( i = 0; i < numRows; i += 2)
	{
		for (j = 0 ; j < xSize; j += 2)
		{
			fwrite(&img[i][j], sizeof(unsigned char), SINGLE_BYTE, outputFile);
		}
	}

	fclose(outputFile);

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

///////////////////////// ReadImage() //////////////////////////////////////////]
/*
PURPOSE
The read image function reads an image from a file of size specified in the
parameters
INPUT
	FILE * openedFile : pointer to the file being opened must be opened priot to
											call
	unsigned xSize	`	:	number of columns in image to be read
	unsigned ySize		: number of rows in image to be read
	unsigned char *** img 	: Pointer to be set to the

OUTPUT

	unsigned char ** img 	: read-in matrix of pixel densities

*/

unsigned char ** ReadImage(
 FILE * openedFile, unsigned xSize, unsigned ySize, unsigned * numRows)
{
	unsigned char ** img;			//pointer to hold the address of the image matrix
	unsigned bytesRead;				// variable counting the bytes read per line
	unsigned i, j;						// COUNTER VARIABLES

	*(numRows) = 0;

	img = (unsigned char **) malloc( sizeof(unsigned char *) * ySize);

	for ( i = 0; i < ySize; i ++)
	{
		img[i] = (unsigned char *) malloc (sizeof(unsigned char) * xSize);
		bytesRead = fread(img[i], sizeof(unsigned char), xSize , openedFile );
		if (bytesRead < xSize)
		{
			if (feof(openedFile))
			{
				printWarning("Reached Unexpected EOF Attemptin padding to size\n");
					for (j = (bytesRead -1); j < xSize; j++)
					{
						img[i][j] = 0;
					}
			}
			else
			{
				printError("Error Reading File");
			}

		}
		*(numRows) += 1;

	}
	return img;
}
