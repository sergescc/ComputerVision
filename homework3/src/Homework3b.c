//////////////////////// homework3b.c /////////////////////////////////////////
/*
By:   Sergio Coronado
        16.484 Computer Vision
        Assignemnt #2
        Part 1

PURPOSE:
  Program reads in an image performs a FFT outputs the spectrum and then
  performs the inverse FFT

USAGE:

  HW2a <input-file> <ospectrum-file> <reverse-file> <xSize> <ySize>

*/

///////////////////////// Includes /////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include "ImageProcessing.h"
#include "CursorCntl.h"
#include <string.h>

///////////////////////// Constants ////////////////////////////////////////////


#define NUM_ARGS 3
#define STANDARD_TIERS 4
#define NOISE_TIERS_SIZE 10 
#define NOISE_TIER_MAX 128 
#define FILE_NAME_BUFFER_SIZE 100
#define X_SIZE 512
#define Y_SIZE 512
#define DATA_BUFFER_SIZE 100


///////////////////////// Main /////////////////////////////////////////////////


int main ( int argc, char ** argv)
{
  unsigned i;
	unsigned char ** standard;
	unsigned char ** sobelStandard;
  unsigned char ** working; 
  unsigned char ** workingSobel;
  unsigned char ** workingBlurred;
  unsigned workingThreshhold;
  int workingBlur;
  int error;
					           // Matrix Holding Original Image Values

  FILE * csvfile;

  char fileNameBuffer[FILE_NAME_BUFFER_SIZE];
  char dataBuffer[DATA_BUFFER_SIZE];

	if (argc <  NUM_ARGS)
	{
  	printError("Usage: HW3b <ResultFileNameRoot> <ImagesNoBlurLocation> <ImagesBlurLocation>\n");
  	exit(0);
	}

	printOK("Creating  Image sobelStandard \n");

	standard = MakeStandard(X_SIZE, Y_SIZE, STANDARD_TIERS);

  printOK("Applying Sobel To Standard\n");

  sobelStandard = ApplySobel(standard, X_SIZE, Y_SIZE);
  MakeBinary(sobelStandard, X_SIZE, Y_SIZE, 30);

  // sprintf(fileNameBuffer, "%s-NoBlur.csv", argv[1]);
  // csvfile = fopen(fileNameBuffer, "w");

  // if (csvfile == NULL)
  // {
  //   printError("Unable to Create CSV File\n");
  //   exit(-1);
  // }

  // sprintf(dataBuffer, "Noise Level, Error, Threshold\n");
  // fwrite(dataBuffer, sizeof(char), strlen(dataBuffer), csvfile);

  // for ( i = NOISE_TIERS_SIZE; i < NOISE_TIER_MAX; i += NOISE_TIERS_SIZE)
  // {
  //   sprintf(dataBuffer, "Calculating Threshold without blur for Noise Level: %d \n", i);
  //   printOK(dataBuffer);
  //   working = CopyImage (standard, X_SIZE, Y_SIZE);
  //   ApplyNoise(working, X_SIZE, Y_SIZE, i);
  //   workingSobel = ApplySobel(working, X_SIZE, Y_SIZE);
  //   workingThreshhold = FindOptimalThreshold(sobelStandard, workingSobel, X_SIZE, Y_SIZE, &error);
  //   MakeBinary(workingSobel, X_SIZE, Y_SIZE, workingThreshhold);
  //   sprintf(fileNameBuffer, "%sNoBlur-Noise-%d-Threshold-%d-Error-%d", argv[2], i, workingThreshhold, error );
  //   OutputImage(fileNameBuffer, workingSobel, X_SIZE, Y_SIZE );
  //   DestroyImage(working, X_SIZE, Y_SIZE);
  //   DestroyImage(workingSobel, X_SIZE, Y_SIZE);
  //   sprintf(dataBuffer, "%d, %d, %d,\n", i, error, workingThreshhold);
  //   fwrite(dataBuffer, sizeof(char), strlen(dataBuffer), csvfile);

  // }
  // fclose(csvfile);


  sprintf(fileNameBuffer, "%s-WithBlur.csv", argv[1]);
  csvfile = fopen(fileNameBuffer, "w");

  if (csvfile == NULL)
  {
    printError("Unable to Create CSV File\n");
    exit(-1);
  }

  sprintf(dataBuffer, "Noise Level, Error, Threshold, Blur\n");
  fwrite(dataBuffer, sizeof(char), strlen(dataBuffer), csvfile);

  for ( i = NOISE_TIERS_SIZE; i < NOISE_TIER_MAX; i += NOISE_TIERS_SIZE)
  {
    sprintf(dataBuffer, "Calculating Threshold with blur for Noise Level: %d \n", i);
    printOK(dataBuffer);

    working = CopyImage (standard, X_SIZE, Y_SIZE);
    ApplyNoise(working, X_SIZE, Y_SIZE, i);

    workingThreshhold = FindOptimalWithBlurThreshold(sobelStandard, working, X_SIZE, Y_SIZE, &error, &workingBlur);

    sprintf(fileNameBuffer, "%sWithBlur-Noise-%d-Threshold-%d-Error-%d-Blur-%d", argv[3], i, workingThreshhold, error, workingBlur);

    workingBlurred = ApplyGaussian(working, X_SIZE, Y_SIZE, workingBlur);
    workingSobel = ApplySobel(workingBlurred, X_SIZE, Y_SIZE);
    MakeBinary(workingSobel, X_SIZE, Y_SIZE, workingThreshhold);
    OutputImage(fileNameBuffer, workingSobel, X_SIZE, Y_SIZE );


    DestroyImage(working, X_SIZE, Y_SIZE);
    DestroyImage(workingSobel, X_SIZE, Y_SIZE);
    DestroyImage(workingBlurred, X_SIZE, Y_SIZE);
    sprintf(dataBuffer, "%d, %d, %d,%d\n", i, error, workingThreshhold, workingBlur);
    fwrite(dataBuffer, sizeof(char), strlen(dataBuffer), csvfile);

  }
  fclose(csvfile);

  sprintf(fileNameBuffer, "%s-Standard", argv[2] );
  OutputImage(fileNameBuffer, standard, X_SIZE, Y_SIZE );

  sprintf(fileNameBuffer, "%s-Standard-Sobel", argv[2] );
  OutputImage(fileNameBuffer, sobelStandard, X_SIZE, Y_SIZE );

	printOK("Clean up\n");

	DestroyImage(standard, X_SIZE, Y_SIZE);
  DestroyImage(sobelStandard, X_SIZE, Y_SIZE);

	exit(0);

}
