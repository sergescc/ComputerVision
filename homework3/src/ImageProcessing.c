//////////////////////// ImageProcessing.c /////////////////////////////////////////
/*
By:   Sergio Coronado
        16.484 Computer Vision
        Assignemnt #2
        Part 1

PURPOSE:
  Module Holding functions for image Manipulations

*/

///////////////////////// Includes /////////////////////////////////////////////

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include "CursorCntl.h"
#include "ImageProcessing.h"

//////////////////////// Macros ////////////////////////////////////////////////


#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

///////////////////////// Constants ////////////////////////////////////////////

#define NUM_THREADS 10      //Number thrad used for simultaneous image processing
#define HEADER_SIZE 20      // Size of BMP header buffer
#define SINGLE_BYTE 1       // Singel Byte Lenght
#define FILTER_ORDER 2      // Order of Butterworth Filter
#define BACKGROUND_SEGMENTS 3
#define SOBEL_FILTER_SIZE 3
#define MAX_INTENSITY 256
#define BLUR_MAX  2
#define BLUR_STEP 1

///////////////////////// Function Protoypes  //////////////////////////////////

void * ProcessRows (void * args);
void * ProcessColumns (void * args);
void * CenterRow (void * args);
void * LFilterRow (void * args);

///////////////////////// ProcArgs ////////////////////////////////////////////
/*
PURPOSER:
Sturcture used to pass arguments to different threads for multithreaded image
processing
*/
typedef struct
{
        unsigned * counter;
        unsigned maxY;
        unsigned maxX;
        pthread_mutex_t * counterLock;
        float ** data;
        int cntl;
}ProcArgs;

///////////////////////// Four1 ////////////////////////////////////////////////
/*
PURPOSE:
Perform one dimensional fourire transform
Function used was aquired form 16.484 course materials

INPUT:
    float data[]:   Array of complex numbers to transform
    int nn:         Number of complex elements in the array
    int isign:      sign of the fourier transfor exponent
                        -1: Inversed FFT
                         1: FFT

OUTPUT:
    float data[]: Array trasnformed
*/


void four1(float data[],int nn, int isign)
{
        int n,mmax,m,j,istep,i;
        double wtemp,wr,wpr,wpi,wi,theta;
        float tempr,tempi;
        n=nn << 1;
        j=1;
        for (i=1;i<n;i+=2) {
                if (j > i) {
                        SWAP(data[j],data[i]);
                        SWAP(data[j+1],data[i+1]);
                }
                m=n >> 1;
                while (m >= 2 && j > m) {
                        j -= m;
                        m >>= 1;
                }
                j += m;
        }
        mmax=2;
        while (n > mmax) {
                istep=2*mmax;
                theta=6.28318530717959/(isign*mmax);
                wtemp=sin(0.5*theta);
                wpr = -2.0*wtemp*wtemp;
                wpi=sin(theta);
                wr=1.0;
                wi=0.0;
                for (m=1;m<mmax;m+=2) {
                        for (i=m;i<=n;i+=istep) {
                                j=i+mmax;
                                tempr=wr*data[j]-wi*data[j+1];
                                tempi=wr*data[j+1]+wi*data[j];
                                data[j]=data[i]-tempr;
                                data[j+1]=data[i+1]-tempi;
                                data[i] += tempr;
                                data[i+1] += tempi;
                        }
                        wr=(wtemp=wr)*wpr-wi*wpi+wr;
                        wi=wi*wpr+wtemp*wpi+wi;
                }
                mmax=istep;
        }
        return;

}

///////////////////////// ReadImage ////////////////////////////////////////////////
/*
PURPOSE:
Read in an image from specified file

INPUT:
    Char * filename:    Filename and path  of image
    unsigned xSize:     Horizontal size of image
    unsigned ySize:     Vertical Size of Image
    unsigned * numRows: Variable that spefies number of row read in

OUTPUT:
    unsigned char **: Pointer to matrix of read in image
*/


unsigned char ** ReadImage(
 char * fileName, unsigned xSize, unsigned ySize, unsigned * numRows)
{

        unsigned char ** img;   //pointer to hold the address of the image matrix
        unsigned bytesRead;     // variable counting the bytes read per line
        unsigned i, j;          // COUNTER VARIABLES

        FILE * openedFile;      // File Handle

        openedFile = fopen(fileName, "r");

        if (openedFile == NULL)
        {
            printError("Error Could Not Open file\n");
            return NULL;
        }

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

        fclose(openedFile);

        return img;
}

///////////////////////// Fourier2D ////////////////////////////////////////////////
/*
PURPOSE:
Performs 2d Fourier Transform

INPUT:
    float ** img:       Poiinter to the vectorized  image matrix to transform
    unsigned xSize:     Horizontal size of image
    unsigned ySize:     Vertical Size of Image
    int cntl:           Direction of Fourier transform
                            -1: Invered FFT
                             1: FFT

OUTPUT:
    float  ** img : Transformed image
*/

void Fourier2D(
 float ** img, unsigned xSize, unsigned ySize, int cntl)
{

        pthread_mutex_t counterLock;    // lock to handle synchronize rows
        unsigned counter, i;            // row counter and loop counter
        pthread_t threads[NUM_THREADS]; // Array of threads used to process image

        ProcArgs threadArgs;            // Struct to pass arguments
        pthread_mutex_init(&counterLock, NULL);
        counter = 0;

        threadArgs.counter = &counter;
        threadArgs.counterLock = &counterLock;
        threadArgs.maxY = ySize;
        threadArgs.maxX = xSize;
        threadArgs.data = img;
        threadArgs.cntl = cntl;



    //For Inverse FFt
	if (cntl < 0 )
	{
                for ( i = 0; i < NUM_THREADS; i ++)
                {
                        pthread_create(&threads[i], NULL, &ProcessColumns, &threadArgs);
                }
                for ( i = 0; i < NUM_THREADS; i++)
                {
                        pthread_join(threads[i], NULL);
                }

                counter = 0;

                for ( i = 0; i < NUM_THREADS; i ++)
                {
                        pthread_create(&threads[i], NULL, &ProcessRows, &threadArgs);
                }
                for ( i = 0; i < NUM_THREADS; i++)
                {
                        pthread_join(threads[i], NULL);
                }
	}
    //For Forward FFT
	else
	{
                for ( i = 0; i < NUM_THREADS; i ++)
                {
                        pthread_create(&threads[i], NULL, &ProcessRows, &threadArgs);
                }
                for ( i = 0; i < NUM_THREADS; i++)
                {
                        pthread_join(threads[i], NULL);
                }

                counter = 0;

                for ( i = 0; i < NUM_THREADS; i ++)
                {
                        pthread_create(&threads[i], NULL, &ProcessColumns, &threadArgs);
                }
                for ( i = 0; i < NUM_THREADS; i++)
                {
                        pthread_join(threads[i], NULL);
                }
	}

}

///////////////////////// ProcessRows ////////////////////////////////////////////////
/*
PURPOSE:
Thread function fro processing thread rows

INPUT:
    void * args: Pointer to the Process thread arguments structures

OUTPUT:
    1D FFT of processed Rows
*/


void * ProcessRows ( void * args)
{
        unsigned row;

        ProcArgs * rows = (ProcArgs *) args;

        while (1)
        {
                pthread_mutex_lock(rows->counterLock);
                row = *(rows->counter);
                *(rows->counter) += 1;
                pthread_mutex_unlock(rows->counterLock);
                if (row >= rows->maxY)
                {
                    pthread_exit(0);
                }
                else
                {
                    four1(rows->data[row] - 1, rows->maxX, rows->cntl);
                }

        }
}

///////////////////////// ProcessColumns  ////////////////////////////////////////////////
/*
PURPOSE:
Thread function for processing colums performs a 1d fft on a column

INPUT:
    void * args: Pointer to the Process thread arguments structures

OUTPUT:
    1D FFT of processed Columns
*/

void * ProcessColumns (void * args)
{
        unsigned column, i, n;
        float * temp;
        unsigned row;

        ProcArgs * columns = (ProcArgs *) args;

        n = columns->maxY << 1;

        while (1)
        {
                pthread_mutex_lock(columns->counterLock);
                column = *(columns->counter) << 1;
                *(columns->counter) += 1;
                pthread_mutex_unlock(columns->counterLock);
                if ( column >= (columns->maxX << 1))
                {
                        pthread_exit(0);
                }
                else
                {

                        temp = (float * ) malloc(sizeof(float) * n);
                        if (temp == NULL)
                        {
                            printError("Insuffificient Memory \n ");
                            exit(0);
                        }
                        for (i = 0 ; i < n; i += 2)
                        {
                                row = i >> 1;
                                temp[i] = columns->data[row][column];
                                temp[i +1] = columns->data[row][column + 1 ];
                        }
                        four1(temp - 1, columns->maxY, columns->cntl);
                        for (i = 0 ; i < n; i += 2)
                        {
                                row = i >> 1;
                                columns->data[row][column] = temp[i];
                                columns->data[row][column + 1] = temp[i + 1];
                        }
                        free(temp);
                }

        }
}

///////////////////////// Vectorize Image ////////////////////////////////////////////////
/*
PURPOSE:
Convert real image to real and imaginary vectored image

INPUT:
    unsigned char *** img:  Image to convert
    unsigned xSize:         Horizontal Size
    unsigned ySize:         Vertical Size

OUTPUT:
    float ** vectorized image
*/

float ** VectorizeImage (unsigned char ** img, unsigned  xSize, unsigned ySize)
{
        unsigned i, j, n;

        float ** newImg;
        newImg = (float **) malloc( sizeof(float *) * ySize);
        n = xSize << 1;
        for ( i = 0; i < ySize; i ++)
        {
            newImg[i] = (float *) malloc (sizeof(float) * n);
            for ( j = 0; j < n; j +=2)
            {
                newImg[i][j] = img[i][j >> 1];
                newImg[i][j+1] = 0;
            }
        }
        return newImg;
}

///////////////////////// Normalize Image ////////////////////////////////////////////////
/*
PURPOSE:
Calculate magnitude of complex image vectors and nromalize image

INPUT:
    float ** img:           Image to
    unsigned xSize:         Horizontal Size
    unsigned ySize:         Vertical Size

OUTPUT:
    unsigned char ** Normalized Image
*/

unsigned char ** NormalizeImage ( float ** img, unsigned xSize, unsigned ySize)
{
        unsigned i, j, n;
        float val;
        float min = FLT_MAX;
        float max = FLT_MIN;
        float ** tempImg;
        unsigned char ** newImg;
        tempImg = (float ** ) malloc ( sizeof(float *) * ySize);
        n = xSize << 1;
        for ( i = 0; i < ySize; i++)
        {
                tempImg[i] = (float * ) malloc (sizeof(float) * xSize);
                for ( j = 0; j < n ; j += 2)
                {
                        val = (tempImg[i][j >> 1] = sqrt(pow(img[i][j],2) + pow (img[i][j+1],2)));
                        if ( val > max)
                        {
                                max = val;
                        }
                        else if (val < min)
                        {
                                min = val;
                        }
                }
        }
        newImg = (unsigned char **) malloc (sizeof(unsigned char *) * ySize);
        for (i = 0; i < ySize; i ++)
        {
                newImg[i] = (unsigned char * ) malloc (sizeof(unsigned char) * xSize);
                for ( j =0 ; j < xSize; j++)
                {
                        val = newImg[i][j] = roundf((tempImg[i][j]-min) * (UCHAR_MAX / (max - min)));
                }
        }

        DestroyFloatImage( tempImg, xSize, ySize);

        return newImg;
}


///////////////////////// Destroy Float  Image ////////////////////////////////////////////////
/*
PURPOSE:
Destroy float image matrix

INPUT:
    float ** img:           Image to
    unsigned xSize:         Horizontal Size
    unsigned ySize:         Vertical Size

OUTPUT:
    void
*/

void DestroyFloatImage (float **  img, unsigned xSize, unsigned ySize)
{
        unsigned i;
        if ( img != NULL)
        {
                for (i = 0; i < ySize; i++)
                {       if (img[i] != NULL)
                        {
                                free(img[i]);
                        }
                }
                free(img);
        }
        img = NULL;
}

///////////////////////// Destroy Int  Image ////////////////////////////////////////////////
/*
PURPOSE:
Destroy float image matrix

INPUT:
    float ** img:           Image to
    unsigned xSize:         Horizontal Size
    unsigned ySize:         Vertical Size

OUTPUT:
    void
*/

void DestroyIntImage (int **  img, unsigned xSize, unsigned ySize)
{
        unsigned i;
        if ( img != NULL)
        {
                for (i = 0; i < ySize; i++)
                {       if (img[i] != NULL)
                        {
                                free(img[i]);
                        }
                }
                free(img);
        }
        img = NULL;
}

///////////////////////// Destroy Integer Image ////////////////////////////////////////////////
/*
PURPOSE:
Destroy unsigned char  image matrix

INPUT:
    float ** img:           Image to
    unsigned xSize:         Horizontal Size
    unsigned ySize:         Vertical Size

OUTPUT:
    void
*/


void DestroyImage (unsigned char **  img, unsigned xSize, unsigned ySize)
{
        unsigned i;
        if ( img != NULL)
        {
                for (i = 0; i < ySize; i++)
                {       if (img[i] != NULL)
                        {
                                free(img[i]);
                        }
                }
                free(img);
        }
        img = NULL;
}



///////////////////////// CenterRow ////////////////////////////////////////////////
/*
PURPOSE:
Thread routine that goed row by row aplying a spatial transform to an image to
center the spectrum

INPUT:
    void * args:    // Pointer to the Processing thread structure that holds
                    // processing information
OUTPUT:
    transformed rwo to center spectrum
*/

void * CenterRow(void * args)
{
        int row;
        unsigned column;
        ProcArgs * rows = (ProcArgs * ) args;

        while(1)
        {
                pthread_mutex_lock (rows->counterLock);
                row = *(rows->counter);
                *(rows->counter) += 1;
                pthread_mutex_unlock ( rows->counterLock);
                if (row >= rows->maxY)
                {
                        pthread_exit(0);
                }
                else if (row % 2 == 0)
                {
                        for ( column =  2 ; column < (rows->maxX << 1); column += 4 )
                        {
                                rows->data[row][column] *= -1;
                                rows->data[row][column + 1] *= -1;
                        }
                }
                else
                {
                        for ( column = 0 ; column < (rows->maxX << 1); column +=4)
                        {
                                rows->data[row][column] *= -1;
                                rows->data[row][column + 1] *= -1;
                        }
                }
        }

}
///////////////////////// CenterSpectrum ////////////////////////////////////////////////
/*
PURPOSE:
Multi-threaded routing to center spectrum

INPUT:
    float ** img:       // vectorised image to center spectrum
    unsigned xSize:     // horizontal size of image
    unsgined ySize      // vertical size of Image

OUTPUT:
    Image with centeres spectrum

*/


void CenterSpectrum (float ** img, unsigned xSize, unsigned ySize )
{
        pthread_t tid[NUM_THREADS];
        pthread_mutex_t counterLock;
        unsigned counter, i;

        ProcArgs threadArgs;
        pthread_mutex_init(&counterLock, NULL);
        counter = 0;

        threadArgs.counter = &counter;
        threadArgs.counterLock = &counterLock;
        threadArgs.maxY = ySize;
        threadArgs.maxX = xSize;
        threadArgs.data = img;
        threadArgs.cntl = 0;

        for (i = 0 ; i < NUM_THREADS; i ++)
        {
                pthread_create(&tid[i], NULL, &CenterRow, &threadArgs);
        }
        for (i = 0 ; i < NUM_THREADS; i ++)
        {
                pthread_join(tid[i], NULL);
        }

}

///////////////////////// CenterSpectrum ////////////////////////////////////////////////
/*
PURPOSE:
Multi-threaded routing to Apply Lowpass filter

INPUT:
    float ** img:       // vectorised image to center spectrum
    unsigned xSize:     // horizontal size of image
    unsgined ySize      // vertical size of Image
    int cutOff          // cutoff frequency of filter

OUTPUT:
    Filtered image

*/

void ApplyLowPass (float ** img, unsigned xSize, unsigned ySize, int cutOff)
{
        pthread_t tid[NUM_THREADS];
        pthread_mutex_t counterLock;
        unsigned counter, i;

        ProcArgs threadArgs;
        pthread_mutex_init(&counterLock, NULL);
        counter = 0;

        threadArgs.counter = &counter;
        threadArgs.counterLock = &counterLock;
        threadArgs.maxY = ySize;
        threadArgs.maxX = xSize;
        threadArgs.data = img;
        threadArgs.cntl = cutOff;

        for (i = 0 ; i < NUM_THREADS; i ++)
        {
                pthread_create(&tid[i], NULL, &LFilterRow, &threadArgs);
        }
        for (i = 0 ; i < NUM_THREADS; i ++)
        {
                pthread_join(tid[i], NULL);
        }
}

///////////////////////// CenterSpectrum ////////////////////////////////////////////////
/*
PURPOSE:
Thread routine to apply filter row by row

INPUT:
    void * args:        // pointer to Process thread arguments

OUTPUT:
    Filtered Rows of image

*/

void * LFilterRow (void * args)
{
        unsigned row, column;
        ProcArgs * rows = (ProcArgs *) args;
        float u, v,f;

        while (1)
        {
                pthread_mutex_lock(rows->counterLock);
                row = *(rows->counter);
                *(rows->counter) += 1;
                pthread_mutex_unlock(rows->counterLock);
                if (row >= rows->maxY)
                {
                        pthread_exit(0);
                }
                for ( column = 0 ; column < (rows->maxX << 1); column += 2)
                {
                        u = row - rows->maxY/2.0;
                        v = column - rows->maxX/2.0;
                        f = sqrt(pow(u,2) + pow (v,2));
                        rows->data[row][column] *= (1 / (1 + pow((f/rows->cntl), 2 * FILTER_ORDER )));
                        rows->data[row][column + 1] *= (1/ (1 + pow((f/rows->cntl), 2 * FILTER_ORDER)));

                }

        }
}

///////////////////////// Ourput Image ////////////////////////////////////////////////
/*
PURPOSE:
Output image to a file

INPUT:
    char * filename:        filename and path to output image to
    unsigned char ** img:   Nromalized image ot output
    unsigned xSize:         horizontal Size of image
    unsigned ySize:         vertical sizxe of image

OUTPUT:
    File with image

*/

void OutputImage ( char * fileName, unsigned char ** img, unsigned xSize, unsigned ySize)
{
        unsigned i, j;
        char headerBuffer[HEADER_SIZE];
        FILE * outputFile;

        outputFile = fopen(fileName, "w");

        if (outputFile == NULL)
        {
                printError("Could Not Open or create ouput file\n");
        }

        printOK("Building Image\n");

        sprintf(headerBuffer, "P5\n%d %d\n255\n", xSize, ySize);

        fwrite(headerBuffer, sizeof(unsigned char), strlen(headerBuffer), outputFile);

        for ( i = 0; i < ySize; i++)
        {
                for (j = 0 ; j < xSize; j++)
                {
                        fwrite(&img[i][j], sizeof(unsigned char), SINGLE_BYTE, outputFile);
                }
        }

        fclose(outputFile);
}

unsigned char ** NormalizeIntImage ( int ** img, unsigned xSize, unsigned ySize)
{
        unsigned i, j;
        int val;
        int min = INT_MAX;
        int max = INT_MIN;
        unsigned char ** newImg;
        for ( i = 0; i < ySize; i++)
        {
                for ( j = 0; j < xSize ; j += 2)
                {
                        val = img[i][j];
                        if ( val > max)
                        {
                                max = val;
                        }
                        else if (val < min)
                        {
                                min = val;
                        }
                }
        }

        newImg = (unsigned char **) malloc (sizeof(unsigned char *) * ySize);
        for (i = 0; i < ySize; i ++)
        {
            newImg[i] = (unsigned char * ) malloc (sizeof(unsigned char) * xSize);
            for ( j =0 ; j < xSize; j++)
            {

                newImg[i][j] = roundf((img[i][j] - min) * (UCHAR_MAX / ((float)(max-min))));

            }
        }

        return newImg;
}

unsigned char ** NormalizeFloatImage ( float ** img, unsigned xSize, unsigned ySize)
{
        unsigned i, j;
        float val;
        float min = INT_MAX;
        float max = INT_MIN;
        unsigned char ** newImg;
        for ( i = 0; i < ySize; i++)
        {
                for ( j = 0; j < xSize ; j += 2)
                {
                        val = img[i][j];
                        if ( val > max)
                        {
                                max = val;
                        }
                        else if (val < min)
                        {
                                min = val;
                        }
                }
        }

        newImg = (unsigned char **) malloc (sizeof(unsigned char *) * ySize);
        for (i = 0; i < ySize; i ++)
        {
            newImg[i] = (unsigned char * ) malloc (sizeof(unsigned char) * xSize);
            for ( j =0 ; j < xSize; j++)
            {

                newImg[i][j] = roundf((img[i][j] - min) * (UCHAR_MAX / (max-min)));

            }
        }

        return newImg;
}

int ** AddBuffer ( unsigned char ** img, unsigned xSize, unsigned ySize, unsigned bfrSize)
{
    unsigned i,j;

    int ** newImage;

    newImage =  (int ** ) malloc(sizeof(int *) * (ySize + (bfrSize << 1)));

    for ( i = 0; i < (ySize + (bfrSize << 1)); i++)
    {
        newImage[i] = (int *) malloc(sizeof(int) *(xSize + (bfrSize << 1)));

        if (i < bfrSize)
        {
            for ( j = 0; j < (xSize + (bfrSize << 1)); j ++)
            {
                if ( j < bfrSize)
                {
                    newImage[i][j] = img[0][0];
                }
                else if ( j >= (bfrSize + xSize) )
                {
                    newImage[i][j] = img[0][xSize-1];
                }
                else
                {
                    newImage[i][j] = img[0][j-bfrSize];
                }
            }
        }
        else if (i >= (bfrSize + ySize))
        {
            for ( j = 0; j < (xSize + (bfrSize << 1)); j ++)
            {
                if ( j < bfrSize)
                {
                    newImage[i][j] = img[ySize-1][0];
                }
                else if ( j >= (bfrSize + xSize) )
                {
                    newImage[i][j] = img[ySize-1][xSize-1];
                }
                else
                {
                    newImage[i][j] = img[ySize-1][j-bfrSize];
                }
            }
        }
        else
        {
            for ( j = 0; j < (xSize + (bfrSize << 1)); j ++)
            {
                if ( j < bfrSize)
                {
                    newImage[i][j] = img[i-bfrSize][0];
                }
                else if ( j >= (bfrSize + xSize) )
                {
                    newImage[i][j] = img[i-bfrSize][xSize-1];
                }
                else
                {
                    newImage[i][j] = img[i-bfrSize][j-bfrSize];
                }
            }
        }

    }

    return newImage;

}

unsigned char ** ApplySobel ( unsigned char ** img, unsigned xSize, unsigned ySize)
{

    float verticalSobel[SOBEL_FILTER_SIZE][SOBEL_FILTER_SIZE] =  {{1, 0, -1},
                                                                  {2, 0, -2},
                                                                  {1, 0, -1}};
    float horizontalSobel[SOBEL_FILTER_SIZE][SOBEL_FILTER_SIZE] =  {{ 1, 2, 1},
                                                                    {0, 0, 0},
                                                                    {-1,-2,-1}};
    unsigned i,j;
    float ** filtered_h;
    float ** filtered_v;
    int ** filtered;
    unsigned char ** normalized;
    float ** filterVert;
    float ** filterHoriz;
    filterVert = (float **) malloc (sizeof(float *)* SOBEL_FILTER_SIZE);
    filterHoriz = (float **) malloc (sizeof(float *)* SOBEL_FILTER_SIZE);

    for (i = 0; i < SOBEL_FILTER_SIZE; i++)
    {
        filterVert[i] = (float *) malloc (sizeof(float) * SOBEL_FILTER_SIZE);
        filterHoriz[i] = (float *) malloc (sizeof(float)* SOBEL_FILTER_SIZE);
        for (j = 0; j < SOBEL_FILTER_SIZE; j++)
        {
            filterVert[i][j] = verticalSobel[i][j];
            filterHoriz[i][j] = horizontalSobel[i][j];
        }
    }
    filtered_h = ApplyConvolution (img, xSize, ySize, filterHoriz, SOBEL_FILTER_SIZE);

    filtered_v = ApplyConvolution (img, xSize, ySize, filterVert, SOBEL_FILTER_SIZE);

    filtered = (int ** ) malloc (sizeof(int * )* ySize);
    for (i =0 ; i < ySize; i++)
    {
        filtered[i] = (int *) malloc (sizeof (int) * xSize);
        for (j = 0; j < xSize; j++)
        {
            filtered[i][j] = roundf(sqrt(pow(filtered_h[i][j],2) + pow(filtered_v[i][j],2)));
        }
    }
    normalized = NormalizeIntImage(filtered, xSize, ySize);
    DestroyIntImage(filtered, xSize, ySize);
    DestroyFloatImage(filterVert, SOBEL_FILTER_SIZE, SOBEL_FILTER_SIZE);
    DestroyFloatImage(filterHoriz, SOBEL_FILTER_SIZE, SOBEL_FILTER_SIZE);

    return normalized;

}

void MakeBinary ( unsigned char ** img, unsigned xSize, unsigned ySize, unsigned threshold)
{
    unsigned i,j;

    for ( i = 0; i < ySize; i++)
    {
        for (j = 0; j < xSize; j ++)
        {
            if (img[i][j] > threshold)
            {
                img[i][j] = 255;
            }
            else
            {
                img[i][j] = 0;
            }
        }
    }
}

void ApplyNoise ( unsigned char ** img, unsigned xSize, unsigned ySize, unsigned char intensity)
{
    unsigned i,j;
    float noiseMultiplier;
    int noisyValue;

    srandom(time(NULL));

    for (i = 0 ; i < ySize; i++ )
    {
        for (j =0; j < xSize; j++)
        {
            noiseMultiplier = (float)(random() - (RAND_MAX / 2))/(RAND_MAX/2);
            noisyValue = img[i][j] + roundf(noiseMultiplier * intensity);
            if (noisyValue < 0)
            {
                img [i][j] = 0;
            }
            else if (noisyValue < UCHAR_MAX)
            {
                img[i][j] = noisyValue;
            }
            else
            {
                img[i][j] = UCHAR_MAX;
            }
        }
    }


}

unsigned char ** MakeStandard ( unsigned xSize, unsigned ySize, unsigned tiers)
{
    int i, j;
    int rMax, binSize;
    int binScale, r;
    unsigned char ** img;
    int p1, p2;
    unsigned totalSegments;

    totalSegments = tiers + BACKGROUND_SEGMENTS;

    img  =  (unsigned char  **) malloc ( sizeof(unsigned char *) * ySize);

    rMax = (int) sqrt(pow(xSize >> 1,2) + pow(ySize >> 1,2));

    binScale = UCHAR_MAX / totalSegments;

    binSize = rMax / totalSegments;


    for ( i = 0 ; i < ySize; i++)
    {
        img[i] = (unsigned char * ) malloc (sizeof(unsigned char ) * xSize);
        p2 = pow(abs(i-(ySize >> 1)),2);
        for (j = 0; j < xSize; j ++)
        {
            p1 = pow(abs(j-(xSize >> 1)),2);
            r = sqrt(p1 + p2);

            img[i][j] = UCHAR_MAX - (binScale * (r/binSize));
            if (img[i][j] <= (BACKGROUND_SEGMENTS * binScale))
            {
                img[i][j] = UCHAR_MAX - (binScale * (totalSegments - BACKGROUND_SEGMENTS));
            }
        }
    }

    return img;

}

float ** ApplyConvolution ( unsigned char ** img, unsigned xSize, unsigned ySize, float ** filter, unsigned fSize)
{

    int i,j,k,m, maxOffset;
    int ** buffered;
    float ** convoluted;
    float value;

    maxOffset = (fSize - 1)/2;

    buffered = AddBuffer(img, xSize, ySize, maxOffset);

    convoluted = (float **) malloc ( sizeof(float *)*ySize);

    for (i = maxOffset; i < ySize + maxOffset; i++)
    {
        convoluted[i - maxOffset] = (float *) malloc(sizeof(float) * xSize);


        for ( j = maxOffset ; j < xSize + maxOffset; j++)
        {
            value = 0;
            for (k = (maxOffset * -1) ; k <= maxOffset; k++)
            {

                for ( m = (maxOffset * -1); m <=  maxOffset; m++)
                {
                    value += buffered[i + k][j + m] * filter[k + maxOffset][m + maxOffset];
                }
            }
            convoluted[i - maxOffset][j - maxOffset] = value;

        }

    }

    DestroyIntImage(buffered,xSize + 2* maxOffset, ySize + 2 * maxOffset);

    return convoluted;
}

unsigned char ** ApplyGaussian (unsigned char ** img, unsigned xSize, unsigned ySize, unsigned scale)
{
    unsigned i,j;
    float ** gaussFilter;
    float ** filtered;
    unsigned char ** normalized;
    unsigned width;

    width = scale << 1;

    if (width % 2 == 0)
    {
        width ++;
    }


    gaussFilter = (float **) malloc (sizeof(float *) * width);


    for ( i = 0; i < width; i ++)
    {
        gaussFilter[i] =  (float *) malloc (sizeof(float) * width);
        for ( j = 0; j < width; j ++)
        {
            gaussFilter[i][j] =  .5 * M_1_PI * (1/pow(scale,2)) * pow(M_E, (-1 * ((pow(j,2) + pow(i,2))/ (2 * pow(scale,2)))));
        }
    }

    filtered = ApplyConvolution(img, xSize, ySize, gaussFilter, width);

    normalized = NormalizeFloatImage(filtered, xSize, ySize);

    DestroyFloatImage(filtered, xSize, ySize);
    DestroyFloatImage(gaussFilter, width, width);

    return normalized;
}



typedef struct ErrorCalculatorArgs
{
    unsigned  * threshold;
    unsigned  * errorArray;
    pthread_mutex_t * threshold_lock;
    unsigned char ** img;
    unsigned char ** standard;
    unsigned xSize;
    unsigned ySize;
}ErrArgs;

void * CalculateErrors ( void * args )
{
    ErrArgs * err_a =  args;
    unsigned threshold;
    unsigned char ** working;
    unsigned * errorArray = err_a->errorArray;

    pthread_mutex_lock(err_a->threshold_lock);
    while (1)
    {
        if (*(err_a->threshold) >= UCHAR_MAX)
        {
            break;
        }
        threshold = *(err_a->threshold);
        
        *(err_a->threshold) += 1;
       
        pthread_mutex_unlock(err_a->threshold_lock);

        working = CopyImage(err_a->img, err_a->xSize, err_a->ySize);

        MakeBinary(working, err_a->xSize, err_a->ySize, threshold);

        errorArray[threshold]= CompareImage(err_a->standard, working, err_a->xSize, err_a->ySize);

        DestroyImage( working, err_a->xSize, err_a->ySize);

        pthread_mutex_lock(err_a->threshold_lock);

    }

    pthread_mutex_unlock(err_a->threshold_lock);

    pthread_exit(0);

}

unsigned FindOptimalThreshold(unsigned char ** standard, unsigned char ** img, unsigned xSize, unsigned ySize, int * errorValue)
{
    unsigned i;
    unsigned min = UINT_MAX;
    pthread_t threads[NUM_THREADS];
    unsigned errorArray[MAX_INTENSITY];
    pthread_mutex_t threshold_lock;
    unsigned threshold =0;
    ErrArgs args;
    pthread_mutex_init(&threshold_lock, NULL);
    args.threshold = &threshold;
    args.errorArray = errorArray;
    args.threshold_lock = &threshold_lock;
    args.img = img;
    args.standard = standard;
    args.xSize = xSize;
    args.ySize = ySize;

    for (i = 0; i < NUM_THREADS; i++)
    {
        pthread_create(&threads[i],NULL, &CalculateErrors, &args);
    }

    for ( i = 0; i < NUM_THREADS; i ++)
    {
        pthread_join(threads[i], NULL);
    }

    for (i = 0; i < MAX_INTENSITY; i ++)
    {
        if (errorArray[i] < min)
        {
            min = errorArray[i];
            threshold = i;
        }
    }

    *(errorValue) = min;

    return threshold;
}


unsigned CompareImage ( unsigned char ** standard, unsigned char ** img, unsigned xSize, unsigned ySize)
{
    unsigned errors;
    unsigned i,j;

    errors = 0;
    for (i = 0; i < ySize; i++)
    {
        for (j= 0 ; j < xSize; j ++)
        {
            if (img[i][j] != standard[i][j])
            {
                errors++;
            }
        }
    }

    return errors;
}

unsigned char ** CopyImage ( unsigned char ** img, unsigned xSize, unsigned ySize)
{
    unsigned char ** copied;
    unsigned i, j;

    copied = (unsigned char **) malloc (sizeof(unsigned char *)* ySize);
    for (i = 0; i < ySize; i++)
    {
        copied[i] = (unsigned char * ) malloc (sizeof( unsigned char ) * xSize);
        for (j=0; j < xSize; j++)
        {
            copied[i][j] = img[i][j];
        }
    }

    return copied;
}
typedef struct BlurErrors
{
    unsigned  ** errorMatrix;
    unsigned char ** standard;
    unsigned char ** img;
    unsigned xSize;
    unsigned ySize;
    pthread_mutex_t * blur_lock;
    int * blur;

}BErrArgs;

void * CalculateErrorsWithBlur( void * argsIn)
{
    unsigned i;
    BErrArgs * err_a = argsIn;
    unsigned ** errorMatrix = err_a->errorMatrix;
    pthread_t threads[NUM_THREADS];
    pthread_mutex_t threshold_lock;
    ErrArgs args;
    pthread_mutex_init(&threshold_lock, NULL);
    unsigned char ** blurred;
    unsigned char ** blurredSobel;
    int blur;
    unsigned threshold;
    threshold = 0;
    args.threshold = &threshold;
    args.threshold_lock = &threshold_lock;
    args.standard = err_a->standard;
    args.xSize = err_a->xSize;
    args.ySize = err_a->ySize;

    pthread_mutex_lock(err_a->blur_lock);
    while (1)
    {
        threshold = 0;
        if (*(err_a->blur) >= BLUR_MAX)
        {
            break;
        }
        blur = *(err_a->blur);
        *(err_a->blur) += BLUR_STEP;

        args.errorArray = errorMatrix[(blur/BLUR_STEP)- 1];

        pthread_mutex_unlock(err_a->blur_lock);

        blurred = ApplyGaussian(err_a->img , err_a->xSize, err_a->ySize, blur);
        blurredSobel = ApplySobel(blurred, err_a->xSize, err_a->ySize);
        args.img = blurredSobel;

        for (i = 0; i < NUM_THREADS; i++)
        {
            pthread_create(&threads[i],NULL, &CalculateErrors, &args);
        }

        for ( i = 0; i < NUM_THREADS; i ++)
        {
            pthread_join(threads[i], NULL);
        }

        DestroyImage(blurred, err_a->xSize, err_a->ySize);
        DestroyImage(blurredSobel, err_a->xSize, err_a->ySize);

        pthread_mutex_lock(err_a->blur_lock);
    }
    pthread_mutex_unlock(err_a->blur_lock);
    pthread_exit(0);

}


unsigned FindOptimalWithBlurThreshold(unsigned char ** standard, unsigned char ** img, unsigned xSize, unsigned ySize, int * errorValue, int *blurScale )
{


    unsigned i, j;
    unsigned min = UINT_MAX;
    pthread_t threads[NUM_THREADS];
    unsigned ** errorMatrix;
    pthread_mutex_t blur_lock;
    unsigned minThreshold;
    int blur = BLUR_STEP;
    BErrArgs args;
    pthread_mutex_init(&blur_lock, NULL);

    errorMatrix = (unsigned ** ) malloc (sizeof(unsigned *) * (BLUR_MAX/BLUR_STEP));
    for (i = BLUR_STEP ; i < BLUR_MAX; i += BLUR_STEP)
    {
        errorMatrix[(i/BLUR_STEP) - 1] = (unsigned *) malloc (sizeof(unsigned) * MAX_INTENSITY);
    }
    args.blur = &blur;
    args.errorMatrix = errorMatrix;
    args.blur_lock = &blur_lock;
    args.standard = standard;
    args.xSize = xSize;
    args.ySize = ySize;
    args.img = img;

    for (i = 0; i < NUM_THREADS; i++)
    {
            pthread_create(&threads[i],NULL, &CalculateErrorsWithBlur, &args);
    }

    for ( i = 0; i < NUM_THREADS; i ++)
    {
        pthread_join(threads[i], NULL);
    }

    for (i = BLUR_STEP ; i < BLUR_MAX; i += BLUR_STEP)
    {
        for (j = 0 ; j < MAX_INTENSITY; j ++)
        {
            if (errorMatrix[(i/BLUR_STEP) -1 ][j] < min)
            {
                min = errorMatrix[(i/BLUR_STEP) -1 ][j];
                minThreshold = j;
                blur = i;
            }
        }
    }

    *(errorValue) = min;
    *(blurScale) = blur;

    return minThreshold;
}