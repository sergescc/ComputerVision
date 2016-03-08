
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

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define NUM_THREADS 1
#define HEADER_SIZE 20
#define SINGLE_BYTE 1
#define FILTER_ORDER 2

void * ProcessRows (void * args);
void * ProcessColumns (void * args);
void * CenterRow (void * args);
void * LFilterRow (void * args);


typedef struct
{
        unsigned * counter;
        unsigned maxY;
        unsigned maxX;
        pthread_mutex_t * counterLock;
        float ** data;
        int cntl;
}ProcArgs;




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

unsigned char ** ReadImage(
 char * fileName, unsigned xSize, unsigned ySize, unsigned * numRows)
{

        unsigned char ** img;                   //pointer to hold the address of the image matrix
        unsigned bytesRead;                             // variable counting the bytes read per line
        unsigned i, j;
                                              // COUNTER VARIABLES
        FILE * openedFile;

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


void Fourier2D(
 float ** img, unsigned xSize, unsigned ySize, int cntl)
{

        pthread_mutex_t counterLock;
        unsigned counter, i;
        pthread_t threads[NUM_THREADS];

        ProcArgs threadArgs;
        pthread_mutex_init(&counterLock, NULL);
        counter = 0;

        threadArgs.counter = &counter;
        threadArgs.counterLock = &counterLock;
        threadArgs.maxY = ySize;
        threadArgs.maxX = xSize;
        threadArgs.data = img;
        threadArgs.cntl = cntl;




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
                    four1(&rows->data[row][0]-1 , rows->maxX, rows->cntl);
                }

        }
}

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
                column = *(columns->counter);
                *(columns->counter) += 2;
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
                        four1(temp -1 , columns->maxY, columns->cntl);
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
}


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
}

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
