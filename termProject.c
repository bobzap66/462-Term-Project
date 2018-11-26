#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>

double r2()
{
    double rv;
    rv = rand() / (double)RAND_MAX ;
    rv = rv - rand()%2;
    return rv;
}

void print_matrix(int n, double m[n][n])
{
    int i, j; // i = row; j = column;
    printf("\nMatrix Given\n");
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
            printf("%f", m[i][j]);
        printf("\n");
    }
}

void print_vector(int n, double v[n])
{
    int i;
    printf("\nVector Given\n");
    for (i=0; i<n; i++)
        printf("%f", v[i]);
    printf("\n");
}


void matrix_x_vector(int n, double y[n], double x[n][n], double A[n][n])
{
    int i, j; // i = row; j = column;
  //  printf("\nResulted Matrix of [M]*[V]\n");
    // Load up A[n][n]
    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            A[j][0] += x[j][i] * y[i];
       //     printf("%f", A[j][0]);
        }
      //  printf("\n");
    }
   // printf("\n");
    // Print A[n][n]
   // for (i=0; i<n; i++)
   // printf("%f\n", A[i][0]);
}



int main(int argc, char* argv[])
{
	double t1, t2;
    	int n = 256;
    	int i, j;
    	double m1[256][256];
    	double m2[256][256];
    	double result[256][256];
    	int rank, numprocs, count, remainder, myRowSize;
    	int* sendcounts = NULL;
    	int* senddispls = NULL;
    	int* recvcounts = NULL;
    	int* recvdispls = NULL;	
    	float elapsed;
    	srand(1);
    	MPI_Init (&argc, &argv);
    	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if(rank == 0)
	{	
    		for(i = 0; i < 256; i++)
    		{
        		for(j = 0; j < 256; j++)
        		{
            			m1[i][j] = r2();
            			m2[i][j] = r2();
				result[i][j] = 0;
        		}
    		}/*
		if(numprocs == 1)
    		{
			t1 = MPI_Wtime();
    			matrix_x_vector(n, vector, matrix, result);
			t2 = MPI_Wtime();
    			elapsed = (t2 - t1)*1000;
    			printf("Time = %f milliseconds on one processor\n", elapsed);
		}*/
		sendcounts = (int*)malloc(sizeof(int)*numprocs); 
		senddispls = (int*)malloc(sizeof(int)*numprocs); 
		recvcounts = (int*)malloc(sizeof(int)*numprocs); 
		recvdispls = (int*)malloc(sizeof(int)*numprocs); 
		count = 256/numprocs;
		remainder = 256 - count * numprocs;
		int prefixSum = 0;
		for(i = 0; i < numprocs; i++)
		{
			recvcounts[i] = (i < remainder) ? count+1 : count;
		        sendcounts[i] = recvcounts[i] * 256;
			recvdispls[i] = prefixSum;
			senddispls[i] = prefixSum * 256;
			prefixSum += recvcounts[i];
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	MPI_Bcast(vector, 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//	printf("My rank is %d and the first element in my vector is %f\n", rank, vector[0]);
	if (rank != 0)
	{
		count = 256 / numprocs;
		remainder = 256 - count * numprocs;
	}
	myRowSize = rank < remainder ? count + 1: count;
	double* matrixPart = (double*)malloc(sizeof(double) * myRowSize * 256);
	MPI_Scatterv(matrix, sendcounts, senddispls, MPI_DOUBLE, matrixPart, myRowSize * 256, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank == 0)
	{
		free(sendcounts);
		free(senddispls);
	}
	double* resultPart = (double*)malloc(sizeof(double) * myRowSize);
//#pragma omp parallel for
	for (i = 0; i < myRowSize; i++)
	{
		resultPart[i] = 0;
		for (j = 0; j < 256; j++)
		{
			resultPart[i] += matrixPart[i * 256 + j] * vector[j];
		}
	}
	free(matrixPart);
	MPI_Gatherv(resultPart, myRowSize, MPI_DOUBLE, result, recvcounts, recvdispls, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	free(resultPart);
	if(rank == 0)
	{
		free(recvcounts);
		free(recvdispls);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	t2 = MPI_Wtime();
	MPI_Finalize();
    	elapsed = (t2 - t1)*1000;
	if( rank == 0)
	{
		printf("Results:\n");
		for( i = 0; i < 256; i++)
			printf("%f ", result[0][i]);
		printf("\n");
    		printf("Time = %f milliseconds on %d processors\n", elapsed, numprocs);
	}
	return 0;
}
   


