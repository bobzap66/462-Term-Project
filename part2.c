#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 256


//returns a random double between -1 and 1
double r2()
{
	double rv;
	rv = rand() / (double)RAND_MAX;
	rv = rv - rand()%2;
	return rv;
}
//prints the matrix
void printMatrix(double m[], int n)
{
	int i, j;
	printf("\nMatrix Given");
	for(i = 0; i < n*n; i++)
	{
		j = i%n;
		if(j == 0)
			printf("\n");
		printf("%f ", m[i]);
	}
	printf("\n");
}
//multiplies two square matrices when given a size
void matrixMatrix(double* m1, double* m2, double* result, int n)
{
	int i, j, k;//index variables
	double sum;
	sum = 0;
	for (i = 0; i < n; i++) 
	{
      		for (j = 0; j < n; j++) 
		{
       			for (k = 0; k < n; k++) 
			{
				result[i*n+j] += m1[i*n+k]*m2[k*n+j];
        		}
		}
	}
}
//initializes the matrices for this example
void initializeMatrices (double* m1, double* m2, double* result)
{
	int i, j; // index variables
	for (i=0; i<N; i++)
	{
		for (j=0; j<N; j++) 
		{
			m1[i*N+j] = r2();
			m2[i*N+j] = r2();
 			result[i*N+j] = 0;
		}
	}
} 



int main(int argc, char* argv[])
{
	int gridSize;
	double t1, t2;
	int i, j, k;
	double elapsed;
	srand(1);
	double *m1;
	double *m2;
	double *result;
	int dimensions[2];
	int periodic[2];
	int subDims[2];
	int rank, numprocs, sourceProcess;
	double* m1Block;
	double* m2Block;
	double* resultBlock;
	double* m1ABlock;
	double* resultRow;
	int blockSize;
	double* mRow;
	int pivot;
	int nextProc;
	int prevProc;
	MPI_Comm GridComm;
	MPI_Status stat;
	MPI_Comm colComm;
	MPI_Comm rowComm;
	int gridCoords[2];
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//hard coded in, could not get the sqrt() function from math.h to work
	if(numprocs == 4)
		gridSize = 2;
	else if(numprocs == 16)
		gridSize = 4;
	else if(numprocs == 64)
		gridSize = 8;
	else
		gridSize = 16;

	//initialize everything
	dimensions[0] = dimensions[1] = gridSize;
	periodic[0] = periodic[1] = 0;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periodic, 1, &GridComm);
	MPI_Cart_coords(GridComm, rank, 2, gridCoords);
	subDims[0] = 0;
	subDims[1] = 1;
	MPI_Cart_sub(GridComm, subDims, &rowComm);
	subDims[0] = 1;
	subDims[1] = 0;
	MPI_Cart_sub(GridComm, subDims, &colComm);
	blockSize = N/gridSize;
	m1Block = (double*) malloc(sizeof(double)*blockSize*blockSize);
	m2Block = (double*) malloc(sizeof(double)*blockSize*blockSize);
	resultBlock = (double*) malloc(sizeof(double)*blockSize*blockSize);
	m1ABlock = (double*) malloc(sizeof(double)*blockSize*blockSize);
	for(i = 0; i < blockSize*blockSize; i++)
		resultBlock[i] = 0;
	if(rank == 0) //initialize matrices
	{
		m1 = (double*) malloc(sizeof(double)*N*N);
		m2 = (double*) malloc(sizeof(double)*N*N);
		result = (double*) malloc(sizeof(double)*N*N);
		initializeMatrices(m1, m2, result);
		t1 = MPI_Wtime();
		if(numprocs == 1)
		{
			matrixMatrix(m1, m2, result, N);
			t2 = MPI_Wtime();
			elapsed = (t2 - t1) * 1000;
			printMatrix(result, N);
			printf("Time = %f milliseconds on one processor\n", elapsed);
		}
	}
	//distribute data
	mRow = (double*) malloc(sizeof(double)*blockSize*N);
	if(gridCoords[1] == 0)
	{
		MPI_Scatter(m1, blockSize*N, MPI_DOUBLE, mRow, blockSize*N, MPI_DOUBLE, 0, colComm);
	}
	for(i = 0; i < blockSize; i++)
		MPI_Scatter(&mRow[i*N], blockSize, MPI_DOUBLE, &(m1ABlock[i*blockSize]), blockSize, MPI_DOUBLE, 0, rowComm);
	if(gridCoords[1] == 0)
	{
		MPI_Scatter(m2, blockSize*N, MPI_DOUBLE, mRow, blockSize*N, MPI_DOUBLE, 0, colComm);
	}
	for(i = 0; i < blockSize; i++)
		MPI_Scatter(&mRow[i*N], blockSize, MPI_DOUBLE, &(m2Block[i*blockSize]), blockSize, MPI_DOUBLE, 0, rowComm);
	//multiply matrices and switch data sets among processors
	for(i = 0; i < gridSize; i++)
	{
		pivot = (gridCoords[0] + i)%gridSize;
		if(gridCoords[1] == pivot)
		{
			for(j = 0; j < blockSize*blockSize; j++)
			{
				m1Block[j] = m1ABlock[j];
			}
		}
		MPI_Bcast(m1Block, blockSize*blockSize, MPI_DOUBLE, pivot, rowComm);
		
		matrixMatrix(m1Block, m2Block, resultBlock, blockSize);
		nextProc = gridCoords[0] + 1;
		if(gridCoords[0] == gridSize -1)
			nextProc = 0;
		prevProc = gridCoords[0]-1;
		if(gridCoords[0] == 0)
			prevProc = gridSize -1;

		MPI_Sendrecv_replace(m2Block, blockSize*blockSize, MPI_DOUBLE, nextProc, 0, prevProc, 0, colComm, &stat);
	}
	//gather data back to root
	resultRow = (double*) malloc(sizeof(double)*N*blockSize);
	for(i = 0; i < blockSize; i++)
		MPI_Gather(&resultBlock[i*blockSize], blockSize, MPI_DOUBLE, &resultRow[i*N], blockSize, MPI_DOUBLE, 0, rowComm);
	if(gridCoords[1] == 0)
		MPI_Gather(resultRow, blockSize*N, MPI_DOUBLE, result, blockSize*N, MPI_DOUBLE, 0, colComm);
	if(rank == 0)//print final matrix and time
	{
		t2 = MPI_Wtime();
		elapsed = (t2 - t1) * 1000;
		printMatrix(result, N);
		printf("Time = %f milliseconds on %d processors\n", elapsed, numprocs);
	}
	//cleanup
	if(rank == 0)
	{
		free(m1);
		free(m2);
		free(result);
	}
	free(mRow);
	free(m1ABlock);
	free(m1Block);
	free(m2Block);
	free(resultRow);
	free(resultBlock);
	MPI_Finalize();
	return 0;
}

