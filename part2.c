#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 256

int matrixToArray(int i, int j, int n)
{
	int k;
	k = i * n + j;
	return k;
}

int getI(int k, int n)
{
	int i;
	i = k / n;
	return i;
}

int getJ(int k, int n)
{
	int j;
	j = k%n;
	return j;
}
double r2()
{
	double rv;
	rv = rand() / (double)RAND_MAX;
	rv = rv - rand()%2;
	return rv;
}

void printMatrix(double m[], int n)
{
	int i, j;
	printf("\nMatrix Given");
	for(i = 0; i < n*n; i++)
	{
		j = getJ(i, n);
		if(j == 0)
			printf("\n");
		printf("%f ", m[i]);
	}
	printf("\n");
}

void matrixMatrix(double m1[], double m2[], double result[], int n)
{
	int i, j, k, x, y;
	double sum;
	sum = 0;
	for (i = 0; i < n; i++) 
	{
      		for (j = 0; j < n; j++) 
		{
       			for (k = 0; k < n; k++) 
			{
				x = matrixToArray(i, k, n);
				y = matrixToArray(k, j, n);
          			sum = sum + m1[x]*m2[y];
        		}
		x = matrixToArray(i, j, n);
		result[x] = sum;
		sum = 0;
		}
	}
}
int proc_map(int i, int size)
{
    	size = size - 1;
	int r;
	r = N/size;
    	int proc = i / r;
    	return proc + 1;
}

void initializeMatrices (double* m1, double* m2, double* result)
{
	int i, j; // Loop variables
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

	int rank, numprocs, sourceProcess;
	MPI_Comm GridComm;
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(numprocs == 4)
		gridSize = 2;
	else if(numprocs == 16)
		gridSize = 4;
	else if(numprocs == 64)
		gridSize = 8;
	else
		gridSize = 16;
	dimensions[0] = dimensions[1] = gridSize;
	periodic[0] = periodic[1] = 1;
//	MPI
	

	if(rank == 0)
	{
		m1 = (double*) malloc(sizeof(double)*N*N);
		m2 = (double*) malloc(sizeof(double)*N*N);
		result = (double*) malloc(sizeof(double)*N*N);
		initializeMatrices(m1, m2, result);
		printMatrix(m1, N);
		printMatrix(m2, N);
		t1 = MPI_Wtime();
		if(numprocs == 1)
			matrixMatrix(m1, m2, result, N);
		else
			printf("This is program part 1. It is only designed to work for 1 processor.\n");
		t2 = MPI_Wtime();
		elapsed = (t2 - t1) * 1000;
		printMatrix(result, N);
		printf("Time = %f milliseconds on one processor\n", elapsed);
	}
	MPI_Finalize();
	return 0;
}

