/*This program contains a lot of unnecessary functions that were skeletons for the second part of this program.*/

#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 256
//converts a i and j into an array index
int matrixToArray(int i, int j)
{
	int k;
	k = i * N + j;
	return k;
}
//retrieves the i value from an array index
int getI(int k)
{
	int i;
	i = k / N;
	return i;
}
//retrieves the j index from an array index
int getJ(int k)
{
	int j;
	j = k%N;
	return j;
}
//returns a random number between -1 and 1
double r2()
{
	double rv;
	rv = rand() / (double)RAND_MAX;
	rv = rv - rand()%2;
	return rv;
}
//prints a matrix
void printMatrix(double m[N*N])
{
	int i, j;
	printf("\nMatrix Given");
	for(i = 0; i < N*N; i++)
	{
		j = getJ(i);
		if(j == 0)
			printf("\n");
		printf("%f ", m[i]);
	}
	printf("\n");
}
//performs linear matrix matrix multiplication
void matrixMatrix(double m1[N*N], double m2[N*N], double result[N*N])
{
	int i, j, k, x, y;
	double sum;
	sum = 0;
	for (i = 0; i < N; i++) 
	{
      		for (j = 0; j < N; j++) 
		{
       			for (k = 0; k < N; k++) 
			{
				x = matrixToArray(i, k);
				y = matrixToArray(k, j);
          			sum = sum + m1[x]*m2[y];
        		}
		x = matrixToArray(i, j);
		result[x] = sum;
		sum = 0;
		}
	}
}



int main(int argc, char* argv[])
{
	double t1, t2;
	int i, j, k;
	double elapsed;
	srand(1);

	int rank, numprocs, sourceProcess;
	MPI_Status stat;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank == 0)
	{
		double m1[N*N];
		double m2[N*N];
		double result[N*N];
		for( i = 0; i < N*N; i++)
		{
			m1[i] = r2();
			m2[i] = r2();
			result[i] = 0;
		}
	//	printMatrix(m1);
	//	printMatrix(m2);
		t1 = MPI_Wtime();
		if(numprocs == 1)
			matrixMatrix(m1, m2, result);
		else
			printf("This is program part 1. It is only designed to work for 1 processor.\n");
		t2 = MPI_Wtime();
		elapsed = (t2 - t1) * 1000;
		printMatrix(result);
		printf("Time = %f milliseconds on one processor\n", elapsed);
	}
	MPI_Finalize();
	return 0;
}

