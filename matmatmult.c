#include <stdlib.h>
#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define N 256

double r2()
{
	double rv;
	rv = rand() / (double)RAND_MAX;
	rv = rv - rand()%2;
	return rv;
}

void printMatrix(double m[N][N])
{
	int i, j;
	printf("\nMatrix Given\n");
	for(i = 0; i < N; i++)
	{
		for(j = 0; j < N; j++)
		{
			printf("%f ", m[i][j]);
		}
		printf("\n");
	}
}

void matrixMatrix(double m1[N][N], double m2[N][N], double result[N][N])
{
	int i, j, k;
	double sum;
	sum = 0;
	for (i = 0; i < N; i++) 
	{
      		for (j = 0; j < N; j++) 
		{
       			for (k = 0; k < N; k++) 
			{
          			sum = sum + m1[i][k]*m2[k][j];
        		}
		result[i][j] = sum;
		sum = 0;
		}
	}
}
int proc_map(int i, int size)
{
    	size = size - 1;
	double result;
	int r;
	result = N/size;
	if(fmod(result, 1) == 0)
		r = (int)result;
	else
		r = (int) result + 1;
    	int proc = i / r;
    	return proc + 1;
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
		double m1[N][N];
		double m2[N][N];
		double result[N][N];
		for( i = 0; i < N; i++)
		{
			for(j = 0; j < N; j++)
			{
				m1[i][j] = r2();
				m2[i][j] = r2();
				result[i][j] = 0;
			}
		}
		if(numprocs == 1)
			matrixMatrix(m1, m2, result);

	//printMatrix(m1);
	//printMatrix(m2);
		MPI_Barrier(MPI_COMM_WORLD);
		t1 = MPI_Wtime();
		for (i=1; i< numprocs; i++)
          	{
              		for (j =0; j < N; j++)
			{
				MPI_Send(m2[j], N, MPI_DOUBLE, i, 1000 + j, MPI_COMM_WORLD);
			}
		}
	
	
		for(i = 0; i < N; i++)
		{
			sourceProcess = proc_map(i, numprocs);
			MPI_Recv(result[i], N, MPI_DOUBLE, sourceProcess, i, MPI_COMM_WORLD, &stat);
		}	
		MPI_Barrier(MPI_COMM_WORLD);
		t2 = MPI_Wtime();
		elapsed = (t2 - t1)*1000;
		printMatrix(result);
		printf("Time = %f milliseconds on %d processor(s_\n", elapsed, numprocs);
	}
	else
	{
		double m2[N][N];
		for( i = 0; i < N; i++)
		{
			MPI_Recv(m2[i], N, MPI_DOUBLE, 0, 100 + i, MPI_COMM_WORLD, &stat);
		}
		double result[N];
		int processor = proc_map(i, numprocs);
		if(rank == processor)
		{
			double buffer[N];
			MPI_Recv(buffer, N, MPI_DOUBLE, 0, (100*(i+1)), MPI_COMM_WORLD, &stat);
			for(j = 0; j < N; j++)
			{
				double sum = 0; 
				for(k = 0; k < N; k++)
					sum = sum +(buffer[k]*m2[k][j]);
				result[j] = sum;
			}
			MPI_Send(result, N, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
		} 
	}
	MPI_Finalize();
	return 0;
}

