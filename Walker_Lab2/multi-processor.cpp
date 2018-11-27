/* This file generates a 128x128 matrix and a 128x1 vector and multiplies them together, for part two of the project (works for 4, 16, and 64 processors)
 * unix time commands will be used to time this program.
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

int main() {
	int size, rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int A[128][128];
	int x[128];
	int y[128];
	int partialA[128 / size][128];
	int partialY[128 / size];
	int i, j;
	int sum;

	if(rank == 0) {
		for(i = 0; i < 128; i++) {
			for(j = 0; j < 128; j++) {
				A[i][j] = rand() % 2;
				if(A[i][j] == 1 && (i + j) % 2 == 1) A[i][j] = -1;
			}

			x[i] = rand() % 2;
			if(x[i] == 1 && i % 2 == 1) x[i] = -1;
		}
	}	
	MPI_Bcast(x, 128, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(A, (128 / size) * 128, MPI_INT, partialA, (128 / size) * 128, MPI_INT, 0, MPI_COMM_WORLD);

	for(i = 0; i < 128 / size; i++) {
		sum = 0;
		for(j = 0; j < 128; j++) sum += partialA[i][j] * x[j];
		partialY[i] = sum;
	}

	MPI_Gather(partialY, 128 / size, MPI_INT, y, 128 / size, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
