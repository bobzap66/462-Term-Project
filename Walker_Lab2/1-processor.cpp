/* This file generates a 128x128 matrix and a 128x1 vector and multiply them together, for part one a the project, unix time commands will 
 * be used to time this program.
 */

#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
using namespace std;

int main() {
	int size, rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int A[128][128];
	int x[128];
	int y[128];
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
	
		for(i = 0; i < 128; i++) {
			sum = 0;			

			for(j = 0; j < 128; j++) sum += (A[i][j] * x[j]);

			y[i] = sum;
		}
	}	
	
	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return 0;
}
