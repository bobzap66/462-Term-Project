#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
using namespace std;

const int MATRIX_SIZE = 256;

int main() {
	int size, rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int i, j, k, l, m;
	int temp = 0;
	int destination = 0;
	int source = 0;
	int blockSize = (MATRIX_SIZE * MATRIX_SIZE) / size;
	int rowB = sqrt(blockSize);
	int columnB = rowB;
	double A[MATRIX_SIZE][MATRIX_SIZE];
	double B[MATRIX_SIZE][MATRIX_SIZE];
	double C[MATRIX_SIZE][MATRIX_SIZE];
	double a[rowB][columnB];
	double b[rowB][columnB];
	double c[rowB][columnB];
	double recv[rowB][columnB];
	double bigArray [(MATRIX_SIZE * MATRIX_SIZE)/blockSize][rowB][columnB];
	double bigBrray [(MATRIX_SIZE * MATRIX_SIZE)/blockSize][rowB][columnB];
	double bigCrray [(MATRIX_SIZE * MATRIX_SIZE)/blockSize][rowB][columnB];

	if(rank == 0) {
		for(i = 0; i < MATRIX_SIZE; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) {
			        double rv;
        			rv = rand() / (double)RAND_MAX;
        			rv = rv - rand()%2;
				A[i][j] = rv;

        			rv = rand() / (double)RAND_MAX;
        			rv = rv - rand()%2;
				B[i][j] = rv;
			}
		}

		for(i = 0; i < size; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) {
				for(k = 0; k < MATRIX_SIZE; k++) {
					temp = ((k / rowB)) + ((j / columnB) * sqrt(size)); 
					bigArray[temp][j-((j/columnB)*columnB)][k-((k/rowB)*rowB)] = A[j][k];
					bigBrray[temp][j-((j/columnB)*columnB)][k-((k/rowB)*rowB)] = B[j][k];
				}
			}
		}
	}	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatter(bigArray, blockSize, MPI_DOUBLE, a, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Scatter(bigBrray, blockSize, MPI_DOUBLE, b, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	for(i = 0; i < columnB; i++) {
		for(j = 0; j < rowB; j++) c[i][j] = 0;
	}

	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			for(k = 0; k < size; k++) {
				c[i][j] += a[i][k] * b[k][j]; 
			}
		}
	}

	for(m = 0; m < (MATRIX_SIZE / columnB) - 1; m++) {
		for(l = 0; l < MATRIX_SIZE / columnB; l++) {
			if(rank / (MATRIX_SIZE / columnB) == l) {
				if(l == 0) destination = size - ((MATRIX_SIZE / columnB) - rank);
				else destination = rank - (MATRIX_SIZE / columnB);

				MPI_Send(b, blockSize, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
			}
			else if(l == 0 && (rank / (MATRIX_SIZE / columnB)) == (MATRIX_SIZE / columnB) - 1) {
				source = rank - ((MATRIX_SIZE / columnB) * ((MATRIX_SIZE / columnB) - 1));			
			
				MPI_Recv(recv, blockSize, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else if(rank / (MATRIX_SIZE / columnB) == l - 1 ) {
				source = rank + (MATRIX_SIZE / columnB);
		
				MPI_Recv(recv, blockSize, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		for(i = 0; i < size; i++) {
			for(j = 0; j < size; j++) b[i][j] = recv[i][j];
		}

		for(l = 0; l < MATRIX_SIZE / columnB; l++) {
			if(rank % (MATRIX_SIZE / columnB) == l) {
				if(rank % (MATRIX_SIZE / columnB) == 0) destination = rank + ((MATRIX_SIZE / columnB) - 1);
				else destination = rank - 1;


				MPI_Send(a, blockSize, MPI_DOUBLE, destination, 0, MPI_COMM_WORLD);
			}
			else if(l == 0 && (rank % (MATRIX_SIZE / columnB)) == (MATRIX_SIZE / columnB) - 1) {
				source = rank - ((MATRIX_SIZE / columnB) - 1);			
			
				MPI_Recv(recv, blockSize, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
			else if(rank % (MATRIX_SIZE / columnB) == l - 1 ) {
				source = rank + 1;
		
				MPI_Recv(recv, blockSize, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}

		for(i = 0; i < size; i++) {
			for(j = 0; j < size; j++) a[i][j] = recv[i][j];
		}
			
		for(i = 0; i < size; i++) {
			for(j = 0; j < size; j++) {
				for(k = 0; k < size; k++) {
					c[i][j] += a[i][k] * b[k][j]; 
				}
			}
		}
	}

	MPI_Gather(c, blockSize, MPI_DOUBLE, bigCrray, blockSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if(rank == 0) {
		for(i = 0; i < size; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) {
				for(k = 0; k < MATRIX_SIZE; k++) {
					temp = ((k / rowB)) + ((j / columnB) * sqrt(size)); 
					C[j][k] = bigCrray[temp][j-((j/columnB)*columnB)][k-((k/rowB)*rowB)];
				}
			}
		}
	}

	MPI_Finalize();

	return 0;
}
