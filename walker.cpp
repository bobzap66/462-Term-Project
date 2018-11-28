#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
using namespace std;

const int MATRIX_SIZE = 8;

int main() {
	int size, rank;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	int i, j, k;
	int temp = 0;
	int rowB = (MATRIX_SIZE * MATRIX_SIZE) / (size * size);
	int columnB = (MATRIX_SIZE * MATRIX_SIZE) / (size * size);
	int blockSize = rowB * columnB;
	int A[MATRIX_SIZE][MATRIX_SIZE];
	int B[MATRIX_SIZE][MATRIX_SIZE];
	int C[MATRIX_SIZE][MATRIX_SIZE];
	int a[rowB][columnB];
	int b[rowB][columnB];
	int c[rowB][columnB];
	int bigArray [(MATRIX_SIZE * MATRIX_SIZE)/blockSize][rowB][columnB];
	int bigBrray [(MATRIX_SIZE * MATRIX_SIZE)/blockSize][rowB][columnB];

	if(rank == 0) {
		for(i = 0; i < MATRIX_SIZE; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) {
				A[i][j] = rand() % 2;
				if(A[i][j] == 1 && (i + j) % 2 == 1) A[i][j] = -1;

				B[i][j] = rand() % 2;
				if(B[i][j] == 1 && (i + j) % 2 == 1) B[i][j] = -1;
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

		cout << "Num Processors: " << size << endl
		     << "Rows per Block: " << rowB << endl
		     << "Columnss per Block: " << columnB << endl
		     << "Block Size: " << blockSize << endl << endl;

		cout << "Matrix A - " << endl;
		for(i = 0; i < MATRIX_SIZE; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) cout << A[i][j] << " ";
			cout << endl;
		}
		cout << endl;

		for(i = 0; i < size; i++) {
			cout << "Submatrix of A for processor " << i  << " -" << endl;
			for(j = 0; j < rowB; j++) {
				for(k = 0; k < columnB; k++) cout << bigArray[i][j][k] << " ";
				cout << endl;
			}
			cout << endl;
		}

		cout << "Matrix B - " << endl;
		for(i = 0; i < MATRIX_SIZE; i++) {
			for(j = 0; j < MATRIX_SIZE; j++) cout << B[i][j] << " "; 
			cout << endl;	
		}
		cout << endl << endl;
		
		for(i = 0; i < size; i++) {
			cout << "Submatrix of B for processor " << i  << " -" << endl;
			for(j = 0; j < rowB; j++) {
				for(k = 0; k < columnB; k++) cout << bigBrray[i][j][k] << " ";
				cout << endl;
			}
			cout << endl;
		}
	}	

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatter(bigArray, blockSize, MPI_INT, a, blockSize, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(bigBrray, blockSize, MPI_INT, b, blockSize, MPI_INT, 0, MPI_COMM_WORLD);
	
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

	if(rank == 0) {
		cout << "On Processer " << rank << ", a * b =" << endl;
		for(i = 0; i < columnB; i++) {
			for(j = 0; j < rowB; j++) cout << c[i][j] << " ";
			cout << endl;
		}	
		cout << endl << endl;
	}

	MPI_Finalize();

	return 0;
}
