/*GAUSS-LEGENDRE QUADRATURE
 * Generate the abscissas and weigts for Gauss-Legendre quadrature by
 * solving an eigenvalue proglem for a symmetric, tridiagonal matrix.
 *
 * REFERENCES
 *  Golub and Welsch, Mathematics of Computation 23-106 (1969)
 * 
 * PARAMETERS
 *  n	[input]		number of points
 *  X	[output]	abcsissas
 *  W	[output]	weights
 */

// this will go into the interface file eventually...
//#ifndef GAUSSLEG_H
//#define GAUSSLEG_H

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

// function implementation
void gauleg(int n, double* X, double* W){
	// declare variables
	const double eps = 1.0e-14;
	double *D, *E, *Z, *T;
	int ldz = n;
	double beta, *tau;
	int i, j, info;
	
	// allocate memory
	X = (double*) calloc(n,sizeof(double));
	W = (double*) calloc(n,sizeof(double));
	D = (double*) calloc(n,sizeof(double));
	E = (double*) calloc(n-1,sizeof(double));
	Z = (double*) calloc(ldz*n,sizeof(double));
	T = (double*) calloc(n*n,sizeof(double));
	tau = (double*) calloc(n,sizeof(double));

	// generate subdiagonal elements of positive definite tridiagonal matrix
	for (i = 0; i < n-1; i++){
		beta = 0.5/sqrt(1.0 - pow(2.0*(i + 1.0), -2.0));
		E[i] = beta;
	}
		
	// generate symmetric positive definite tridiagonal matrix
	for (i = 0; i < n-1; i++){
		beta = 0.5/sqrt(1.0 - pow(2.0*(i + 1.0), -2.0));
		T[i*(n+1)+1] = beta;
		T[(i+1)*n+i] = beta;
	}

	// print D, E, and T
	printf("\nD = ");
	for (i = 0; i < n-1; i++)
		printf("%.6f ",D[i]);
	printf("\n");

	printf("\nE = ");
	for (i = 0; i < n-1; i++)
		printf("%.6f ",E[i]);
	printf("\n");
	
	printf("\nT = \n");
	for (j = 0; j < n; j++){
		for (i = 0; i < n; i++){
			printf("%.6f ", T[i + j*n]);
		}
		printf("\n");
	}

	// compute eigenvalues and first component of the orthonormalized eigenvectors
	info = LAPACKE_dstedc(LAPACK_ROW_MAJOR, 'I', n, D, E, Z, ldz);
	printf("%i",info);

	// print D, E, T, and Z
	printf("\nD = ");
	for (i = 0; i < n-1; i++)
		printf("%.6f ",D[i]);
	printf("\n");

	printf("\nE = ");
	for (i = 0; i < n-1; i++)
		printf("%.6f ",E[i]);
	printf("\n");
	
	printf("\nT = \n");
	for (j = 0; j < n; j++){
		for (i = 0; i < n; i++){
			printf("%.6f ", T[i + j*n]);
		}
		printf("\n");
	}

	printf("\nZ = \n");
	for (j = 0; j < n; j++){
		for (i = 0; i < ldz; i++){
			printf("%.6f ", Z[i + j*n]);
		}
		printf("\n");
	}
	// determine abscissas and associated weights */
}

int main(){
	int n;
	double * X;
	double * W;
	n = 5;

	// note: have to make sure the range of integration is from -1 to 1, see Numerical Recipes p. 207 (p. 184 on pdf)

	gauleg(n, X, W);
}
