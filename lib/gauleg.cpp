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
	double *D, *E, *Z;
	int ldz = n;
	int i, j, info;
	
	// allocate memory
	D = (double*) calloc(n,sizeof(double));
	E = (double*) calloc(n-1,sizeof(double));
	Z = (double*) calloc(ldz*n,sizeof(double));

	// generate subdiagonal elements of positive definite tridiagonal matrix
	for (i = 0; i < n-1; i++){
		E[i] = 0.5/sqrt(1.0 - pow(2.0*(i + 1.0), -2.0));
	}
		
//	// print D and E
//	printf("\nD = ");
//	for (i = 0; i < n; i++)
//		printf("%.6f ",D[i]);
//	printf("\n");
//
//	printf("\nE = ");
//	for (i = 0; i < n-1; i++)
//		printf("%.6f ",E[i]);
//	printf("\n");
	
	// compute eigenvalues and eigenvectors
	info = LAPACKE_dstedc(LAPACK_ROW_MAJOR, 'I', n, D, E, Z, ldz);

//	// print D, E, and Z
//	printf("\nD = ");
//	for (i = 0; i < n; i++)
//		printf("%.6f ",D[i]);
//	printf("\n");
//
//	printf("\nE = ");
//	for (i = 0; i < n-1; i++)
//		printf("%.6f ",E[i]);
//	printf("\n");
//	
//	printf("\nZ = \n");
//	for (j = 0; j < n; j++){
//		for (i = 0; i < ldz; i++){
//			printf("%.6f ", Z[i + j*n]);
//		}
//		printf("\n");
//	}

	// determine abscissas and associated weights on the interval [-1, 1]
	for (i = 0; i < n; i++){
		X[i] = D[i];
		W[i] = 2*pow(Z[i], 2);
	}

//	// print X and W
//	printf("\nX = ");
//	for (i = 0; i < n; i++)
//		printf("%.6f ",X[i]);
//	printf("\n");
//
//	printf("\nW = ");
//	for (i = 0; i < n; i++)
//		printf("%.6f ",W[i]);
//	printf("\n");

	// release memory
	free(D);
	free(E);
	free(Z);
}

int main(){
	int n;
	double *X, *W;
	n = 5;

	// note: have to make sure the range of integration is from -1 to 1, see Numerical Recipes p. 207 (p. 184 on pdf)

	X = (double*) calloc(n,sizeof(double));
	W = (double*) calloc(n,sizeof(double));
	
	gauleg(n, X, W);

	// print abscissas and weights
	int i;
	printf("\nX = ");
	for (i = 0; i < n; i++)
		printf("%.6f ",X[i]);
	printf("\n");

	printf("\nW = ");
	for (i = 0; i < n; i++)
		printf("%.6f ",W[i]);
	printf("\n");

	return(0);
}
