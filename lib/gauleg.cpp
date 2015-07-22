/* GAUSS-LEGENDRE QUADRATURE
 *  Generate the abscissas and weigts for Gauss-Legendre quadrature by
 *  solving an eigenvalue proglem for a symmetric, tridiagonal matrix.
 *
 * REFERENCES
 *  Golub and Welsch, Mathematics of Computation 23-106 (1969)
 * 
 * PARAMETERS
 *  n	[input]		number of points
 *  X	[output]	abcsissas
 *  W	[output]	weights
 */

/* HEADER FILES */
#include "gauleg.h"

/* IMPLEMENTATIONS */
/* Generate abscissas X and weights W for an nth order
 * Gauss-Legendre quadrature on the interval [-1, 1] */
void gauleg(int n, double * X, double * W){
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
		
	// compute eigenvalues and eigenvectors
	info = LAPACKE_dstedc(LAPACK_ROW_MAJOR, 'I', n, D, E, Z, ldz);

	// determine abscissas and associated weights on the interval [-1, 1]
	for (i = 0; i < n; i++){
		X[i] = D[i];
		W[i] = 2*pow(Z[i], 2);
	}

	// release memory
	free(D);
	free(E);
	free(Z);
}

