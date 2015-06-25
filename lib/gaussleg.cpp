/* GAUSS-LEGENDRE QUADRATURE
 * Generate the abscissas and weigts for Gauss-Legendre quadrature by
 * solving an eigenvalue proglem for a symmetric, tridiagonal matrix.
 *
 * REFERENCES
 *  Golub and Welsch, Mathematics of Computation 23-106 (1969)
 * 
 * PARAMETERS
 *  n	[input]		number of points
 *  X   [output]	abcsissas
 *  W   [output]	weights
 */

// this will go into the interface file eventually...
//#ifndef GAUSSLEG_H
//#define GAUSSLEG_H

#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>

// function
void gaussleg(int n, double* X, double* W){
	// generate symmetric tridiagonal matrix

	// compute the eigenvalues and first component of the orthonormalized eigenvectors
}

int main(){
	int n;
	double * X;
	double * W;

	gaussleg(n, X, W);
}
