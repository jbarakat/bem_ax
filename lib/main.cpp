/* MAIN PROGRAM
 * Execute library functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "gauleg.h"
#include "bessel.h"
#include "ellint.h"
//#include <boost/math/special_functions>
#include "grnfcn.h"
#include "bedisc.h"
#include "interp.h"
#include <gsl/gsl_sf_trig.h>

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

void testInterp();

int main(){
	testInterp();

	return(0);
}

void testInterp(){
	int i, j, k;
	double *X, *Y, *T;
	double x, y;
	double a, b;
	double *A, *B, *C;
	int N = 10;
	
	// allocate memory
	X = (double*) malloc((N+1) * sizeof(double));
	Y = (double*) malloc((N+1) * sizeof(double));
	T = (double*) malloc((N+1) * sizeof(double));
	A = (double*) malloc( N    * sizeof(double));
	B = (double*) malloc((N+1) * sizeof(double));
	C = (double*) malloc( N    * sizeof(double));
	
	// define parameters of an ellipse
	a = 1.;
	b = 1.;
	for (i = 0; i < N+1; i++){
		T[i] = i*M_PI/N;
		X[i] = a*gsl_sf_cos(T[i]);
		Y[i] = b*gsl_sf_sin(T[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	// interpolate (Lagrange)
	x = 0.9;
	lagrange(N, X, Y, x, y);
	printf("y = %.4f\n", y);

	// interpolate (cubic spline)
	double slope1 = 0;
	double slope2 = 0;
	spline_clamped(N, X, Y, slope1, slope2, A, B, C);
	printf("A        B        C\n");
	for (i = 0; i < N+1; i++){
		printf("%.4f %.4f %.4f\n", A[i], B[i], C[i]);
	}
	
	free(X);
	free(Y);
	free(T);
	free(A);
	free(B);
	free(C);
}
