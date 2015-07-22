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
#include <gsl/gsl_sf_trig.h>

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

void testLagrange();

int main(){
	testLagrange();

	return(0);
}

void testLagrange(){
	int i, j, k;
	double *X, *Y, *T;
	double x, y;
	double a, b;
	int N = 100;
	
	// allocate memory
	X = (double*) malloc(N * sizeof(double));
	Y = (double*) malloc(N * sizeof(double));
	T = (double*) malloc(N * sizeof(double));
	
	// define parameters of an ellipse
	a = 1.;
	b = 1.;
	for (i = 0; i < N; i++){
		T[i] = i*M_PI/(N-1);
		X[i] = a*gsl_sf_cos(T[i]);
		Y[i] = b*gsl_sf_sin(T[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	// interpolate 
	x = 0.6;
	be_lagrange(N, X, Y, x, y);
	printf("y = %.4f\n", y);

}
