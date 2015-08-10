/* RUN PROGRAM
 *  Set up and execute boundary element simulation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
//#include <boost/math/special_functions>
#include "interp.h"
#include "quad.h"
#include "surface.h"
#include "solver.h"
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif


void testSingleLayer();
void testQuadrature();
void testLagrange();
void testGeom();
void testInterp();

int main(){
	testSingleLayer();
//	testQuadrature();
//	testLagrange();
//	testGeom();
//	testInterp();

	return(0);
}

void testSingleLayer(){
	int i, j, k;
	int IGF = 0;
	int N = 2;
	int M = 1;
	int nquad;
	double lamb = 1;
	double *x, *r, *thet;
	double a, b;
	double *v;

	// interface parameters
	double gamm = 1;
	double ES = 0.;
	double ED = 0.;
	double EB = 0.;
	double ET = 0.;

	// constitutive model
	int model = 0;

	// get number of collocation points
	int ncoll = N*M + 1;

	// get number of quadrature points
	printf("nquad = ");
	scanf("%u", &nquad);

	// allocate memory
	x     = (double*) malloc((N+1)         * sizeof(double));
	r     = (double*) malloc((N+1)         * sizeof(double));
	thet  = (double*) malloc((N+1)         * sizeof(double));
	v     = (double*) malloc(2*ncoll       * sizeof(double));

	// define coordinates on an ellipse
	a = 1.;
	b = 1.;
	for (i = 0; i < N+1; i++){
		thet[i] = (N-i)*M_PI/N;
		x[i] = a*gsl_sf_cos(thet[i]);
		r[i] = b*gsl_sf_sin(thet[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	surface spheroid(model, N, M, lamb, gamm,
	                 ES, ED, EB, ET,
									 x, r);

	singleLayer(IGF, nquad, spheroid, v);

	free(x);
	free(r);
	free(thet);
	free(v);
}
}
