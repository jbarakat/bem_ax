/* RUN PROGRAM
 *  Set up and execute boundary element simulation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "gauleg.h"
#include "bessel.h"
#include "ellint.h"
//#include <boost/math/special_functions>
#include "grnfcn.h"
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

void ellipse(int, double, double, double *, double *);

int main(){
	
  /*-----------------------------------------------------*/
  /*----------------------- SETUP -----------------------*/
  /*-----------------------------------------------------*/	

	// declare variables
	int     i,  j,  k;
	int     nelem, nnode;
	double  a,  b;
	double *x, *r;

	// get number of boundary elements
//	printf("Number of boundary elements = ");
//	scanf("%u\n",&nelem);
	nelem = 10;
	nnode = nelem + 1;

	// allocate memory
	x = (double*) malloc( nnode * sizeof(double));
	r = (double*) malloc( nnode * sizeof(double));


  /*-----------------------------------------------------*/
  /*----------------- INITIAL CONDITION -----------------*/
  /*-----------------------------------------------------*/	

	// choose major and minor radii
	a = 1.;
	b = 1.;

	// generate initial contour


	return(0);
}

/* Generate the coordinates of an ellipse */
void ellipse(int n, double a, double b, double *x, double *y){
	// declare variables
	double *thet;

	// allocate memory 
	thet = (double*) malloc((n+1) * sizeof(double));
	
	
	
}
