/* RUN PROGRAM
 *  Set up and execute boundary element simulation.
 */

#include <iostream>
//#include <cstdlib>
//#include <cstdio>
#include <complex.h>
//#include <boost/math/special_functions/bernoulli.hpp>
#include "solver.h"
#include <math.h>
#include <vector>
#include <string>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

void ellipse(int, double, double, vector<double> &, vector<double> &);

int main(){
	
  /*-----------------------------------------------------*/
  /*----------------------- SETUP -----------------------*/
  /*-----------------------------------------------------*/	

	// declare variables
	int             i,  j,  k;
	int             nelem, nnode, nlocl;
	int             nstep;
	int             nquad;
	int             model;
	double          dt;
	double          a,  b;
	vector<double>  x, r;

	double          lamb, gamm;
	double          ES, ED, EB, ET;

	string          opath = "output/";

  /*-----------------------------------------------------*/
  /*----------------- INITIAL CONDITION -----------------*/
  /*-----------------------------------------------------*/	

	// NOTE: IN THE FUTURE, ALL OF THE STUFF IN THIS SECTION
	// WILL BE PART OF AN INPUT FILE.

	// choose number of timesteps and size of timestep
	nstep = 5000000;
	dt    = 0.1;

	// choose number of quadrature points
	nquad = 6;

	// get number of boundary elements
	nelem = 20;
	nnode = nelem + 1;

	// get number of native elements
	nlocl = 2;

	// allocate memory (for performance)
	x.reserve(nnode);
	r.reserve(nnode);

	// choose major and minor radii
	a = 1.;
	b = 1.;
	
	// generate initial contour
	ellipse(nelem, a, b, x, r);
	
	// choose constitutive model and assign parameters
	model = 0;
	lamb  = 1.0;
	gamm  = 0.001;
	ES    = 0.0;
	ED    = 0.0;
	EB    = 0.0;
	ET    = 0.0;

	// initialize boundary
	surface drop(model, nelem, nlocl,
	             lamb, gamm, ES, ED, EB, ET, 
				       x.data(), r.data());

//	double *xp, *rp;
//	xp = &x[0];
//	rp = &r[0];
//	
//	surface drop(model, nelem, nlocl,
//	             lamb, gamm, ES, ED, EB, ET, 
//				       xp, rp);
	

  /*-----------------------------------------------------*/
  /*------------------- TIME EVOLUTION ------------------*/
  /*-----------------------------------------------------*/

	timeInt(nstep, nquad, dt, drop, opath);
	
	return(0);
}

/* Generate (n+1) nodes on an ellipse */
void ellipse(int n, double a, double b, vector<double> &x, vector<double> &y){
	// declare variables
	int    i;
	double thet;
	double xp, yp;

	// assign coordinates on ellipse
	for (i = 0; i < n+1; i++){
		thet = i*M_PI/n;
		xp = a*gsl_sf_cos(thet);
		yp = b*gsl_sf_sin(thet);

		x.push_back(xp);
		y.push_back(yp);
	}
}
