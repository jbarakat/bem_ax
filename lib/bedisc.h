/* BOUNDARY ELEMENT DISCRETIZATION
 *  Discretize the meridional arc into boundary elements (BEs),
 *  approximate boundary functions over each element, and
 *  approximate integrals over the elements using quadrature
 *  rules.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [Ch. 6]
 *  
 * PARAMETERS
 *  XG,YG  [input]		global coordinates of the BE end nodes
 *  XE,YE  [input]		global coordinates of the native BE nodes
 */

#ifndef BEDISC_H
#define BEDISC_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include "gauleg.h"
#include <math.h>
#include <gsl/gsl_sf_trig.h>


/* PROTOTYPES */
void be_straight();
void be_straight(double, double, double, double,
                 double, double&, double&);


/* IMPLEMENTATIONS */
/***********************************************************************/
// straight line segments
void be_straight(double x1, double y1, double x2, double y2,
                 double xi, double &x, double &y){
	// check that -1 <= xi <= 1
	if (xi > 1 || xi < -1){
			printf("Error: xi not in [-1, 1].");
			return;
	}
	else
		x = 0.5*(x2 + x1) + 0.5*(x2 - x1)*xi;
		y = 0.5*(y2 + y1) + 0.5*(y2 - y1)*xi;
}


/* Determine native element nodes:
 *  M = 0	- uniform elements
 *  M = 1 - linear basis
 *  M = 2 - quadratic basis
 */
void be_native(int M, int N, double *XG, double *YG, int n1, int n2){
	// check value of M
	if (M != 0 && M != 1 && M != 2){
		printf("Error: M must equal 0, 1, or 2");
		return;
	}
	
	// uniform elements
	if (M == 0){
		
	}
	
	// linear basis in Lagrange polynomials
	if (M == 1){

	}

	// quadratic basis in Lagrange polynomials
	if (M == 2){

	}
}


/***********************************************************************/

#endif
