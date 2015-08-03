/* INTERPOLATION
 *  Generate interpolants using cubic splines or Lagrange polynomials.
 *
 * REFERENCES
 *  Moin, Cambridge University Press (2010) (Ch. 1)
 *  Pozrikidis, Chapman & Hall/CRC (2002) (Ch. 3)
 *  
 * PARAMETERS
 *  x,y    [input]		set of N+1 grid points
 *  xi,yi  [output]		interpolated grid point
 *  N      [input]		number of segments
 *  a      [output]		spline coefficient of 3rd derivative
 *  b      [output]		spline coefficient of 2nd derivative
 *  c      [output]		spline coefficient of 1st derivative
 *  L      [output]		Lagrange polynomial
 */

#ifndef INTERP_H
#define INTERP_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_trig.h>

/* PROTOTYPES */
void spline(const int, double*, double*,
            double, double, double, double,
						double*, double*, double*,
						double*, double*, double*);
//void spline(const int, double*, double*, double, double,
//            double, double&);
void spline(const int, double*, double*, double, double,
            double*, double*, double*);

void lagrange(const int, double*  , double , double*         );
void lagrange(const int, double*  , double , double*, double*);
void lagrange(const int, const int, double*, double*, double*);
void lagrange(const int, const int, double*, double*, double*, double*);
void lagrange(const int, double*  , double*, double , double&);

#endif
