/* ELLIPTIC INTEGRALS
 *  Evaluate incomplete and complete elliptic integrals.
 *
 * REFERENCES
 *  Abramowitz and Stegun, Dover Publications (1964)
 *
 * PARAMETERS
 *  k   [input]         elliptic modulus
 *  phi [input]         amplitude
 *  n   [input]         characteristic
 *  F   [output]        incomplete elliptic integral of the first kind
 *  K   [output]        complete elliptic integral of the first kind
 *  E   [output]        (in)complete elliptic integral of the second kind
 *  Pi  [output]        (in)complete elliptic integral of the third kind
 */

#ifndef ELLINT_H
#define ELLINT_H

/* HEADER FILES */
#include <lapacke.h>
#include <cblas.h>
#include <gsl/gsl_sf_ellint.h>
#include <math.h>

/* PROTOTYPES */
double ellintF(double, double);
double ellintE(double, double);
//double ellintPi(double, double);
double ellintK(double);
double ellintE(double);
//double ellintPi(double);

#endif
