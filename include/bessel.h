/* BESSEL FUNCTIONS
 *  Evaluate Bessel functions and modified Bessel functions.
 *
 * REFERENCES
 *  Abramowitz and Stegun, Dover Publications (1964)
 *
 * PARAMETERS
 *  x	[input]		field point
 *  n	[input]		integer order
 *  nu	[input]		fractional order
 *  J	[output]	Bessel function of the first kind
 *  Y	[output]	Bessel function of the second kind
 *  I	[output]	modified Bessel function of the first kind
 *  K	[output]	modified Bessel function of the second kind
 */

#ifndef BESSEL_H
#define BESSEL_H

/* HEADER FILES */
#include <lapacke.h>
#include <cblas.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>

/* PROTOTYPES */
double besselJ(int, double);
double besselJ(double, double);
void besselJArray(int, int, double, double*);

double besselY(int, double);
double besselY(double, double);
void besselYArray(int, int, double, double*);

double besselI(int, double);
double besselI(double, double);
void besselIArray(int, int, double, double*);

double besselK(int, double);
double besselK(double, double);
void besselKArray(int, int, double, double*);

#endif
