/* BESSEL FUNCTIONS
 *  Evaluate Bessel functions and modified Bessel functions.
 *
 * REFERENCES
 *  Abramowitz & Stegun, Dover Publications (1964)
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

/* HEADER FILES */
#include <lapacke.h>
#include <cblas.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>

/* PROTOTYPES */
double besselJ(int, double);
double besselJ(double, double);
double* besselJArray(int, int, double);

double besselY(int, double);
double besselY(double, double);
double* besselYArray(int, int, double);

double besselI(int, double);
double besselI(double, double);
double* besselIArray(int, int, double);

double besselK(int, double);
double besselK(double, double);
double* besselKArray(int, int, double);
