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
 *
 * NOTE
 *  GSL provides a basic implementation of Bessel functions, but
 *  does not support Bessel functions with complex arguments. The
 *  latter is handled by SLATEC subroutines.
 */

#ifndef BESSEL_H
#define BESSEL_H

/* HEADER FILES */
#include <lapacke.h>
#include <cblas.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>

/* COMPLEX DATAT TYPE */
#ifndef lapack_double_complex
#define lapack_double_complex double complex
#endif

/* SLATEC FORTRAN SUBROUTINES (NETLIB/AMOS LIBRARY) */
extern "C" void dbesj_(double*, double*, int*, double*, int*);
extern "C" void dbesy_(double*, double*, int*, double*);
extern "C" void dbesi_(double*, double*, int*, int*, double*, int*);
extern "C" void dbesk_(double*, double*, int*, int*, double*, int*);

extern "C" void zbesj_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern "C" void zbesy_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern "C" void zbesi_(double*, double*, double*, int*, int*, double*, double*, int*, int*);
extern "C" void zbesk_(double*, double*, double*, int*, int*, double*, double*, int*, int*);

/* PROTOTYPES */
double besselJ(int, double);
double besselJ(double, double);
double complex besselJ(int, double complex);
double complex besselJ(double, double complex);
void besselJArray(int, int, double, double*);
void besselJArray(int, int, double complex, double complex*);

double besselY(int, double);
double besselY(double, double);
double complex besselY(int, double complex);
double complex besselY(double, double complex);
void besselYArray(int, int, double, double*);
void besselYArray(int, int, double complex, double complex*);

double besselI(int, double);
double besselI(double, double);
double complex besselI(int, double complex);
double complex besselI(double, double complex);
void besselIArray(int, int, double, double*);
void besselIArray(int, int, double complex, double complex*);

double besselK(int, double);
double besselK(double, double);
double complex besselK(int, double complex);
double complex besselK(double, double complex);
void besselKArray(int, int, double, double*);
void besselKArray(int, int, double complex, double complex*);

#endif
