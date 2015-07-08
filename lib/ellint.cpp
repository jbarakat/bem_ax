/* ELLIPTIC INTEGRALS
 *  Evaluate incomplete and complete elliptic integrals.
 *
 * REFERENCES
 *  Abramowitz and Stegun, Dover Publications (1964)
 *
 * PARAMETERS
 *  k	[input]		elliptic modulus
 *  phi	[input]		amplitude
 *  n	[input]		characteristic
 *  F	[output]        incomplete elliptic integral of the first kind
 *  K	[output]        complete elliptic integral of the first kind
 *  E	[output]        (in)complete elliptic integral of the second kind
 *  Pi	[output]        (in)complete elliptic integral of the third kind
 */

/* HEADER FILES */
#include "ellint.h"

/* IMPLEMENTATIONS */
// Incomplete elliptic integral of the first kind
double ellintF(double phi, double k){
	const int mode = 0;
	double F;
	F = gsl_sf_ellint_F(phi, k, mode);
	return(F);
}


// Incomplete elliptic integral of the second kind

double ellintE(double phi, double k){
	const int mode = 0;
	double E;
	E = gsl_sf_ellint_E(phi, k, mode);
	return(E);
}


//// Incomplete elliptic integral of the third kind
//double ellintPi(double phi, double k, double n){
//	const int mode = 0;
//	double Pi;
//	Pi = gsl_sf_ellint_P(phi, k, n, mode);
//	return(Pi);
//}


// Complete elliptic integr4al of the first kind
double ellintK(double k){
	const int mode = 0;
	double K;
	K = gsl_sf_ellint_Kcomp(k, mode);
	return(K);
}


// Complete elliptic integral of the second kind
double ellintE(double k){
	const int mode = 0;
	double E;
	E = gsl_sf_ellint_Ecomp(k, mode);
	return(E);
}



//// Complete elliptic integral of the third kind
//double ellintPi(double k, double n){
//	const int mode = 0;
//	double Pi;
//	Pi = gsl_sf_ellint_Pcomp(k, n, mode);
//	return(Pi);
//}


