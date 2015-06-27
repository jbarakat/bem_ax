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
 *  J	[output]        Bessel function of the first kind
 *  Y	[output]        Bessel function of the second kind
 *  I	[output]        modified Bessel function of the first kind
 *  K	[output]        modified Bessel function of the second kind
 */

/* HEADER FILES */
#include "bessel.h"

/* IMPLEMENTATIONS */
// Bessel functions of the first kind
double besselJ(int n, double x){
	double Jn;
	
	if (n == 0)
		Jn = gsl_sf_bessel_J0(x);
	else if (n == 1)
		Jn = gsl_sf_bessel_J1(x);
	else
		Jn = gsl_sf_bessel_Jn(n, x);
	
	return(Jn);
}

double besselJ(double nu, double x){
	double Jnu;
	
	Jnu = gsl_sf_bessel_Jnu(nu, x);

	return(Jnu);
}

double* besselJArray(int nmin, int nmax, double x){
	double *Jn;
	int info;
	
	info = gsl_sf_bessel_Jn_array(nmin, nmax, x, Jn);

	return(Jn);	
}

// Bessel functions of the second kind
double besselY(int n, double x){
	double Yn;
	
	if (n == 0)
		Yn = gsl_sf_bessel_Y0(x);
	else if (n == 1)
		Yn = gsl_sf_bessel_Y1(x);
	else
		Yn = gsl_sf_bessel_Yn(n, x);
	
	return(Yn);
}

double besselY(double nu, double x){
	double Ynu;
	
	Ynu = gsl_sf_bessel_Ynu(nu, x);

	return(Ynu);
}

double* besselYArray(int nmin, int nmax, double x){
	double *Yn;
	int info;
	
	info = gsl_sf_bessel_Yn_array(nmin, nmax, x, Yn);

	return(Yn);	
}

// Modified Bessel functions of the first kind
double besselI(int n, double x){
	double In;
	
	if (n == 0)
		In = gsl_sf_bessel_I0(x);
	else if (n == 1)
		In = gsl_sf_bessel_I1(x);
	else
		In = gsl_sf_bessel_In(n, x);
	
	return(In);
}

double besselI(double nu, double x){
	double Inu;
	
	Inu = gsl_sf_bessel_Inu(nu, x);

	return(Inu);
}

double* besselIArray(int nmin, int nmax, double x){
	double *In;
	int info;
	
	info = gsl_sf_bessel_In_array(nmin, nmax, x, In);

	return(In);	
}

// Modified bessel functions of the second kind
double besselK(int n, double x){
	double Kn;
	
	if (n == 0)
		Kn = gsl_sf_bessel_K0(x);
	else if (n == 1)
		Kn = gsl_sf_bessel_K1(x);
	else
		Kn = gsl_sf_bessel_Kn(n, x);
	
	return(Kn);
}

double besselK(double nu, double x){
	double Knu;
	
	Knu = gsl_sf_bessel_Knu(nu, x);

	return(Knu);
}

double* besselKArray(int nmin, int nmax, double x){
	double *Kn;
	int info;
	
	info = gsl_sf_bessel_Kn_array(nmin, nmax, x, Kn);

	return(Kn);	
}
