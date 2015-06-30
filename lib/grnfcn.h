/* GREEN'S FUNCTIONS
 *  Evaluate Green's functions for Stokes flow.
 *
 * REFERENCES
 *  -- Green's function for a Stokeslet in a circular tube --
 *  Liron and Shahar, Journal of Fluid Mechanics 86-4 (1978)
 *  
 * PARAMETERS
 */

#ifndef GRNFCN_H
#define GRNFCN_H

/* HEADER FILES */
#include "bessel.h"
#include <gsl/gsl_sf_log.h>

/* PROTOTYPES */
// Green's function for a Stokeslet in a circular tube


/* IMPLEMENTATIONS */
/***********************************************************************/
/* Green's function for a Stokeslet in a circular tube
 *  NOMENCLATURE: LS78-X denotes Eq. (X) in Liron and Shahar (1978). */

/* MAIN FUNCTION */
void gf_tube(){
	//

}

/* AUXILIARY FUNCTIONS */
/* calculate the denominator of the harmonic functions Dk(s) (see LS78-5.1)
 * and its gradient wrt s, Dk'(s) == dDk/ds
 *  From LS78-5.1, 
 *  
 *   Dk(s) = s*Ik(s)*[I(k-1)(s)*I(k+1)(s)]' - 2*[s*Ik(s)]'*I(k-1)(s)*I(k+1)(s)
 *
 *  which can be rewritten, using derivative and recursion relations for
 *  modified Bessel functions, in the following form:
 *
 *   Dk(s) = s*[Ik(s)^2 - I(k-1)(s)*I(k+1)(s)]*[I(k-1)(s) + I(k+1)(s)]
 *            - 4*I(k-1)(s)*Ik(s)*I(k+1)(s)
 * 
 *  The gradient of Dk(s) is determined analytically to be
 *
 *   Dk'(s) = Ik(s)*{Ik(s) + s*[I(k-1)(s) + I(k+1)(s)]}
 *            + I(k-1)(s)*I(k+1)(s)*[(2 - k - 2*s)*I(k-1)(s)
 *                                  + (2 + k - 2*s)*I(k+1)(s)]
 *            - s*Ik(s)*[I(k+1)(s)^2 + I(k-1)(s)^2]
 *            + 4*Ik(s)*{I(k-1)(s)*I(k+1)(s) - s*[I(k-1)(s)*Ik(s) 
 *                      + I(k-1)(s)*I(k+1)(s) + Ik(s)I(k+1)(s)]}
 */
void calcDk(int k, double s, double &Dk, double &dDkds){
	// declare variables
	double Ikm1, Ik, Ikp1;
	double Dk1, Dk2;
	double dDkds1, dDkds2, dDkds3, dDkds4;
	double* Ikarray;

	// calculate modified Bessel functions of the first kind
	if (k <= -1){
		Ikarray = (double*) calloc(3,sizeof(double));
		besselIArray(-k-1, -k+1, s, Ikarray);
		Ikp1 = Ikarray[0];
		Ik   = Ikarray[1];
		Ikm1 = Ikarray[2];
	}
	else if (k == 0){
		Ikarray = (double*) calloc(2,sizeof(double));
		besselIArray(k, k+1, s, Ikarray);
		Ikp1 = Ikarray[1];
		Ik   = Ikarray[0];
		Ikm1 = Ikp1;
	}
	else {
		Ikarray = (double*) calloc(3,sizeof(double));
		besselIArray(k-1, k+1, s, Ikarray);
		Ikp1 = Ikarray[2];
		Ik   = Ikarray[1];
		Ikm1 = Ikarray[0];
	}

	// calculate Dk(s)
	Dk1 = s*(pow(Ik, 2) - Ikm1*Ikp1)*(Ikm1 + Ikp1);
	Dk2 = -4*Ikm1*Ik*Ikp1;
	Dk = Dk1 + Dk2;

	// calculate the gradient of Dk(s) wrt s
	dDkds1 = Ik*(Ik + s*(Ikm1 + Ikp1));
	dDkds2 = Ikm1*Ikp1*((2 - k - 2*s)*Ikm1 + (2 + k - 2*s)*Ikp1);
  dDkds3 = -s*Ik*(Ikp1*Ikp1 + Ikm1*Ikm1);
  dDkds4 = 4*Ik*(Ikm1*Ikp1 - s*(Ikm1*Ik + Ikm1*Ikp1 + Ik*Ikp1));
	dDkds = dDkds1 + dDkds2 + dDkds3 + dDkds4;

	// release memory
	free(Ikarray);
}


void calcDkArray(int kmin, int kmax, double s, double *Dk){
}

/* calculate the roots of Dk(s)
 *  The roots of Dk(s) include zero,	s = 0
 *  a complex sequence,								xn = an +i*bn,	n = 0, 1, ...
 *  and a purely imaginary sequence,	yn = i*cn,			n = 0, 1, ...
 *
 *  From LS78-6.4,6.5, the estimates of xn, yn are
 *
 *   xn ~ 0.5*log[(2*n + k + 1)*PI] + (2*n + k + 1)*0.5*PI*i
 *   yn ~ i*(0.25*PI + 0.5*k*PI + n*PI), n = 0, 1
 *  
 *  The exact roots may be computed using the Newton-Raphson method by
 *  inputting the function Dk(s) and its gradient Dk'(s):
 *
 *
 *
 *  The roots of Dk(s) in the first quadrant are determined by solving the
 *  equation (LS78-6.1),
 *
 *   0 = t*[Jk^2(t) - J(k-1)(t)*J(k+1)(t)]*[J(k-1)(t) - J(k+1)(t)] 
 *        - 4*J(k-1)(t)*J(k+1)(t)
 *
 *  where t = i*s. Denote the two sets of roots by sn = xn and s = yn,
 *  n = 0, 1, ..., where
 *
 *   xn = an + i*bn	are the complex roots
 *   yn = i*cn			are the imaginary roots
 *
 *  From LS78-6.4,6.5, the estimates of xn, yn are
 *
 *   xn ~ 0.5*log[(2*n + k + 1)*PI] + (2*n + k + 1)*0.5*PI*i, n = 0, 1, ...
 *   yn ~ i*(0.25*PI + 0.5*k*PI + n*PI), n = 0, 1, ...
 */
void calcDkRoots(int n, int k, double &an, double &bn, double &cn,
	lapack_complex_double *xn, lapack_complex_double *yn){
	// declare variables
	
	// calculate estimates of an, bn, cn
	an = 0.5*gsl_sf_log((2.0*n + k + 1.0)*M_PI);
	bn = (2.0*n + k + 1.0)*0.5*M_PI;
	cn = (0.25 + 0.5*k + n)*M_PI;

	//

}

/***********************************************************************/

#endif
