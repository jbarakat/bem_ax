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
#include <math.h>
#include <gsl/gsl_sf_log.h>

//typedef lapack_complex_float fcmplx;
//typedef lapack_complex_double dcmplx;

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

/* TEMPLATE */
template<class T> 

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
 *  This expression can be rewritten, using derivative and recursion relations
 *  for modified Bessel functions, in the following form:
 *
 *   Dk(s) = s*[Ik(s)^2 - I(k-1)(s)*I(k+1)(s)]*[I(k-1)(s) + I(k+1)(s)]
 *            - 4*I(k-1)(s)*Ik(s)*I(k+1)(s)
 * 
 *  The gradient of Dk(s) is determined analytically to be
 *
 *   Dk'(s) = 2*s*[Ik(s)^3 - I(k-1)(s)*Ik(s)*I(k+1)(s)]
 *            + Ik(s)^2*[(k - 4)*I(k-1)(s) - (k + 4)*I(k+1)(s)]
 *            - k*I(k-1)(s)*I(k+1)(s)*[I(k-1)(s) - I(k+1)(s)]
 *            + 8*s^(-1)*I(k-1)(s)*Ik(s)*I(k+1)(s)
 */
void calcDk(int k, double complex s, double complex &Dk, double complex &dDkds){
	// declare variables
	double complex Ikm1, Ik, Ikp1;
	double complex Dk1, Dk2;
	double complex dDkds1, dDkds2, dDkds3;
	double complex* Ikarray;
	
	// allocate memory
	Ikarray = (double complex*) calloc(3,sizeof(double complex));
	
	// evaluate modified Bessel functions
	besselIArray(k-1, k+1, s, Ikarray);
	Ikm1 = Ikarray[0]; 
	Ik   = Ikarray[1];
	Ikp1 = Ikarray[2];
	
//	printf("\nIkm1 = %.4f + i*%.4f\n", creal(Ikm1), cimag(Ikm1));
//	printf("\nIk = %.4f + i*%.4f\n", creal(Ik), cimag(Ik));
//	printf("\nIkp1 = %.4f + i*%.4f\n", creal(Ikp1), cimag(Ikp1));
	
	// calculate Dk(s)
	Dk1 = s*(Ik*Ik - Ikm1*Ikp1)*(Ikm1 + Ikp1);
	Dk2 = -4*Ikm1*Ik*Ikp1;
	Dk = Dk1 + Dk2;

	// calculate the gradient of Dk(s) wrt s
	dDkds1 = 2*s*(Ik*Ik*Ik - Ikm1*Ik*Ikp1);
	dDkds2 = Ik*Ik*((k - 4)*Ikm1 - (k + 4)*Ikp1) - k*Ikm1*Ikp1*(Ikm1 - Ikp1);
	dDkds3 = 8*Ikm1*Ik*Ikp1/s;
	
	dDkds = dDkds1 + dDkds2 + dDkds3;

	// release memory
	free(Ikarray);
}


void calcDkArray(int kmin, int kmax, double complex s, double complex *Dk, double complex *dDkds){
}

/* calculate the Dk(s) in the first quadrant (see LS78-6.1)
 * and its gradient wrt s, Dk'(s) == dDk/ds
 *  From LS78-6.1, 
 *  
 *   Dk(t) = t*[Jk^2(t) - J(k-1)(t)*J(k+1)(t)]*[J(k-1)(t) - J(k+1)(t)] 
 *        - 4*J(k-1)(t)*J(k+1)(t)
 *
 *  where t = i*s. The gradient of this expression is determined analytically 
 *  to be
 *
 *   Dk'(t) = 2*t*[-Jk(t)^3 + J(k-1)(t)*Jk(t)*J(k+1)(t)]
 *            + Jk(t)^2*[(k - 4)*J(k-1)(t) + (k + 4)*J(k+1)(t)]
 *            - k*J(k-1)(t)*J(k+1)(t)*[J(k-1)(t) + J(k+1)(t)]
 *            + 8*t^(-1)*J(k-1)(t)*Jk(t)*J(k+1)(t)
 */
void calcDk1Q(int k, double complex s, double complex &Dk, double complex &dDkds){
	// declare variables
	double complex Jkm1, Jk, Jkp1;
	double complex Dk1, Dk2;
	double complex dDkds1, dDkds2, dDkds3, dDkds4;
	double complex* Jkarray;
	
	// allocate memory
	Jkarray = (double complex*) calloc(3,sizeof(double complex));
	
	// evaluate modified Bessel functions
	besselJArray(k-1, k+1, s, Jkarray);
	Jkm1 = Jkarray[0];
	Jk   = Jkarray[1];
	Jkp1 = Jkarray[2];

	// calculate Dk(s)
	Dk1 = s*(Jk*Jk - Jkm1*Jkp1)*(Jkm1 - Jkp1);
	Dk2 = -4*Jkm1*Jk*Jkp1;
	Dk = Dk1 + Dk2;
	
	// calculate the gradient of Dk(s) wrt s
	dDkds1 = 2*s*(-Jk*Jk*Jk + Jkm1*Jk*Jkp1);
	dDkds2 = Jk*Jk*((k - 4)*Jkm1 + (k + 4)*Jkp1) - k*Jkm1*Jkp1*(Jkm1 + Jkp1);
	dDkds3 = 8*Jkm1*Jk*Jkp1/s;
	
	dDkds = dDkds1 + dDkds2 + dDkds3;

	// release memory
	free(Jkarray);
}


/* calculate the roots of Dk(s)
 *  The roots of Dk(s) in the first quadrant are determined by solving the
 *  equation (LS78-6.1),
 *
 *   0 = t*[Jk^2(t) - J(k-1)(t)*J(k+1)(t)]*[J(k-1)(t) - J(k+1)(t)] 
 *        - 4*J(k-1)(t)*J(k+1)(t)
 *
 *  where t = i*s. Clearly, s = 0 is a root. Additionally, let s = x denote
 *  a complex root and s = y denote a purely imaginary root. These additional
 *  roots comprise an infinite sequence:
 *
 *   xn = an + i*bn,   n = 0, 1, ...
 *   yn = i*cn,        n = 0, 1, ...
 *
 *  From LS78-6.4,6.5, the estimates of xn, yn are
 *
 *   xn ~ 0.5*log[(2*n + k + 1)*PI] + (2*n + k + 1)*0.5*PI*i
 *   yn ~ i*(0.25*PI + 0.5*k*PI + n*PI)
 *
 *  where n = 0, 1, ... The exact roots may be computed using the Newton-
 *  Raphson method by inputting the function Dk(s) and its gradient Dk'(s).
 */
void calcDkRoots(int n, int k, double &an, double &bn, double &cn, double complex *xn, double complex *yn){
	// declare variables
	const int MAXIT = 100;
	const double ACC = 1e-14;
	double complex t, dt, Dk, dDkdt;
	double tr, ti, dtr, dti, err, errr, erri; 
	
	// calculate estimates of an, bn, cn
	an = 0.5*gsl_sf_log((2*n + k + 1)*M_PI);
	bn = (2*n + k + 1)*0.5*M_PI;
	cn = (0.25 + 0.5*k + n)*M_PI;
	printf("an = %.8f\n",an);
	printf("bn = %.8f\n",bn);
	printf("cn = %.8f\n",cn);
	
	// Newton iteration for complex roots, x = a + i*b
	tr = bn;
	ti = -an;
	t = tr + I*ti;
	calcDk1Q(k, t, Dk, dDkdt);
	dt = -Dk/dDkdt;
	dtr = creal(dt);
	dti = cimag(dt);
	printf("dtr = %.4f\n",dtr);
	printf("dti = %.4f\n",dti);

	err = sqrt(dtr*dtr + dti*dti);
	printf("%.16f\n",err);
	
//	errr = sqrt(pow(dtr - ACC, 2));
//	erri = sqrt(pow(dti - ACC, 2));
//	printf("%.16f\n",errr);
//	printf("%.16f\n",erri);
	
	while (err > ACC){
		tr += dtr;
		ti += dti;
		t = tr + I*ti;
		calcDk1Q(k, t, Dk, dDkdt);
		dt = -Dk/dDkdt;
		dtr = creal(dt);
		dti = cimag(dt);
		err = sqrt(dtr*dtr + dti*dti);
		printf("%.16f\n",err);
//		errr = sqrt(pow(dtr - ACC, 2));
//		erri = sqrt(pow(dti - ACC, 2));
//		printf("%.4f\n",errr);
//		printf("%.4f\n",erri);
	}
	an = -ti;
	bn = tr;
	
	cn = 0.0;
	xn[0] = 0.0;
	yn[0] = 0.0;
	

}

/***********************************************************************/

#endif
