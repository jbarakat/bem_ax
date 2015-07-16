/* GREEN'S FUNCTIONS
 *  Evaluate Green's functions for Stokes flow.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [pp. 89-91]
 *  Tozeren, Inter. J. Num. Meth. Fluids 4, 159-170 (1984) 
 *  
 * PARAMETERS
 */

#ifndef GRNFCN_H
#define GRNFCN_H

/* HEADER FILES */
#include "bessel.h"
#include "ellint.h"
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

/* PROTOTYPES */
// Green's function for a Stokeslet in a circular tube


/* IMPLEMENTATIONS */
/* Green's function M evaluated at (x,r) due to a ring of point forces
 * located at (x0,r0) in an unbound domain (free space).
 */
void gf_axR(double x, double r, double x0, double r0,
            double &Mxx, double &Mxr, double &Mrx, double &Mrr){
	// declare variables
	double K, E;
	double x2 = x*x;
	double x02 = x0*x0;
	double r2 = r*r;
	double r02 = r0*r0;
	double X = x - x0; 
	double X2 = X*X;
	double R2 = X2 + (r - r0)*(r - r0);
	double R = pow(R2, 0.5);
	double k2 = 4*r*r0/(X2 + (r + r0)*(r + r0));
	double k = pow(k2, 0.5);

	// evaluate complete elliptic integrals
	K = ellintK(k);
	E = ellintE(k);
	
	// calculate components of the Green's function
	Mxx = 2*k*pow(r0/r, 0.5)*(K + (X2/R2)*E);
	Mxr = -k*(X/pow(r0*r, 0.5))*(K - (r2 - r02 + X2)*E/R2);
	Mrx = k*(X/r)*pow(r0/r, 0.5)*(K + (r2 - r02 - X2)*E/R2);
	Mrr = (k/(r0*r))*pow(r0/r, 0.5)*((r02 + r2 + 2*X2)*K - (2*X2*X2 + 3*X2*(r02 + r2) + (r2 - r02)*(r2 - r02))*E/R2);

}

/* velocity u evaluated at (x,r) due to a ring of point forces 
 * of strength f located at (x0,r0) in an unbound domain (free space).
 */
void gf_axR_vel(double x, double r, double x0, double r0,
                double fx, double fr, double &ux, double &ur){
	// declare variables
	double Mxx, Mxr, Mrx, Mrr;
	
	// calculate components of the Green's function
	gf_axR(x, r, x0, r0, Mxx, Mxr, Mrx, Mrr);

	// calculate velocity components
	ux = (Mxx*fx + Mxr*fr)/(8*M_PI);
	ur = (Mrx*fx + Mrr*fr)/(8*M_PI);

}

/* Green's function M = MR + MC evaluated at (x,r) due to a ring of 
 * point forces located at (x0,r0), bounded externally by a cylindrical
 * tube of radius rc.
 *
 * NOTE: The lines commented out comprise the code for evaluating the
 * free-space Green's function, MR, in terms of Fourier integrals.
 * This method is not effective if r < r0, for which the components of
 * B diverge like exp[(r0 - r)*t] and the Fourier integrals become
 * improper. Instead, the closed-form expression for MR in terms of
 * complete elliptic integrals is used (see gf_axR above).
 */
void gf_axT(double x, double r, double x0, double r0, double rc,
            double &Mxx, double &Mxr, double &Mrx, double &Mrr){
	// declare variables
	int it;
	double t, s;
	double X = x - x0; 
	double omeg, omeg0, omegc, omegn;
	double I0, I1, I00, I10, I0c, I1c, I0n, I1n;
	double K0, K1, K00, K10, K0c, K1c, K0n, K1n;
	double Axx, Axr, Arx, Arr;
//	double Bxx, Bxr, Brx, Brr;
	double Bxxc, Bxrc, Brxc, Brrc;
	double Lxx, Lxr, Lrx, Lrr;
	double Fxx, Fxr, Frx, Frr;
	double MRxx, MRxr, MRrx, MRrr;
	double MCxx, MCxr, MCrx, MCrr;
//	double mRxx, mRxr, mRrx, mRrr;
	double mCxx, mCxr, mCrx, mCrr;
	double detL, fc;
	double Xt, cosXt, sinXt;
	double modB, modF;
	
	// integration parameters
	const int MAXIT = 1000000;
	const double TOL = 0.00001;
	double dt = 0.001;
	double ds = 0.0001;
	
	// initialize
//	mRxx = 0.;
//	mRxr = 0.;
//	mRrx = 0.;
//	mRrr = 0.;

	mCxx = 0.;
	mCxr = 0.;
	mCrx = 0.;
	mCrr = 0.;
	
	for (it = 0; it < MAXIT; it++){
		// update step size
		t = it*dt;
		s = it*ds;
		if (t == 0){
			t = 0.00001;
			s = 0.00001;
		}

	// // DELETE THIS LATER
	//	s  = 1;
	//	ds = dt; 
	
		// change of variable
		t = -gsl_sf_log(s);
		if (s >= 1)
			break;
		
		// scale radial coordinate
		omeg = t*r;
		omeg0 = t*r0;
		omegc = t*rc;
		omegn = 2*omegc - omeg - omeg0;

		// evaluate trigonometric functions
		Xt = X*t;
		cosXt = gsl_sf_cos(Xt);
		sinXt = gsl_sf_sin(Xt);
		
		// evaluate modified Bessel functions
		I0  = besselI(0, omeg);
		I1  = besselI(1, omeg);
		I00 = besselI(0, omeg0);
		I10 = besselI(1, omeg0);
		I0c = besselI(0, omegc);
		I1c = besselI(1, omegc);
		I0n = besselI(0, omegn);
		I1n = besselI(1, omegn);
		
		K0  = besselK(0, omeg);
		K1  = besselK(1, omeg);
		K00 = besselK(0, omeg0);
		K10 = besselK(1, omeg0);
		K0c = besselK(0, omegc);
		K1c = besselK(1, omegc);
		K0n = besselK(0, omegn);
		K1n = besselK(1, omegn);

		// calculate components of L matrix in La = 4b
		Lxx =     t*I0c;
		Lxr = omegc*I1c + 2*I0c;
		Lrx =     t*I1c;
		Lrr = omegc*I0c;
		detL = Lxx*Lrr - Lxr*Lrx;
		if (detL == 0)
			printf("Error: matrix L is singular.\n");
		fc = 4/detL; 

		// calculate components of B matrix
//		Bxx  = (-2*K0  +  omeg *K1) *I00 - omeg0*K0 *I10;
//		Bxr  =           -omeg *K0  *I00 + omeg0*K1 *I10;
//		Brx  =           -omeg *K1  *I10 + omeg0*K0 *I00;
//		Brr  = ( 2*K1  +  omeg *K0) *I10 - omeg0*K1 *I00;

		Bxxc = (-2*K0c +  omegc*K1c)*I00 - omeg0*K0c*I10;
		Bxrc =           -omegc*K0c *I00 + omeg0*K1c*I10;
		Brxc =           -omegc*K1c *I10 + omeg0*K0c*I00;
		Brrc = ( 2*K1c +  omegc*K0c)*I10 - omeg0*K1c*I00;
		
		// calculate components of A matrix from inverting La = 4b
		Axx = fc*( Lrr*Bxxc - Lxr*Bxrc);
		Axr = fc*(-Lrx*Bxxc + Lxx*Bxrc);
		Arx = fc*( Lrr*Brxc - Lxr*Brrc);
		Arr = fc*(-Lrx*Brxc + Lxx*Brrc);

		// calculate components of F matrix
		Fxx = t*I0*Axx + (omeg*I1 + 2*I0)*Axr;
		Fxr = t*I0*Arx + (omeg*I1 + 2*I0)*Arr;
		Frx = t*I1*Axx +  omeg*I0        *Axr;
		Frr = t*I1*Arx +  omeg*I0        *Arr;

		// regularize the xx integrand
		Fxx = Fxx + 8*K0n;

	// // DELETE THIS LATER
	//	if (it % 100 == 0)
	//		printf("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", t, s, Fxx, Fxr, Frx, Frr);

		// assign weight for numerical integration (extended trapezoidal rule)
		if (it == 0)
			fc = 0.5;
		else
			fc = 1.0;

		// increment the integrands
//		mRxx +=  fc*Bxx*cosXt;
//		mRxr +=  fc*Brx*sinXt; // note that mRxr is related to Brx
//		mRrx +=  fc*Bxr*sinXt; // and mRrx is related to Bxr
//		mRrr += -fc*Brr*cosXt;

		mCxx +=  fc*Fxx*cosXt/s;
		mCxr +=  fc*Fxr*sinXt/s;
		mCrx +=  fc*Frx*sinXt/s;
		mCrr += -fc*Frr*cosXt/s;

		// evaluate maximum modulus of B and F
		modB = 0.;
//		modB =            pow(pow(Bxx, 2), 0.5);
//		modB = fmax(modB, pow(pow(Bxr, 2), 0.5));
//		modB = fmax(modB, pow(pow(Brx, 2), 0.5));
//		modB = fmax(modB, pow(pow(Brr, 2), 0.5));

		modF =            pow(pow(Fxx, 2), 0.5);
		modF = fmax(modF, pow(pow(Fxr, 2), 0.5));
		modF = fmax(modF, pow(pow(Frx, 2), 0.5));
		modF = fmax(modF, pow(pow(Frr, 2), 0.5));

		// break loop below tolerance
		if (modB < TOL && modF < TOL)
			break;
}
	printf("%d\n", it);
	
	// calculate components of the free-space Green's function
	gf_axR(x, r, x0, r0, MRxx, MRxr, MRrx, MRrr);
//	MRxx = -4*r0*mRxx*dt;
//	MRxr = -4*r0*mRxr*dt;
//	MRrx = -4*r0*mRrx*dt;
//	MRrr = -4*r0*mRrr*dt;

	// calculate components of the complementary Green's function
	MCxx = r0*mCxx*ds - 4*M_PI*r0/pow(X*X + (2*rc - r - r0)*(2*rc - r - r0), 0.5);
	MCxr = r0*mCxr*ds;
	MCrx = r0*mCrx*ds;
	MCrr = r0*mCrr*ds;
	
	// calculate components of the total Green's function
	Mxx = MRxx + MCxx;
	Mxr = MRxr + MCxr;
	Mrx = MRrx + MCrx;
	Mrr = MRrr + MCrr;

}

/* velocity u evaluated at (x,r) due to a ring of point forces 
 * of strength f located at (x0,r0), bounded externally by a
 * cylindrical tube of radius rc.
 */
void gf_axT_vel(double x, double r, double x0, double r0, double rc,
                double fx, double fr, double &ux, double &ur){
	// declare variables
	double Mxx, Mxr, Mrx, Mrr;
	
	// calculate components of the Green's function
	gf_axT(x, r, x0, r0, rc, Mxx, Mxr, Mrx, Mrr);

	// calculate velocity components
	ux = (Mxx*fx + Mxr*fr)/(8*M_PI);
	ur = (Mrx*fx + Mrr*fr)/(8*M_PI);

}

/* AUXILIARY FUNCTIONS */

/***********************************************************************/

#endif
