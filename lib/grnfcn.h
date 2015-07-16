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
void gf_axR(double, double, double, double,
            double&, double&, double&, double&);
void gf_axR_vel(double, double, double, double,
                double, double, double&, double&);
void gf_axT(double, double, double, double, double,
            double&, double&, double&, double&);
void gf_axT_ker(double, double, double, double, double, double,
            double &, double &, double&, double &);
void gf_axT_vel(double, double, double, double, double,
                double, double, double&, double&);

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
 * The numerical integrals are evaluated using an extended midpoint rule
 * with n levels of refinement.
 */
void gf_axT(double x, double r, double x0, double r0, double rc,
            double &Mxx, double &Mxr, double &Mrx, double &Mrr){
	// declare variables
	int i, j, k, m, n;
	double MRxx, MRxr, MRrx, MRrr;
	double MCxx, MCxr, MCrx, MCrr;
	double mCxx, mCxr, mCrx, mCrr;
	double dmCxx, dmCxr, dmCrx, dmCrr;
	double fc, dev;
	
	const int MAXIT = 1000;
	const double TOL = 0.0001;
	int np, nt;
	double s, ds;
	double smin = 0.;
	double smax = 1.;
	double Ds = smax - smin;
	
	// initialize
	mCxx = 0.;
	mCxr = 0.;
	mCrx = 0.;
	mCrr = 0.;

	dmCxx = 0.;
	dmCxr = 0.;
	dmCrx = 0.;
	dmCrr = 0.;
	
	s = 0.5*Ds;
	ds = Ds;
	np = 1;  // number of points added
	nt = 1;  // total number of points

	/* perform integration using the extended midpoint rule,
	 * tripling the number of integration points at each level
	 * of refinement (incremented by n).
	 */
	for (n = 0; n < MAXIT; n++){
		if (n == 0){
			// evaluate kernels at midpoint of the domain
			gf_axT_ker(x, r, x0, r0, rc, s, mCxx, mCxr, mCrx, mCrr);
			
			// calculate change in kernels
			dmCxx -= mCxx*ds; dmCxx = fabs(dmCxx);
			dmCxr -= mCxr*ds; dmCxr = fabs(dmCxr);
			dmCrx -= mCrx*ds; dmCrx = fabs(dmCrx);
			dmCrr -= mCrr*ds; dmCrr = fabs(dmCrr);

			// prepare for next stage of refinement
			np = 2;
			nt = 3;
		}
		else {
			// add 2 * 3^(n-1) additional points
			if (n != 1){
				np *= 3;
				nt *= 3;
			}
			
			// store kernel before adding new points
			dmCxx = mCxx*ds;
			dmCxr = mCxr*ds;
			dmCrx = mCrx*ds;
			dmCrr = mCrr*ds;

			// refine grid spacing
			ds /= 3.0;
			
			// evaluate kernels at additional grid points
			for (k = 0; k < nt; k++){
				if ((2*k + 1) % 3 != 0){ /* avoid double counting
																	  previous grid points */
					s = smin + (double(k) + 0.5)*ds;
					gf_axT_ker(x, r, x0, r0, rc, s, mCxx, mCxr, mCrx, mCrr);
				}
			}

			// calculate change in kernels
			dmCxx -= mCxx*ds; dmCxx = fabs(dmCxx);
			dmCxr -= mCxr*ds; dmCxr = fabs(dmCxr);
			dmCrx -= mCrx*ds; dmCrx = fabs(dmCrx);
			dmCrr -= mCrr*ds; dmCrr = fabs(dmCrr);

		}

		// calculate maximum deviation
		dev = fmax(dmCxx, dmCxr);
		dev = fmax(dev, dmCrx);
		dev = fmax(dev, dmCrr);

		// break loop below tolerance
		if (dev < TOL)
			break;
	}

	//printf("%d stages of refinement\n", n);
	//printf("%d total points\n", nt);
	
	// calculate components of the free-space Green's function
	gf_axR(x, r, x0, r0, MRxx, MRxr, MRrx, MRrr);

	// calculate components of the complementary Green's function
	MCxx = r0*mCxx*ds - 4*M_PI*r0/pow((x - x0)*(x - x0) + (2*rc - r - r0)*(2*rc - r - r0), 0.5);
	MCxr = r0*mCxr*ds;
	MCrx = r0*mCrx*ds;
	MCrr = r0*mCrr*ds;
	
	// calculate components of the total Green's function
	Mxx = MRxx + MCxx;
	Mxr = MRxr + MCxr;
	Mrx = MRrx + MCrx;
	Mrr = MRrr + MCrr;

}

/* Evaluate the kernel of the complementary Green's function, MC.
 *
 * NOTE: The lines commented out comprise the code for evaluating the
 * kernel of the free-space Green's function, MR, in terms of Fourier
 * integrals. This method is not effective if r < r0, for which the
 * components of B diverge like exp[(r0 - r)*t] and the integrals are
 * made improper. Instead, the closed-form expression for MR in terms
 * of complete elliptic integrals is used (see gf_axR above).
 */
void gf_axT_ker(double x, double r, double x0, double r0, double rc, double s,
            double &mCxx, double &mCxr, double &mCrx, double &mCrr){
	// declare variables
	double t;
	double X = x - x0;
	double Xt, cosXt, sinXt;
	double omeg, omeg0, omegc, omegn;
	double I0, I1, I00, I10, I0c, I1c, I0n, I1n;
	double K0, K1, K00, K10, K0c, K1c, K0n, K1n;
	double Axx, Axr, Arx, Arr;
//	double Bxx, Bxr, Brx, Brr;
	double Bxxc, Bxrc, Brxc, Brrc;
	double Lxx, Lxr, Lrx, Lrr;
	double Fxx, Fxr, Frx, Frr;
	double detL, fc;

	// change of integration coordinate
	t = -gsl_sf_log(s);
	
	// evaluate trigonometric functions
	Xt = X*t;
	cosXt = gsl_sf_cos(Xt);
	sinXt = gsl_sf_sin(Xt);
	
	// scale radial coordinate
	omeg = t*r;
	omeg0 = t*r0;
	omegc = t*rc;
	omegn = 2*omegc - omeg - omeg0;

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

//	// DELETE THIS LATER
//	printf("%.6f, %.6f, %.6f, %.6f, %.6f, %.6f\n", t, s, Fxx, Fxr, Frx, Frr);

	// increment the integrands
//	mRxx +=  Bxx*cosXt;
//	mRxr +=  Brx*sinXt; // note that mRxr is related to Brx
//	mRrx +=  Bxr*sinXt; // and mRrx is related to Bxr
//	mRrr += -Brr*cosXt;

	mCxx +=  Fxx*cosXt/s;
	mCxr +=  Fxr*sinXt/s;
	mCrx +=  Frx*sinXt/s;
	mCrr += -Frr*cosXt/s;

//	MRxx = -4*r0*mRxx*dt;
//	MRxr = -4*r0*mRxr*dt;
//	MRrx = -4*r0*mRrx*dt;
//	MRrr = -4*r0*mRrr*dt;

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
