/* GREEN'S FUNCTIONS
 *  Evaluate Green's functions for Stokes flow.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992)
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
#include <gsl/gsl_sf_log.h>

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

/* TEMPLATE */
template<class T> 

/* PROTOTYPES */
// Green's function for a Stokeslet in a circular tube


/* IMPLEMENTATIONS */
/***********************************************************************/

/* MAIN FUNCTIONS */
// total Green's function for a ring of point forces in a tube
void gf_axT(){

}

/* Green's function corresponding to a ring of point forces
 *  x 	axial coordinate
 *  r	cylindrical radial coordinate
 *  R	spherical radial coordinate
 *  M	Green's function tensor
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

	// complete elliptic integrals
	K = ellintK(k);
	E = ellintE(k);
	
	// tensor components of the Green's function
	Mxx = 2*k*pow(r0/r, 0.5)*(K + (X2/R2)*E);
	Mxr = -k*(X/pow(r0*r, 0.5))*(K - (r2 - r02 + X2)*E/R2);
	Mrx = k*(X/r)*pow(r0/r, 0.5)*(K + (r2 - r02 - X2)*E/R2);
	Mrr = (k/(r0*r))*pow(r0/r, 0.5)*((r02 + r2 + 2*X2)*K - (2*X2*X2 + 3*X2*(r02 + r2) + (r2 - r02)*(r2 - r02))*E/R2);
}

/* velocity field at (x,r) due to a ring of point forces at (x0,r0)
 *  f	point force
 *  u	velocity
 */
void gf_axR_vel(double x, double r, double x0, double r0,
                double fx, double fr, double &ux, double &ur){
	// declare variables
	double Mxx, Mxr, Mrx, Mrr;
	
	// calculate Green's function components
	gf_axR(x, r, x0, r0, Mxx, Mxr, Mrx, Mrr);

	// calculate velocity components
	ux = (Mxx*fx + Mxr*fr)/(8*M_PI);
	ur = (Mrx*fx + Mrr*fr)/(8*M_PI);
}

// complementary Green's function to satisfy the no-slip condition on tube walls
void gf_axC(double x, double r, double x0, double r0, double rc,
            double &Mxx, double &Mxr, double &Mrx, double &Mrr){
	// declare variables
	int it, Nt;
	double t, dt;
	double X = x - x0; 
	double omeg, omeg0, omegc, omegn;
	double I0, I1, I00, I10, I0c, I1c, I0n, I1n;
	double K0, K1, K00, K10, K0c, K1c, K0n, K1n;
	double Axx, Axr, Arx, Arr;
	double Bxx, Bxr, Brx, Brr;
	double Lxx, Lxr, Lrx, Lrr;
	double Fxx, Fxr, Frx, Frr;
	double detL, fc;
	double cosXt, sinXt;
	
	// integration parameters (step size, upper bound)
	dt = 0.01;
	Nt = 1000;
	
	for (it = 0; it < Nt; it++){
		t = it*dt;
		if (t == 0)
			t = 0.00001;
		
		// scaled radius
		omeg = t*r;
		omeg0 = t*r0;
		omegc = t*rc;
		omegn = 2*omegc - omeg - omeg0;
		
		// Bessel functions
		I0 = besselI(0, omeg);
		I1 = besselI(1, omeg);
		I00 = besselI(0, omeg0);
		I10 = besselI(1, omeg0);
		I0c = besselI(0, omegc);
		I1c = besselI(1, omegc);
		I0n = besselI(0, omegn);
		I1n = besselI(1, omegn);
		
		K0 = besselK(0, omeg);
		K1 = besselK(1, omeg);
		K00 = besselK(0, omeg0);
		K10 = besselK(1, omeg0);
		K0c = besselK(0, omegc);
		K1c = besselK(1, omegc);
		K0n = besselK(0, omegn);
		K1n = besselK(1, omegn);

		// components of L matrix in La = 4b
		Lxx = t*I0c;
		Lxr = omegc*I1c + 2*I0c;
		Lrx = t*I1c;
		Lrr = omegc*I0c;
		detL = Lxx*Lrr - Lxr*Lrx;
		if (detL == 0)
			printf("Error: matrix L is singular.\n");
		fc = 4/detL; 

		// components of B matrix
		Bxx = (-2*K0c + omegc*K1c)*I00 - K0c*omeg0*I10;
		Bxr =               -omegc*K0c + K1c*omeg0*I10;
		Brx =           -omegc*K1c*I10 + K0c*omeg0*I00;
		Brr = ( 2*K1c + omegc*K0c)*I10 - K1c*omeg0*I00;

		// components of A matrix from inverting La = 4b
		Axx = fc*( Lrr*Bxx - Lxr*Bxr);
		Axr = fc*(-Lrx*Bxx + Lxx*Bxr);
		Arx = fc*( Lrr*Brx - Lxr*Brr);
		Arr = fc*(-Lrx*Brx + Lxx*Brr);

		// components of F matrix
		Fxx = t*I0*Axx + (omeg*I1 + 2*I0)*Axr;
		Fxr = t*I0*Arx + (omeg*I1 + 2*I0)*Arr;
		Frx = t*I1*Axr +  omeg*I0        *Axr;
		Frr = t*I1*Arx +  omeg*I0        *Arr;

		// regularize the xx integrand
		Fxx = Fxx + 8*K0n;


	}


	// evaluate components of B
	
}

/* AUXILIARY FUNCTIONS */

/***********************************************************************/

#endif
