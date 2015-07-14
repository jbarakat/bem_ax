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


//void gf_axR(double xmin, double xmax, const int Nx, double x0,
//            double rmin, double rmax, const int Nr, double r0, double * M){
//	// declare variables
//	int i, j, k;
//	double xarray[Nx], rarray[Nr];
//	double dx = (xmax - xmin)/Nx;
//	double dr = (rmax - rmin)/Nr;
//	double x, r, x2, r2, x02, r02;
//	double X, R, X2, R2, k, k2;
//	double Mxx, Mxr, Mrx, Mrr;
//	double K, E;
//
//	// x and r vectors (equidistant grid points)
//	for (i = 0; i < Nx; i++){
//		x = xmin + i*dx;
//		xarray[i] = x;
//	}
//
//	for (i = 0; i < Nr; i++){
//		r = rmin + i*dr;
//		rarray[i] = r;
//	}
//
//	// START FROM HERE
//	// 
//	double x2 = x*x;
//	double x02 = x0*x0;
//	double r2 = r*r;
//	double r02 = r0*r0;
//	double X = x - x0; 
//	double X2 = X*X;
//	double R2 = X2 + (r - r0)*(r - r0);
//	double R = pow(R2, 0.5);
//	double k2 = 4*r*r0/(X2 + (r + r0)*(r + r0));
//	double k = pow(k2, 0.5);
//
//	// complete elliptic integrals
//	K = ellintK(k);
//	E = ellintE(k);
//	
//	// tensor components of the Green's function
//	Mxx = 2*k*pow(r0/r, 0.5)*(K + (X2/R2)*E);
//	Mxr = -k*(X/pow(r0*r, 0.5))*(K - (r2 - r02 + X2)*E/R2);
//	Mrx = k*(X/r)*pow(r0/r, 0.5)*(K + (r2 - r02 + X2)*E/R2);
//	Mrr = (k/(r0*r))*pow(r0/r, 0.5)*((r02 + r2 + 2*X2)*K - (2*X2*X2 + 3*X2*(r02 + r2) + (r2 - r02)*(r2 - r02))*E/R2);
//
//	// assign tensor components using row-major ordering
//	M[0] = Mxx;
//	M[1] = Mxr;
//	M[2] = Mrx;
//	M[3] = Mrr;
//}



// complementary Green's function to satisfy the no-slip condition on tube walls
void gf_axC(){

}

/* AUXILIARY FUNCTIONS */

/***********************************************************************/

#endif
