/* QUADRATURE
 *  Integrate the Green's function over a spline element using Gauss-
 *  Legendre quadrature. The density functions are interpolated using
 *  Lagrange polynomials over each element.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [Ch. 6]
 *  
 * PARAMETERS
 *  x,r    [input]		source points
 *  xp,rp  [input]		field points
 *  nquad  [input]		number of quadrature points (per element)
 */

#ifndef QUAD_H
#define QUAD_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include "stokes.h"
#include "grnfcn.h"
#include "gauleg.h"
#include <gsl/gsl_sf_trig.h>

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Evaluate the single layer potential over an axisymmetric
 * contour.
 *  IGF = 0  free-space axisymmetric Green's function
 *  IGF = 1  Green's function bounded externally by a
 *            cylindrical tube
 */
void singleLayer(const int IGF, int nquad, stokes Stokes){
	if (IGF != 0 && IGF != 1){
		printf("Error: IGF can only take values of 0 or 1.");
		return;
	}
	
	// declare variables
	int i, j, k;												// for loop indices
	int nelem, ngeom;										// number of geometric elements and nodes
	int nlocl, nglob;										// number of local and global basis nodes
	double area, vlme;									// total area and volume
	double lamb;												// viscosity ratio
	
	double  ax,  bx, cx;								// spline coefficients
	double  ar,  br, cr;
	
	double *zquad, *wquad;							// quadrature abscissas and weights
	double  l  ,  dl ;									// polygonal arc length and differential
	double  w  ,  J  ;									// local weight and Jacobian

	double *L;													// Lagrange polynomials
	double *zlocl;											// local Gauss-Lobatto grid points

	double  xp  ,  rp  ;								// field point coordinates
	double  x   ,  r   ;								// source point coordinates
	double  dx  ,  dr  ;
	double  dxdl,  drdl, dsdl;
	double  dfx ,  dfr ;
	double  dfxp,  dfrp;

	double   x0,   r0;
	double   l0,   s0;
	double dfx0, dfr0;
	double  ks0,  kp0;
	double  tx0,  tr0;
	double  nx0,  nr0;
	
	double   x1,   r1;
	double   l1,   s1;
	double dfx1, dfr1;
	double  ks1,  kp1;
	double  tx1,  tr1;
	double  nx1,  nr1;

	double   lM,   lD;
	double dfxM, dfxD;
	double dfrM, dfrD;
	double  nxM,  nxD;
	double  nrM,  nrD;

	double *A, *df, *v;	// linear system: A*df = v
	double Mxx, Mxr, Mrx, Mrr;

	double vx, vr;
	double cf;
	
	/* get number of boundary elements,
	 * geometric nodes, and basis nodes */
	nelem = Stokes.getNElem();
	ngeom = nelem + 1;
	nlocl = Stokes.getNLocl();
	nglob = Stokes.getNGlob();

	// allocate memory
	A       = (double*) malloc( 4*nglob*nglob * sizeof(double));
	df      = (double*) malloc( 2*nglob       * sizeof(double));
	v       = (double*) malloc( 2*nglob       * sizeof(double));
	L       = (double*) malloc(   nlocl       * sizeof(double));
  zlocl   = (double*) malloc(   nlocl       * sizeof(double));
  zquad   = (double*) malloc(   nquad       * sizeof(double));
  wquad   = (double*) malloc(   nquad       * sizeof(double));

	// get area and volume
	area  = Stokes.getArea();
	vlme  = Stokes.getVlme();

	// get viscosity ratio
	lamb  = Stokes.getVisc();

	// assign Gauss-Lobatto points on the interval [-1,1]
	cf = M_PI/(nlocl-1);
	for (i = 0; i < nlocl; i++){
		zlocl[i] = gsl_sf_cos(cf*i);
	}
	
	
	
	
	// assign native nodes and build connectivity list
	
	
	
	
	
	
	
	// NEED TO GET RADIUS OF TUBE!
	double rc = 2.;
	






	// get Gauss-Legendre abscissas and weights
	gauleg(nquad, zquad, wquad);

	/*--------------------------------------------*/
	/*- single layer potential over the boundary -*/
	/*--------------------------------------------*/

	cf = -1./(8.*M_PI);

	/*- assemble the matrix of influence coefficients -*/
	
	// loop over boundary elements
	for (i = 0; i < nelem+1; i++){
		// loop over local element nodes
		for (j = 0; j < nlocl; j++){
			
			
		} // end of local element nodes
	} // end of boundary elements



	
	/* loop over boundary nodes and evaluate integrals
	 * at field points (xp,rp) */
	for (i = 0; i < ngeom; i++){
		/* get position and jump in traction
		 * at field point */
		Stokes.getNode(i, xp  , rp  );
		Stokes.getTrct(i, dfxp, dfrp);

		// initialize quadrature
		vx = 0.;
		vr = 0.;

		/* loop over boundary elements and carry out
		 * quadrature over source points (x,r) */
		for (j = 0; j < nelem; j++){
			// get parameters for the boundary element
			Stokes.getNode(j,     x0,  r0);
			Stokes.getNode(j+1,   x1,  r1);
			Stokes.getPoly(j,     l0);
			Stokes.getPoly(j+1,   l1);
	
			Stokes.getSpln(j,    ax ,  bx , cx ,
			                     ar ,  br , cr );
	
			Stokes.getNrml(j,    nx0,  nr0);
			Stokes.getNrml(j+1,  nx1,  nr1);

			Stokes.getTrct(j,   dfx0, dfr0);
			Stokes.getTrct(j+1, dfx1, dfr1);
			
			// prepare for quadrature
			  lM = 0.5*(  l1 +   l0);
			  lD = 0.5*(  l1 -   l0);

			dfxM = 0.5*(dfx1 + dfx0);
			dfxD = 0.5*(dfx1 - dfx0);

			dfrM = 0.5*(dfr1 + dfr0);
			dfrD = 0.5*(dfr1 - dfr0);
			  
			 nxM = 0.5*( nx1 +  nx0);
			 nxD = 0.5*( nx1 -  nx0);
	
			 nrM = 0.5*( nr1 +  nr0);
			 nrD = 0.5*( nr1 -  nr0);
	
			// loop over Gauss-Legendre grid points
			for (k = 0; k < nquad; k++){
				// map grid point onto polygonal arc length
				l  = lM + lD*zquad[k];
				dl = l  - l0;
	
				// interpolate to source point
				x    = ((  ax*dl +    bx)*dl + cx)*dl + x0;
				r    = ((  ar*dl +    br)*dl + cr)*dl + r0;
				dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
				drdl = (3.*ar*dl + 2.*br)*dl + cr;
				dsdl = sqrt(dxdl*dxdl + drdl*drdl);

				// define the Jacobian for the line integral
				J = lD*dsdl;

				// get weight for kth node
				w = J*wquad[k];



				// check for singularity
				/* ------- NEED A CHECK FOR i == j --------*/


				
				// evaluate Green's functions
				if (IGF == 0){				/* free-space stokeslet */
					gf_axR(xp, rp, x, r,
					       Mxx, Mxr, Mrx, Mrr);
				}
				else if (IGF == 1){		/* stokeslet bounded externally
															 * by a cylindrical tube */
					gf_axT(xp, rp, x, r, rc,
					       Mxx, Mxr, Mrx, Mrr);
				}

				






				
				
			} // end of Gauss-Legendre grid points
			
		} // end of boundary elements
		
		vx *= cf;
		vr *= cf;

		// set velocity for the boundary node

		//-------- NEED SET FUNCTION HERE -----------//

	} // end of boundary nodes (evaluation at field points)
}


				/* BASICALLY I HAVE SEVERAL OPTIONS HERE...
				 *  For fluid-fluid interfaces...
				 *  1. Free-space GF stokeslet + stresslet (if VISCRAT != 1)
				 *  2. Free-space GF stokeslet only (if VISCRAT == 1)
				 *  3. Tube GF stokeslet only for VISCRAT == 1
				 *  
				 *  For solid boundaries
				 *  1. Probably only stokeslet or stresslet? Figure out
				 *     whether I want to do a completed single / double layer
				 *     formulation here...
				 * 
				 * 
				 */
	
				




/* Determine native element nodes:
 *  M = 0	- uniform elements
 *  M = 1 - linear basis
 *  M = 2 - quadratic basis
 */
void be_native(int M, int N, double *XG, double *YG, int n1, int n2){
	// check value of M
	if (M != 0 && M != 1 && M != 2){
		printf("Error: M must equal 0, 1, or 2");
		return;
	}
	
	// uniform elements
	if (M == 0){
		
	}
	
	// linear basis in Lagrange polynomials
	if (M == 1){

	}

	// quadratic basis in Lagrange polynomials
	if (M == 2){

	}
}


#endif
