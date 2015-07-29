/* QUADRATURE
 *  Integrate the Green's function over a spline element using
 *  Gauss-Legendre quadrature.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [Ch. 6]
 *  
 * PARAMETERS
 *  x,r    [input]		source points
 *  xp,rp  [input]		field points
 *  ngl    [input]		number of quadrature points
 *  nfn    [input]		order of LAGRANGE INTERPOLANT...
 */

#ifndef QUAD_H
#define QUAD_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include "stokes.h"
#include "grnfcn.h"
#include "gauleg.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Evaluate the single layer potential over an axisymmetric
 * contour.
 *  IGF = 0  free-space axisymmetric Green's function
 *  IGF = 1  Green's function bounded externally by a
 *            cylindrical tube
 */
void singleLayer(const int IGF, int ngl, stokes Stokes){
	if (IGF != 0 && IGF != 1){
		printf("Error: IGF can only take values of 0 or 1.");
		return;
	}
	
	// declare variables
	int i, j, k, nnode, nelem;
	double A, V;
	double lamb;
	
	double  ax,  bx, cx;
	double  ar,  br, cr;
	
	double  x   ,  r   ;
	double  dx  ,  dr  ;
	double  dxdl,  drdl, dsdl;
	double  xp  ,  rp  ;
	double  dfx ,  dfr ;
	double  dfxp,  dfrp;

	double  x0,  r0;
	double  l0,  s0;
	double ks0, kp0;
	double tx0, tr0;
	double nx0, nr0;
	
	double  x1,  r1;
	double  l1,  s1;
	double ks1, kp1;
	double tx1, tr1;
	double nx1, nr1;

	double *zgl, *wgl;
	double  l  ,  dl ;
	double  w  ,  J  ;

	double   lM,   lD;
	double dfxM, dfxD;
	double dfrM, dfrD;
	double  nxM,  nxD;
	double  nrM,  nrD;

	double Mxx, Mxr, Mrx, Mrr;

	double vx, vr;
	double cf;
	
	// get number of boundary elements
	nelem = Stokes.getNElem();
	nnode = nelem + 1;

	// allocate memory
  zgl   = (double*) malloc( ngl * sizeof(double));
  wgl   = (double*) malloc( ngl * sizeof(double));

	// get area and volume
	A     = Stokes.getArea();
	V     = Stokes.getVlme();

	// get viscosity ratio
	lamb  = Stokes.getVisc();







	// NEED TO GET RADIUS OF TUBE!
	double rc = 2.;
	






	// get Gauss-Legendre abscissas and weights
	gauleg(ngl, zgl, wgl);

	/*--------------------------------------------*/
	/*- single layer potential over the boundary -*/
	/*--------------------------------------------*/

	cf = -1./(8.*M_PI);	
	
	/* loop over boundary nodes and evaluate integrals
	 * at field points (xp,rp) */
	for (i = 0; i < nnode; i++){
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
			for (k = 0; k < ngl; k++){
				// map grid point onto polygonal arc length
				l  = lM + lD*zgl[k];
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
				w = J*wgl[k];



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
	
				








/* Evaluate single layer potential over the nth spline element */

void singleLayer(int NGL, int n){

}




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
