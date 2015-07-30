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
 *  xf,rf  [input]		field points
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
	
	/*---------------------------------------------*/
	/*------------------- SETUP -------------------*/
	/*---------------------------------------------*/
	
	// declare variables
	int i, j, k, m, n;									// for loop indices
	int nelem, ngeom;										// number of geometric elements and nodes
	int nlocl, nglob;										// number of local and global basis nodes
	int ising;													// indicator for singular element
	double area, vlme;									// total area and volume
	double lamb;												// viscosity ratio
	
	double *xg, *rg;										// geometric points
	double *xc, *rc;										// collocation points
	
	double  ax,  bx, cx;								// spline coefficients
	double  ar,  br, cr;
	
	double *zquad, *wquad;							// quadrature abscissas and weights
	double  l    ,  dl   ;							// polygonal arc length and differential
	double  w    ,  J    ;							// local weight and Jacobian

	double *L    ,  Lq   ;							// Lagrange interpolating polynomials
	double *zlocl;											// local Gauss-Lobatto grid points

	double  xf  ,  rf  ;								// field point coordinates
	double  xq  ,  rq  ;								// source point coordinates
	double  dx  ,  dr  ;
	double  dxdl,  drdl, dsdl;
	double  dfx ,  dfr ;
	double  dfxf,  dfrf;

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
	double sum;
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
	A       = (double*) calloc( 4*nglob*nglob , sizeof(double));
	df      = (double*) malloc( 2*nglob       * sizeof(double));
	v       = (double*) malloc( 2*nglob       * sizeof(double));
	L       = (double*) malloc(   nlocl*nquad * sizeof(double));
  zlocl   = (double*) malloc(   nlocl       * sizeof(double));
  zquad   = (double*) malloc(   nquad       * sizeof(double));
  wquad   = (double*) malloc(   nquad       * sizeof(double));
  xg      = (double*) malloc(   ngeom       * sizeof(double));
  rg      = (double*) malloc(   ngeom       * sizeof(double));
  xc      = (double*) malloc(   nglob       * sizeof(double));
  rc      = (double*) malloc(   nglob       * sizeof(double));

	// get area and volume
	area  = Stokes.getArea();
	vlme  = Stokes.getVlme();

	// get viscosity ratio
	lamb  = Stokes.getVisc();

	// calculate Gauss-Lobatto points on the interval [-1,1]
	cf = M_PI/(nlocl-1);
	for (i = 0; i < nlocl; i++){
		zlocl[nlocl - i - 1] = gsl_sf_cos(cf*i);
	}
	
	/* get Gauss-Legendre abscissas and weights on the
	 * interval [-1,1] */
	gauleg(nquad, zquad, wquad);
	
	/* get Lagrange interpolating polynomials on the
	 * interval [-1,1] */
	lagrange(nlocl-1, nquad, zlocl, zquad, L);

	/*  NOTE: L[i,j] = L[nquad*j + i] is the
	 *  jth polynomial of the Gauss-Lobatto
	 *  grid evaluated at the ith quadrature
	 *  (Gauss-Legendre) grid point, defined
	 *  on the interval [-1,1], */

	/* calculate coordinates of collocation points
	 * (i.e., the global basis nodes */
	for (i = 0; i < nelem; i++){			// loop over global elements
		// get geometric parameters
		Stokes.getNode(i,   x0 , r0);
		Stokes.getNode(i+1, x1 , r1);
		Stokes.getPoly(i,   l0); 
		Stokes.getPoly(i+1, l1); 
	
		Stokes.getSpln(i,   ax , bx , cx ,
		                    ar , br , cr );
		
		lM = 0.5*(l1 + l0);
		lD = 0.5*(l1 - l0);

		for (j = 0; j < nlocl; j++){		// loop over local element nodes,
			if (j == nlocl-1 && i != nelem-1) // avoid double counting
				continue;

			/* get global basis node corresponding
			 * to local (i,j) basis */
			n = i*(nlocl - 1) + j;
	
			if (j == 0) { // if x0 coincides with geometric node
				xc[n] = x0;
				rc[n] = r0;
			}
			else if (j == nlocl-1 && i == nelem-1) { // last geometric node
				xc[n] = x1;
				rc[n] = r1;
			}
			else {
				// map Gauss-Lobatto point onto polygonal arc length
				l  = lM + lD*zlocl[j];
				dl = l  - l0;
				
				// interpolate to collocation point
				xc[n] = ((ax*dl + bx)*dl + cx)*dl + x0;
				rc[n] = ((ar*dl + br)*dl + cr)*dl + r0;
			}
		}
	}
	
	
	
	
	
	// NEED TO GET RADIUS OF TUBE!
	// (NEED A SMARTER WAY TO GET THIS)
	double rtube = 2.;
	






	/*---------------------------------------------*/
	/*- ASSEMBLE MATRIX OF INFLUENCE COEFFICIENTS -*/
	/*---------------------------------------------*/

	for (m = 0; m < nglob; m++){				/* ROWS: loop over global element
																			 * nodes (collocation points) */
		// get field point
		xf = xc[m];
		rf = rc[m];

		for (i = 0; i < nelem+1; i++){		/* COLUMNS: loop over boundary
		                                   * elements */
			// get geometric parameters
			Stokes.getNode(i,    x0, 	 r0);
			Stokes.getNode(i+1,  x1, 	 r1);
			Stokes.getPoly(i,    l0);
			Stokes.getPoly(i+1,  l1);
	
			Stokes.getSpln(i,    ax ,  bx , cx ,
			                     ar ,  br , cr );

			// check if xf,rf lies in the boundary element
			if (m >= i*(nlocl-1) && m <= (i+1)*(nlocl-1))
				ising = 1;
			else
				ising = 0;

			for (j = 0; j < nlocl; j++){		/* COLUMNS: loop over local element
			                                 * nodes */
				/* get global basis node corresponding
				 * to local (i,j) basis */
				n = i*(nlocl - 1) + j;
				
				
				
				
				
				// evaluate quadrature     //// THIS WILL CHANGE IF THERE'S A SINGULARITY ising == 1
				sum = 0.;
				for (k = 0; k < nquad; k++){	// loop over quadrature points
					// get Lagrange polynomial
					Lq = L[nquad*j + k];
					
					// map quadrature point onto polygonal arc length
					l  = lM + lD*zquad[k];
					dl = l  - l0;
		
					// interpolate to source point
					xq   = ((  ax*dl +    bx)*dl + cx)*dl + x0;
					rq   = ((  ar*dl +    br)*dl + cr)*dl + r0;
					dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
					drdl = (3.*ar*dl + 2.*br)*dl + cr;
					dsdl = sqrt(dxdl*dxdl + drdl*drdl);

					// define the Jacobian for the line integral
					J = lD*dsdl;

					// get weight for kth node
					w = J*wquad[k];
					
					
				
					/// FIRST THING TO DO IS JUST EVALUATE THE REGULAR MATRIX ELEMENTS,
					/// INCLUDING THE NODES SHARED BY TWO BEs, THEN
					/// EVALUATE THE SINGULAR ELEMENTS LATER


					// check for singularity
					// THIS WILL AFFECT THE CALCULATION OF THE GREEN'S FUNCTIONS
					
					//
					//  START FROM HERE
					//
					
					// evaluate Green's functions
					if (IGF == 0){				/* free-space stokeslet */
						gf_axR(xf, rf, xq, rq,
						       Mxx, Mxr, Mrx, Mrr);
					}
					else if (IGF == 1){		/* stokeslet bounded externally
																 * by a cylindrical tube */
						gf_axT(xf, rf, xq, rq, rtube,
						       Mxx, Mxr, Mrx, Mrr);
					}
				}
				
				// THIS IS WHERE WE AVOID DOUBLE COUNTING;
				// AFTER DOING ALL THAT WORK, DECIDE TO ADD
				// THE RESULT TO THE PREVIOUS ELEMENT OR TO
				// START A NEW ELEMENT

//				if (j == 0 && i != 0)
//					// ADD ADDITIONAL INTEGRAL (2 INSTEAD OF 1)
//					n;


				
		//		if (j == nlocl-1 && i != nelem-1) // avoid double counting
		//			continue;
				
				
			} // end of local element nodes
		} // end of boundary elements
	} // end of global element nodes

	/*  FOOTNOTE: A[m,n] corresponds to the matrix element 
	 *  associated with the integral evaluated at the mth
	 *  global boundary node over the density basis
	 *  function at the nth global basis node, where
	 *  n = i*(nlocl - 1) + j. */












	/*--------------------------------------------*/
	/*- single layer potential over the boundary -*/
	/*--------------------------------------------*/

	cf = -1./(8.*M_PI);



	// THE STUFF BELOW WILL BE DELETED... SUPPLANTED BY THE
	// ASSEMBLY OF THE MATRIX OF INFLUENCE COEFFICIENTS ABOVE
	
	/* loop over boundary nodes and evaluate integrals
	 * at field points (xf,rf) */
	for (i = 0; i < ngeom; i++){
		/* get position and jump in traction
		 * at field point */
		Stokes.getNode(i, xf  , rf  );
		Stokes.getTrct(i, dfxf, dfrf);

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
				// map quadrature point onto polygonal arc length
				l  = lM + lD*zquad[k];
				dl = l  - l0;
	
				// interpolate to source point
				xq   = ((  ax*dl +    bx)*dl + cx)*dl + x0;
				rq   = ((  ar*dl +    br)*dl + cr)*dl + r0;
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
					gf_axR(xf, rf, xq, rq,
					       Mxx, Mxr, Mrx, Mrr);
				}
				else if (IGF == 1){		/* stokeslet bounded externally
															 * by a cylindrical tube */
					gf_axT(xf, rf, xq, rq, rtube,
					       Mxx, Mxr, Mrx, Mrr);
				}

				






				
				
			} // end of Gauss-Legendre grid points
			
		} // end of boundary elements
		
		vx *= cf;
		vr *= cf;

		// set velocity for the boundary node

		//-------- NEED SET FUNCTION HERE -----------//

	} // end of boundary nodes (evaluation at field points)

//	// release memory
//	free(A);
//	free(df);   
//	free(v); 
//	free(L);  
//	free(zlocl);
//	free(zquad);
//	free(wquad);
//	free(xg);
//	free(rg); 
//	free(xc); 
//	free(rc); 
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
