/* QUADRATURE
 *  Integrate the Green's function over a spline element using Gauss-
 *  Legendre quadrature. The density functions are interpolated using
 *  Lagrange polynomials over each element.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [Ch. 6]
 *  Pozrikidis, Chapman & Hall/CRC (2002) [Ch. 3]
 *  
 * PARAMETERS
 *  xq,rq  [input]		source points
 *  x ,r   [input]		field points
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
#include "gaulog.h"
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Evaluate the single layer potential over an axisymmetric
 * contour.
 *  IGF = 0  free-space axisymmetric Green's function
 *  IGF = 1  Green's function bounded externally by a
 *            cylindrical tube
 */
void singleLayer(const int IGF, int nquad, stokes Stokes, double* A){

	/* NOTE: A is a 2*nglob x 2*nglob matrix */
	
	if (IGF != 0 && IGF != 1){
		printf("Error: IGF can only take values of 0 or 1.\n");
		return;
	}

	if (nquad % 2 != 0){
		printf("Error: choose an even number of quadrature points to avoid singularities.\n");
		return;
	}
	
	if (A == NULL){
		printf("Error: no memory allocated for A.\n");
		return;
	}
	
	
	// NEED TO GET RADIUS OF TUBE!
	// (NEED A SMARTER WAY TO GET THIS)
	double rtube = 2.;
	





	
	/*---------------------------------------------*/
	/*------------------- SETUP -------------------*/
	/*---------------------------------------------*/
	
	// declare variables
	int i, j, k, m, n;									// for loop indices
	int nelem, ngeom;										// number of geometric elements and nodes
	int nlocl, nglob;										// number of local and global basis nodes

	double area, vlme;									// total area and volume
	double lamb;												// viscosity ratio
	
	double *xg, *rg;										// geometric points
	double *xc, *rc;										// collocation points
	
	double  ax,  bx, cx;								// spline coefficients
	double  ar,  br, cr;
	
	double *zquad, *wquad;							// quadrature abscissas and weights
	double  l    ,  dl   ;							// polygonal arc length and differential
	double  z    ,  w    ;							// local abscissa and weight
	double  h    ;											// metric coefficient
	
	double *L    ,  Lq   ;							// Lagrange interpolating polynomials
	double *zlocl;											// local Gauss-Lobatto grid points

	double  x   ,  r   ;								// field point coordinates
	double  xq  ,  rq  ;								// source point coordinates
	double  dx  ,  dr  ;
	double  dxdl,  drdl, dsdl;
	double  dfx ,  dfr ;

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

	double *df, *v;									// linear system: A*df = v

	double Axx , Axr , Arx , Arr ;		 	// block components of A
	double Mxx , Mxr , Mrx , Mrr ;			// total Green's function
	double MRxx, MRxr, MRrx, MRrr;			// free-space (axisymmetric) Green's function
	double MCxx, MCxr, MCrx, MCrr;			// complementary Green's function
	
	double  Ising;											// singular integral subtracted off
	double  I1   ,  I2   ,  I3 ,  I4 ;	// four integrals (sums to the singular integral)
	double *L1   , *L2   , *L3 , *L4 ;	// Lagrange polynomials for additional integrals
	double *z1   , *z2   , *z3 , *z4 ;	// parameters for additional integrals
	double  zsp  ,  zsm  ; 							// 1 + zsing and 1 - zsing
	double *zqsng, *wqsng; 							// singular quadrature abscissas and weights
	double  cf1  ,  cf2  ,  cf3,  cf4;	// coefficients of singular integrals

	int     ising;											// indicator for singular element
	int     nsing;											// number of quadrature points for singular kernels
	double  zsing;											// singular point on the [-1,1] interval
	double  lsing;											// singular point on the polygonal interval
	double  logZ ;

	double vx, vr;											// velocity components
	double cf;													// coefficient
	
	/* get number of boundary elements,
	 * geometric nodes, and basis nodes */
	nelem = Stokes.getNElem();
	ngeom = nelem + 1;
	nlocl = Stokes.getNLocl();
	nglob = Stokes.getNGlob();

	/* choose number of points for singular quadrature
	 * (maximum is 6 points) */
	nsing = 6;

	// allocate memory
	//A       = (double*) calloc( 4*nglob*nglob , sizeof(double));
	df      = (double*) malloc( 2*nglob       * sizeof(double));
	v       = (double*) malloc( 2*nglob       * sizeof(double));
  zquad   = (double*) malloc(   nquad       * sizeof(double));
  wquad   = (double*) malloc(   nquad       * sizeof(double));
  zqsng   = (double*) malloc(   nsing       * sizeof(double));
  wqsng   = (double*) malloc(   nsing       * sizeof(double));
	L       = (double*) malloc(   nlocl*nquad * sizeof(double));
	L1      = (double*) malloc(   nlocl*nquad * sizeof(double));
	L2      = (double*) malloc(   nlocl*nquad * sizeof(double));
	L3      = (double*) malloc(   nlocl*nsing * sizeof(double));
	L4      = (double*) malloc(   nlocl*nsing * sizeof(double));
	z1      = (double*) malloc(   nquad       * sizeof(double));
	z2      = (double*) malloc(   nquad       * sizeof(double));
	z3      = (double*) malloc(   nsing       * sizeof(double));
	z4      = (double*) malloc(   nsing       * sizeof(double));
  zlocl   = (double*) malloc(   nlocl       * sizeof(double));
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

	/* get Gauss-logarithmic abscissas and weights on
	 * the interval [0,1] */
	gaulog(nsing, zqsng, wqsng);
	
	/* get Lagrange interpolating polynomials on the
	 * interval [-1,1] */
	lagrange(nlocl-1, nquad, zlocl, zquad, L);

	/*  NOTE: L[i,j] = L[nquad*i + j] is the
	 *  ith polynomial of the Gauss-Lobatto
	 *  grid evaluated at the jth quadrature
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
	


	/*---------------------------------------------*/
	/*- ASSEMBLE MATRIX OF INFLUENCE COEFFICIENTS -*/
	/*---------------------------------------------*/

	for (m = 0; m < nglob; m++){				/* ROWS: loop over global element
																			 * nodes (collocation points) */
		// get field point
		x  = xc[m];
		r  = rc[m];

		for (i = 0; i < nelem; i++){			/* COLUMNS: loop over boundary
		                                   * elements */
			// get geometric parameters
			Stokes.getNode(i,    x0, 	 r0);
			Stokes.getNode(i+1,  x1, 	 r1);
			Stokes.getPoly(i,    l0);
			Stokes.getPoly(i+1,  l1);
	
			Stokes.getSpln(i,    ax ,  bx , cx ,
			                     ar ,  br , cr );
			
			lM = 0.5*(l1 + l0);
			lD = 0.5*(l1 - l0);

			// check if (x,r) is in the boundary element
			if (m >= i*(nlocl-1) && m <= (i+1)*(nlocl-1)){
				ising = 1; // singular element

				/* When the element is singular, the diagonal components
				 * of the free-space Green's function contain a logarithmic
				 * singularity:
				 * 
				 *  Mxx ~ Mrr ~ -2*log{[x - xq)^2 + (r - rq)^2]^(1/2)}
				 *
				 * This singularity is subtracted off inside
				 * the kernel and the remainder is evaluated separately.
				 *
				 * The singular integral is split up into four separate
				 * integrals, two of which are regular (use Gauss-Legendre
				 * quadrature on the [-1,1] interval) and the remaining
				 * two are singular (use a special logarithmic quadrature
				 * on the [0,1] interval). For reference, see Pozrikidis
				 * (2002), Ch. 3, pp. 71-73.
				 */

				// get singular point
				j     = m - (i*(nlocl-1));
				zsing = zlocl[j];
				lsing = lM + lD*zsing;

				/* prepare for singular quadrature
				 *  integrals 1,2 are regular
				 *  integrals 3,4 are singular */
				zsp   = 1 + zsing;
				zsm   = 1 - zsing;

				if (fabs(zsp) < 1e-8){
					cf1 = 0;
					cf3 = 0;
				}
				else {
					cf1 = -  zsp*gsl_sf_log(zsp);
					cf3 = -2*zsp;
				}

				if (fabs(zsm) < 1e-8) {
					cf2 = 0;
					cf4 = 0;
				}
				else {
					cf2 = -  zsm*gsl_sf_log(zsm);
					cf4 = -2*zsm;
				}

				for (k = 0; k < nquad; k++){ // regular quadrature
					z1[k] = 0.5*(-zsm + zsp*zquad[k]);
					z2[k] = 0.5*( zsp + zsm*zquad[k]);
				}

				for (k = 0; k < nsing; k++){ // singular quadrature
					z3[k] = zsing - zsp*zqsng[k];
					z4[k] = zsing + zsm*zqsng[k];
				}

				lagrange(nlocl-1, nquad, zlocl, z1, L1);
				lagrange(nlocl-1, nquad, zlocl, z2, L2);
				lagrange(nlocl-1, nsing, zlocl, z3, L3);
				lagrange(nlocl-1, nsing, zlocl, z4, L4);
			}
			else
				ising = 0; // non-singular element

			for (j = 0; j < nlocl; j++){		/* COLUMNS: loop over local element
			                                 * nodes */
				/* get global basis node corresponding
				 * to local (i,j) basis */
				n = i*(nlocl - 1) + j;

				/* NOTE: Nodes shared by adjacent elements
				 * have the same global index. The corresponding
				 * matrix elements are summed. */
				
				// prepare for quadrature
				Axx = 0;
				Axr = 0;
				Arx = 0;
				Arr = 0;
				
				// evaluate quadrature
				for (k = 0; k < nquad; k++){	// loop over quadrature points
					// get Lagrange interpolating polynomial
					Lq = L[nquad*j + k];
					
					// map quadrature point onto polygonal arc length
					z  = zquad[k];
					l  = lM + lD*z;
					dl = l  - l0;
		
					// interpolate to source point
					xq   = ((  ax*dl +    bx)*dl + cx)*dl + x0;
					rq   = ((  ar*dl +    br)*dl + cr)*dl + r0;
					dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
					drdl = (3.*ar*dl + 2.*br)*dl + cr;
					dsdl = sqrt(dxdl*dxdl + drdl*drdl);

					// define the metric coefficient for the line integral
					h = lD*dsdl;

					// get weight for kth node
					w = h*wquad[k];
					
					// evaluate free-space Green's function
					gf_axR(x   , r   , xq  , rq  ,
					       MRxx, MRxr, MRrx, MRrr);
					
				//	if (m == 0){
				//	double R2 = (x - xq)*(x - xq) + (r - rq)*(r - rq);
				//	double R  = sqrt(R2);
				//	double mxx = 2*M_PI*rq/R*(1 + (x - xq)*(x - xq)/R2);
				//	double mxr = - M_PI*(x - xq)/R*(1 - ((x - xq)*(x - xq) - rq*rq)/R2);
				//	//printf("%.14f %.4f\n", MRrx, mxr);
				//	}

					// evaluate complementary Green's function
					if (IGF == 1){				/* stokeslet bounded externally
					                       * by a cylindrical tube */
						gf_axC(x   , r   , xq  , rq  , rtube,
						       MCxx, MCxr, MCrx, MCrr);
					}
					else {
						MCxx = 0;
						MCxr = 0;
						MCrx = 0;
						MCrr = 0;
					}

					/* regularize diagonal components of the free-space
					 * Green's function if near the singular point */
					if (ising == 1){
						logZ  = gsl_sf_log(fabs(z - zsing));

						MRxx += 2*logZ;
						MRrr += 2*logZ;

					}

					// calculate total Green's function
					Mxx = MRxx + MCxx;
					Mxr = MRxr + MCxr;
					Mrx = MRrx + MCrx;
					Mrr = MRrr + MCrr;
					
					// increment the integrals
					Axx = w*Mxx*Lq;
					Axr = w*Mxr*Lq;
					Arx = w*Mrx*Lq;
					Arr = w*Mrr*Lq;

				} // end of quadrature points

				/* add back singular part of the integral
				 * using logarithmic quadrature */
				if (ising == 1){
					// prepare for quadrature
					I1 = 0; // regular integral
					I2 = 0; // regular integral
					I3 = 0; // singular integral
					I4 = 0; // singular integral

					// evaluate regular quadrature
					for (k = 0; k < nquad; k++){ // loop over quadrature points
						/*--------- first integral, I1 ----------*/
						if (fabs(cf1) > 1e-8){
							// get Lagrange interpolating polynomial
							Lq = L1[nquad*j + k];
							
							// get metric coefficient
							z  = z1[k];
							l  = lM + lD*z;
							dl = l  - l0;
		
							dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
							drdl = (3.*ar*dl + 2.*br)*dl + cr;
							dsdl = sqrt(dxdl*dxdl + drdl*drdl);

							h = lD*dsdl;

							// get weight for kth node
							w = h*wquad[k];

							// increment integral
							I1 += w*Lq;
						}
						
						/*-------- second integral, I2 ----------*/
						if (fabs(cf2) > 1e-8){
							// get Lagrange interpolating polynomial
							Lq = L2[nquad*j + k];

							// get metric coefficient
							z  = z2[k];
							l  = lM + lD*z;
							dl = l  - l0;
		
							dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
							drdl = (3.*ar*dl + 2.*br)*dl + cr;
							dsdl = sqrt(dxdl*dxdl + drdl*drdl);

							h = lD*dsdl;

							// get weight for kth node
							w = h*wquad[k];

							// increment integral
							I2 += w*Lq;
						}
					} // end of regular quadrature points

					// evaluate singular quadrature
					for (k = 0; k < nsing; k++){ // loop over quadrature points
						// get singular kernel
						z    = zqsng[k];
						logZ = gsl_sf_log(z);

						/*--------- third integral, I3 ----------*/
						if (fabs(cf3) > 1e-8){
							// get Lagrange interpolating polynomial
							Lq = L3[nsing*j + k];
							
							// get metric coefficient
							z  = z3[k];
							l  = lM + lD*z;
							dl = l  - l0;
		
							dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
							drdl = (3.*ar*dl + 2.*br)*dl + cr;
							dsdl = sqrt(dxdl*dxdl + drdl*drdl);

							h = lD*dsdl;

							// get weight for kth node
							w = h*wquad[k];

							// increment integral
							I3 += w*logZ*Lq;
						}
						
						/*-------- second integral, I2 ----------*/
						if (fabs(cf4) > 1e-8){
							// get Lagrange interpolating polynomial
							Lq = L4[nsing*j + k];

							// get metric coefficient
							z  = z4[k];
							l  = lM + lD*z;
							dl = l  - l0;
		
							dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
							drdl = (3.*ar*dl + 2.*br)*dl + cr;
							dsdl = sqrt(dxdl*dxdl + drdl*drdl);

							h = lD*dsdl;

							// get weight for kth node
							w = h*wquad[k];

							// increment integral
							I4 += w*logZ*Lq;
						}
					} // end of singular quadrature points
					
					I1 *= cf1;
					I2 *= cf2;
					I3 *= cf3;
					I4 *= cf4;

					Ising = I1 + I2 + I3 + I4;

					Axx += Ising;
					Arr += Ising;

				} // end of if statement

				// add 2x2 block of A
				A[2*nglob*(2*m)   + (2*n  )] += Axx;
				A[2*nglob*(2*m)   + (2*n+1)] += Axr;
				A[2*nglob*(2*m+1) + (2*n  )] += Arx;
				A[2*nglob*(2*m+1) + (2*n+1)] += Arr;
				
			} // end of local element nodes
		} // end of boundary elements
	} // end of global element nodes

	/*  FOOTNOTE: A[m,n] corresponds to the matrix element 
	 *  associated with the integral evaluated at the mth
	 *  global boundary node over the density basis
	 *  function at the nth global basis node, where
	 *  n = i*(nlocl-1) + j. */



	printf("\n got to the finish line!\n\n");





	/*--------------------------------------------*/
	/*- single layer potential over the boundary -*/
	/*--------------------------------------------*/

	cf = -1./(8.*M_PI);



//	// THE STUFF BELOW WILL BE DELETED... SUPPLANTED BY THE
//	// ASSEMBLY OF THE MATRIX OF INFLUENCE COEFFICIENTS ABOVE
//	
//	/* loop over boundary nodes and evaluate integrals
//	 * at field points (x,r) */
//	for (i = 0; i < ngeom; i++){
//		/* get position and jump in traction
//		 * at field point */
//		Stokes.getNode(i, x   , r   );
//		Stokes.getTrct(i, dfxf, dfrf);
//
//		// initialize quadrature
//		vx = 0.;
//		vr = 0.;
//
//		/* loop over boundary elements and carry out
//		 * quadrature over source points (x,r) */
//		for (j = 0; j < nelem; j++){
//			// get parameters for the boundary element
//			Stokes.getNode(j,     x0,  r0);
//			Stokes.getNode(j+1,   x1,  r1);
//			Stokes.getPoly(j,     l0);
//			Stokes.getPoly(j+1,   l1);
//	
//			Stokes.getSpln(j,    ax ,  bx , cx ,
//			                     ar ,  br , cr );
//	
//			Stokes.getNrml(j,    nx0,  nr0);
//			Stokes.getNrml(j+1,  nx1,  nr1);
//
//			Stokes.getTrct(j,   dfx0, dfr0);
//			Stokes.getTrct(j+1, dfx1, dfr1);
//			
//			// prepare for quadrature
//			  lM = 0.5*(  l1 +   l0);
//			  lD = 0.5*(  l1 -   l0);
//
//			dfxM = 0.5*(dfx1 + dfx0);
//			dfxD = 0.5*(dfx1 - dfx0);
//
//			dfrM = 0.5*(dfr1 + dfr0);
//			dfrD = 0.5*(dfr1 - dfr0);
//			  
//			 nxM = 0.5*( nx1 +  nx0);
//			 nxD = 0.5*( nx1 -  nx0);
//	
//			 nrM = 0.5*( nr1 +  nr0);
//			 nrD = 0.5*( nr1 -  nr0);
//	
//			// loop over Gauss-Legendre grid points
//			for (k = 0; k < nquad; k++){
//				// map quadrature point onto polygonal arc length
//				l  = lM + lD*zquad[k];
//				dl = l  - l0;
//	
//				// interpolate to source point
//				xq   = ((  ax*dl +    bx)*dl + cx)*dl + x0;
//				rq   = ((  ar*dl +    br)*dl + cr)*dl + r0;
//				dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
//				drdl = (3.*ar*dl + 2.*br)*dl + cr;
//				dsdl = sqrt(dxdl*dxdl + drdl*drdl);
//
//				// define the metric coefficient for the line integral
//				h = lD*dsdl;
//
//				// get weight for kth node
//				w = h*wquad[k];
//
//
//
//				// check for singularity
//				/* ------- NEED A CHECK FOR i == j --------*/
//
//
//				
//				// evaluate Green's functions
//				if (IGF == 0){				/* free-space stokeslet */
//					gf_axR(x , rf, xq, rq,
//					       Mxx, Mxr, Mrx, Mrr);
//				}
//				else if (IGF == 1){		/* stokeslet bounded externally
//															 * by a cylindrical tube */
//					gf_axT(x , rf, xq, rq, rtube,
//					       Mxx, Mxr, Mrx, Mrr);
//				}
//
//				
//
//
//
//
//
//
//				
//				
//			} // end of Gauss-Legendre grid points
//			
//		} // end of boundary elements
//		
//		vx *= cf;
//		vr *= cf;
//
//		// set velocity for the boundary node
//
//		//-------- NEED SET FUNCTION HERE -----------//
//
//	} // end of boundary nodes (evaluation at field points)

	// release memory
	free(df);   
	free(v); 
	free(L);  
	free(zlocl);
	free(zquad);
	free(wquad);
	free(xg);
	free(rg); 
	free(xc); 
	free(rc); 
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
