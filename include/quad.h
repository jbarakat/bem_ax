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
#include "surface.h"
#include "grnfcn.h"
#include "gauleg.h"
#include "gaulog.h"
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>




/* NOTES: Probably can streamline the evaluation of the singular integrals
 * (fewer if statements, evaluate the coefficients later rather than first,
 * etc.). */












/* PROTOTYPES */


/* IMPLEMENTATIONS */

/* Evaluate the single layer potential over an axisymmetric
 * contour.
 *  IGF = 0  free-space axisymmetric Green's function
 *  IGF = 1  Green's function bounded externally by a
 *            cylindrical tube
 */
void singleLayer(const int IGF, int nquad, surface Surface,
                 double *v, double *Df){
	
	if (IGF != 0 && IGF != 1){
		printf("Error: IGF can only take values of 0 or 1.\n");
		return;
	}

	if (nquad % 2 != 0){
		printf("Error: choose an even number of quadrature points to avoid coinciding with the collocation points.\n");
		return;
	}





	
	// NEED TO GET RADIUS OF TUBE!
	// (NEED A SMARTER WAY TO GET THIS)
	double rtube = 2.;





	
	/*---------------------------------------------*/
	/*------------------- SETUP -------------------*/
	/*---------------------------------------------*/
	
	// declare variables
	int    i, j, k, m, n;								// for loop indices
	int    nelem, ngeom;								// number of geometric elements and nodes
	int    nlocl, nglob;								// number of local and global basis nodes

	double area, vlme;									// total area and volume
	double lamb;												// viscosity ratio
	
	double *xc, *rc;										// collocation points
	double *ks, *kp;										// principal curvatures
	double *tx, *tr;										// tangent components
	double *nx, *nr;										// normal components
	
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

	double   x0,   r0;									// parameters for lower geometric node
	double   l0,   s0;
	double  ks0,  kp0;
	double  tx0,  tr0;
	double  nx0,  nr0;
	
	double   x1,   r1;									// parameters for upper geometric node
	double   l1,   s1;
	double  ks1,  kp1;
	double  tx1,  tr1;
	double  nx1,  nr1;

	double   lM,   lD;									// mean and deviator of arc length segment

	double Axx , Axr , Arx , Arr ;		 	// block components of A
	double Mxx , Mxr , Mrx , Mrr ;			// total Green's function
	double MRxx, MRxr, MRrx, MRrr;			// free-space (axisymmetric) Green's function
	double MCxx, MCxr, MCrx, MCrr;			// complementary Green's function
	
	double  Ising;											// singular integral subtracted off
	double  IR0  ,  IR1  ,  IS0,  IS1;	// four integrals (sums to the singular integral)
	double *LR0  , *LR1  , *LS0, *LS1;	// Lagrange polynomials for additional integrals
	double *zR0  , *zR1  , *zS0, *zS1;	// parameters for additional integrals
	double *lR0  , *lR1  , *lS0, *lS1;	// parameters for additional integrals
	double *zqsng, *wqsng; 							// singular quadrature abscissas and weights
	double  cfR0 ,  cfR1 , cfS0, cfS1;	// coefficients of singular integrals

	int     ising;											// indicator for singular element
	int     nsing;											// number of quadrature points for singular kernels
	double  zsing;											// singular point on the [-1,1] interval
	double  lsing;											// singular point on the polygonal interval
	double  lD0  ,  lD1  ;							/* positive deviation of singular point from
																			 * lower and upper bound */
	double  logz ;											// logarithms for singular quadrature
	double  logl ;

	double *A   ;												// single layer operator and density
	double *Dfx , *Dfr ;								// traction components
	double  vx  ,  vr  ;								// velocity components

	double  cf;													// coefficient
	
	/* get number of boundary elements,
	 * geometric nodes, and basis nodes */
	nelem = Surface.getNElem();
	ngeom = nelem + 1;
	nlocl = Surface.getNLocl();
	nglob = Surface.getNGlob();

	/* choose number of points for singular quadrature
	 * (refer to GAULOG.H) */
	nsing = 5;

	// allocate memory
	A       = (double*) calloc( 4*nglob*nglob , sizeof(double));

  xc      = (double*) malloc(   nglob       * sizeof(double));
  rc      = (double*) malloc(   nglob       * sizeof(double));
	Dfx     = (double*) malloc(   nglob       * sizeof(double));
	Dfr     = (double*) malloc(   nglob       * sizeof(double));
	ks      = (double*) malloc(   nglob       * sizeof(double));
	ks      = (double*) malloc(   nglob       * sizeof(double));
	kp      = (double*) malloc(   nglob       * sizeof(double));
	tx      = (double*) malloc(   nglob       * sizeof(double));
	tr      = (double*) malloc(   nglob       * sizeof(double));
	nx      = (double*) malloc(   nglob       * sizeof(double));
	nr      = (double*) malloc(   nglob       * sizeof(double));

  zlocl   = (double*) malloc(   nlocl       * sizeof(double));
  zquad   = (double*) malloc(   nquad       * sizeof(double));
  wquad   = (double*) malloc(   nquad       * sizeof(double));
	L       = (double*) malloc(   nlocl*nquad * sizeof(double));

  zqsng   = (double*) malloc(   nsing       * sizeof(double));
  wqsng   = (double*) malloc(   nsing       * sizeof(double));
	LR0     = (double*) malloc(   nlocl*nquad * sizeof(double));
	LR1     = (double*) malloc(   nlocl*nquad * sizeof(double));
	LS0     = (double*) malloc(   nlocl*nsing * sizeof(double));
	LS1     = (double*) malloc(   nlocl*nsing * sizeof(double));
	zR0     = (double*) malloc(   nquad       * sizeof(double));
	zR1     = (double*) malloc(   nquad       * sizeof(double));
	zS0     = (double*) malloc(   nsing       * sizeof(double));
	zS1     = (double*) malloc(   nsing       * sizeof(double));
	lR0     = (double*) malloc(   nquad       * sizeof(double));
	lR1     = (double*) malloc(   nquad       * sizeof(double));
	lS0     = (double*) malloc(   nsing       * sizeof(double));
	lS1     = (double*) malloc(   nsing       * sizeof(double));
	
	// get area and volume
	area  = Surface.getArea();
	vlme  = Surface.getVlme();

	// get viscosity ratio
	lamb  = Surface.getVisc();
	
	
	
	/*-----------------------------------------------*/
	// NOTE: THIS STUFF CAN BE MOVED OUTSIDE THE FUNCTION...
	// SHOULDN'T HAVE TO RECALCULATE THE FUNCTION INTERPOLATION
	// POINTS (GAUSS-LOBATTO GRID), THE GAUSS QUADRATURE
	// POINTS, THE SINGULAR QUADRATURE POINTS, OR THE
	// LAGRANGE POLYNOMIALS ON THE REGULAR QUADRATURE
	// POINTS EEEEVERY TIME WE SOLVE THE BIE!!


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

	/*-----------------------------------------------*/
	
	


	
	/* calculate coordinates of collocation points
	 * (i.e., the global basis nodes */
	for (i = 0; i < nelem; i++){			// loop over global elements
		// get geometric parameters
		Surface.getNode(i,   x0 , r0);
		Surface.getNode(i+1, x1 , r1);
		Surface.getPoly(i,   l0); 
		Surface.getPoly(i+1, l1); 
	
		Surface.getSpln(i,  ax , bx , cx ,
		                    ar , br , cr );
		
		lM = 0.5*(l1 + l0);
		lD = 0.5*(l1 - l0);

		for (j = 0; j < nlocl; j++){				// loop over local element nodes
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
			Surface.getNode(i  ,  x0 ,  r0);
			Surface.getNode(i+1,  x1 ,  r1);
			Surface.getPoly(i  ,  l0);
			Surface.getPoly(i+1,  l1);
			
			Surface.getSpln(i  ,  ax ,  bx , cx ,
			                      ar ,  br , cr );
			
			lM = 0.5*(l1 + l0);
			lD = 0.5*(l1 - l0);

			// check if (x,r) is in the boundary element
			if (m >= i*(nlocl-1) && m <= (i+1)*(nlocl-1)){
				ising = 1; // singular element

				/* As (x,r) --> (xq,rq), the diagonal components of the
				 * free-space Green's function become logarithmically
				 * singular:
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
				 *  integrals R0,R1 are regular
				 *  integrals S0,S1 are singular */
				lD0   = lsing - l0   ;
				lD1   = l1    - lsing;

				cfR0 = 0;
				cfR1 = 0;
				cfS0 = 0;
				cfS1 = 0;

				if (j != 0){ 				// singular point NOT on the lower geometric node
					cfR0 = -  lD0*gsl_sf_log(lD0);
					cfS0 = -2*lD0;
				}

				if (j != nlocl-1){	// singular point NOT on the upper geometric node
					cfR1 = -  lD1*gsl_sf_log(lD1);
					cfS1 = -2*lD1;
				}

				for (k = 0; k < nquad; k++){ // regular quadrature
					lR0[k] = lsing - 0.5*lD0*(1 + zquad[k]);
					lR1[k] = lsing + 0.5*lD1*(1 + zquad[k]);

					zR0[k] = (lR0[k] - lM)/lD;
					zR1[k] = (lR1[k] - lM)/lD;
				}

				for (k = 0; k < nsing; k++){ // singular quadrature
					lS0[k] = lsing - lD0*zqsng[k];
					lS1[k] = lsing + lD1*zqsng[k];

					zS0[k] = (lS0[k] - lM)/lD;
					zS1[k] = (lS1[k] - lM)/lD;
				}

				lagrange(nlocl-1, nquad, zlocl, zR0, LR0);
				lagrange(nlocl-1, nquad, zlocl, zR1, LR1);
				lagrange(nlocl-1, nsing, zlocl, zS0, LS0);
				lagrange(nlocl-1, nsing, zlocl, zS1, LS1);
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

					/* regularize diagonal components of the free-space
					 * Green's function if near the singular point */
					if (ising == 1){
						logl  = gsl_sf_log(fabs(l - lsing));
						
						MRxx += 2*logl;
						MRrr += 2*logl;

					}

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

					// add free-space and complementary Green's functions
					Mxx = MRxx + MCxx;
					Mxr = MRxr + MCxr;
					Mrx = MRrx + MCrx;
					Mrr = MRrr + MCrr;

					// increment the integrals
					Axx += w*Mxx*Lq;
					Axr += w*Mxr*Lq;
					Arx += w*Mrx*Lq;
					Arr += w*Mrr*Lq;

				} // end of quadrature points

				/* add back singular part of the integral
				 * using logarithmic quadrature */
				if (ising == 1){
					// prepare for quadrature
					IR0 = 0; // regular integral
					IR1 = 0; // regular integral
					IS0 = 0; // singular integral
					IS1 = 0; // singular integral

					// evaluate regular quadrature
					for (k = 0; k < nquad; k++){ // loop over quadrature points
						/*--------- first integral, IR0 ----------*/
						// get Lagrange interpolating polynomial
						Lq = LR0[nquad*j + k];
						
						// get metric coefficient
						l  = lR0[k];
						dl = l  - l0;

						dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
						drdl = (3.*ar*dl + 2.*br)*dl + cr;
						dsdl = sqrt(dxdl*dxdl + drdl*drdl);

						h = dsdl;

						// get weight for kth node
						w = h*wquad[k];

						// increment integral
						IR0 += w*Lq;
						
						/*-------- second integral, IR1 ----------*/
						// get Lagrange interpolating polynomial
						Lq = LR1[nquad*j + k];

						// get metric coefficient
						l  = lR1[k];
						dl = l  - l0;
		
						dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
						drdl = (3.*ar*dl + 2.*br)*dl + cr;
						dsdl = sqrt(dxdl*dxdl + drdl*drdl);

						h = dsdl;

						// get weight for kth node
						w = h*wquad[k];

						// increment integral
						IR1 += w*Lq;
					} // end of regular quadrature points

					// evaluate singular quadrature
					for (k = 0; k < nsing; k++){ // loop over quadrature points
						// get singular kernel
						z    = zqsng[k];
						logz = gsl_sf_log(z);
						logz = -1;


						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1
						// NOTE: GET RID OF LOGZ AND REPLACE WITH -1


						/*--------- third integral, IS0 ----------*/
						// get Lagrange interpolating polynomial
						Lq = LS0[nsing*j + k];
						
						// get metric coefficient
						l = lS0[k];
						dl = l  - l0;

						dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
						drdl = (3.*ar*dl + 2.*br)*dl + cr;
						dsdl = sqrt(dxdl*dxdl + drdl*drdl);

						h = dsdl;

						// get weight for kth node
						w = h*wqsng[k];

						// increment integral
						IS0 += w*logz*Lq;

						
						/*-------- fourth integral, IS1 ----------*/
						// get Lagrange interpolating polynomial
						Lq = LS1[nsing*j + k];

						// get metric coefficient
						l = lS1[k];
						dl = l  - l0;
		
						dxdl = (3.*ax*dl + 2.*bx)*dl + cx;
						drdl = (3.*ar*dl + 2.*br)*dl + cr;
						dsdl = sqrt(dxdl*dxdl + drdl*drdl);

						h = dsdl;

						// get weight for kth node
						w = h*wqsng[k];

						// increment integral
						IS1 += w*logz*Lq;
					} // end of singular quadrature points
					
					IR0 *= cfR0;
					IR1 *= cfR1;
					IS0 *= cfS0;
					IS1 *= cfS1;

					Ising = IR0 + IR1 + IS0 + IS1;

					Axx += Ising;
					Arr += Ising;

				} // end of singular integral

				// add 2x2 block of A
				A[2*nglob*(2*m)   + (2*n  )] += Axx;
				A[2*nglob*(2*m)   + (2*n+1)] += Axr;
				A[2*nglob*(2*m+1) + (2*n  )] += Arx;
				A[2*nglob*(2*m+1) + (2*n+1)] += Arr;
					
			} // end of local element nodes (index j)
		} // end of boundary elements (index i)
	} // end of global element nodes (index m)

	
	/*---------------------------------------------*/
	/*------- CALCULATE SINGLE LAYER DENSITY ------*/
	/*---------------------------------------------*/

	Surface.calcForce(Dfx, Dfr,
	                   ks,  kp,
										 tx,  tr,
										 nx,  nr);
	for (i = 0; i < nglob; i++){ // loop over global element nodes
		Df[2*i  ] = Dfx[i];
		Df[2*i+1] = Dfr[i];
	}
	

	/*---------------------------------------------*/
	/*------ CALCULATE SINGLE LAYER POTENTIAL -----*/
	/*---------------------------------------------*/

	cf = -1./(8.*M_PI);

	for (i = 0; i < nglob; i++){		// loop over collocation points
		// initialize
		vx = 0.;
		vr = 0.;

		for (j = 0; j < nglob; j++){	// loop over collocation points
			vx += A[2*nglob*(2*i  ) + (2*j  )]*Df[2*j  ];
			vx += A[2*nglob*(2*i  ) + (2*j+1)]*Df[2*j+1];
			vr += A[2*nglob*(2*i+1) + (2*j  )]*Df[2*j  ];
			vr += A[2*nglob*(2*i+1) + (2*j+1)]*Df[2*j+1];
		}

		// store velocity
		v[2*i  ] = cf*vx;
		v[2*i+1] = cf*vr;
	}




//	/*------ FOR DEBUGGING -------*/
//
//	printf("\nzlocl = ");
//	for (i = 0; i < nlocl; i++){
//    if (zlocl[i] >= 0)
//      printf(" ");
//		printf("%.4f ", zlocl[i]);
//	}
//	printf("\n");
//
//	printf("\nzquad = ");
//	for (i = 0; i < nquad; i++){
//    if (zquad[i] >= 0)
//      printf(" ");
//		printf("%.4f ", zquad[i]);
//	}
//	printf("\n");
//
//	printf("\nwquad = ");
//	for (i = 0; i < nquad; i++){
//    if (wquad[i] >= 0)
//      printf(" ");
//		printf("%.4f ", wquad[i]);
//	}
//	printf("\n");
//	printf("\n");
//
//	printf("L = \n");
//  for (i = 0; i < nlocl; i++){    	// rows = polynomial
//    for (j = 0; j < nquad; j++){    // columns = grid point
//      if (L[i*nquad+j] >= 0)
//        printf(" ");
//      printf("%.4f ", L[i*nquad + j]);
//    }   
//    printf("\n");
//  }
//	printf("\n");
//
//  printf("A = \n");
//  for (i = 0; i < 2*nglob; i++){
//    for (j = 0; j < 2*nglob; j++){
//      if (A[2*nglob*i + j] >= 0 && A[2*nglob*i + j] < 10) 
//        printf(" ");
//      printf("%.2f ", A[2*nglob*i + j]);
//    }   
//    printf("\n");
//  }
//  printf("\n");
//
//  printf("Df = \n");
//  for (i = 0; i < nglob; i++){
//    for (j = 0; j < 2; j++){
//      if (Df[2*i+j] >= 0 && Df[2*i+j] < 10) 
//        printf(" ");
//      printf("%.2f", Df[2*i+j]);
//      printf(" ");
//    }   
//    printf("\n");
//  }
//  printf("\n");
//
//  printf("v = \n");
//  for (i = 0; i < nglob; i++){
//    for (j = 0; j < 2; j++){
//      if (v[2*i+j] >= 0 && v[2*i+j] < 10) 
//        printf(" ");
//      printf("%.2f", v[2*i+j]);
//      printf(" ");
//    }   
//    printf("\n");
//  }
//
//	/*--------------------------*/






	// release memory
	free(A    );

  free(xc   );
  free(rc   );
	free(Dfx  );
	free(Dfr  );
	free(ks   );
	free(kp   );
	free(tx   );
	free(tr   );
	free(nx   );
	free(nr   );

  free(zlocl);
  free(zquad);
  free(wquad);
	free(L    );

  free(zqsng);
  free(wqsng);
	free(LR0  );
	free(LR1  );
	free(LS0  );
	free(LS1  );
	free(zR0  );
	free(zR1  );
	free(zS0  );
	free(zS1  );
	free(lR0  );
	free(lR1  );
	free(lS0  );
	free(lS1  );
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
	
				




///* Determine native element nodes:
// *  M = 0	- uniform elements
// *  M = 1 - linear basis
// *  M = 2 - quadratic basis
// */
//void be_native(int M, int N, double *XG, double *YG, int n1, int n2){
//	// check value of M
//	if (M != 0 && M != 1 && M != 2){
//		printf("Error: M must equal 0, 1, or 2");
//		return;
//	}
//	
//	// uniform elements
//	if (M == 0){
//		
//	}
//	
//	// linear basis in Lagrange polynomials
//	if (M == 1){
//
//	}
//
//	// quadratic basis in Lagrange polynomials
//	if (M == 2){
//
//	}
//}


#endif
