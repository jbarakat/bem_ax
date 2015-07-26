/* GEOMETRY
 *  Shape class that prescribes the boundary elements and geometrical
 *  parameters of an axisymmetric contour with N+1 nodal points (x,r) 
 *  and interpolated using cubic splines.
 *
 * REFERENCES
 *  Item #1
 *  Item #2
 *  
 * PARAMETERS
 *  x,r   [input]			nodal coordinates
 *  N	    [input]			number of boundary elements
 *  s     [output]		meridional arc length
 *  A     [output]		total area
 *  V     [output]		total volume
 *  cs,cp [output]		principal curvatures
 *  t     [output]		meridional tangent vector
 *  n     [output] 		normal vector
 */

#ifndef GEOM_H
#define GEOM_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "interp.h"

class geom {
private:
	int nnode, nelem;
	double *nodex, *noder;
	double *arcl, area, vlme;
	double *curvs, *curvp;
	double *tangx, *tangr;
	double *nrmlx, *nrmlr;

public:
	/* PROTOTYPES */
	
	/* IMPLEMENTATIONS */

	// Constructors
	geom(){
	}

	geom(int N, double *x, double *r){
		// assign nummber of nodes and boundary elements
		nnode = N + 1;
		nelem = N;

		// allocate memory for pointer arrays
		arcl  = (double*) malloc(nnode * sizeof(double));
		nodex = (double*) malloc(nnode * sizeof(double));
		noder = (double*) malloc(nnode * sizeof(double));
		curvs = (double*) malloc(nnode * sizeof(double));
		curvp = (double*) malloc(nnode * sizeof(double));
		tangx = (double*) malloc(nnode * sizeof(double));
		tangr = (double*) malloc(nnode * sizeof(double));
		nrmlx = (double*) malloc(nnode * sizeof(double));
		nrmlr = (double*) malloc(nnode * sizeof(double));
		
		// calculate geometric parameters
		calcParams(N, x, r,
		           arcl, area, vlme, 
							 curvs, curvp, 
		           tangx, tangr,
							 nrmlx, nrmlr);

	}

	// Destructor
//	~geom(){
//	}
	
	// Set functions

	// Get functions
	void getArcl(double *s){
		int i;
		
		if (s == NULL){
			printf("Error: no memory allocated for s.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			s[i] = arcl[i];
		}
	}

	double getArea(){
		double A;
		A = area;
		return(A);
	}
	
	void getArea(double &A){
		A = area;
	}
	
	double getVlme(){
		double V;
		V = vlme;
		return(V);
	}

	void getVlme(double &V){
		V = vlme;
	}
	
	void getCurv(double *cs, double *cp){
		int i;
		
		if (cs == NULL || cp == NULL){
			printf("Error: no memory allocated for cs, cp.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			cs[i] = curvs[i];
			cp[i] = curvp[i];
		}
	}

	void getTang(double *tx, double *tr){
		int i;
		
		if (tx == NULL || tr == NULL){
			printf("Error: no memory allocated for tx, tr.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			tx[i] = tangx[i];
			tr[i] = tangr[i];
		}
	}

	void getNrml(double *nx, double *nr){
		int i;
		
		if (nx == NULL || nr == NULL){
			printf("Error: no memory allocated for nx, nr.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			nx[i] = nrmlx[i];
			nr[i] = nrmlr[i];
		}
	}

	/* Function to calculate geometric parameters for a given set of N+1
	 * nodal coordinates (x,r) */
	void calcParams(int N, double *x, double *r, 				/* <--- inputs  */
									double *s, double &A, double &V,		/* <--- outputs */
									double *cs, double *cp,             /* <------|     */
									double *tx, double *tr,							/* <------|     */
									double *nx, double *nr){						/* <------|     */
		// declare variables
		int i,j, n;
		int na, nt;
		double *l, lj;
		double *ax, *bx, *cx;
		double *ar, *br, *cr;
		double axi, bxi, cxi;
		double ari, bri, cri;
		double xi, ri, xj, rj;
		double dx, dr, ds, dA, dV, dl;
		double dxdl, drdl, dsdl, dAdl, dVdl;
		double d2xdl2, d2rdl2;
		double ssum, Asum, Vsum;
		const int MAXIT = 20;
		double stol, Atol, Vtol;
		double fc;

		// allocate memory
		l  = (double*) malloc((N+1) * sizeof(double));
		ax = (double*) malloc( N    * sizeof(double));
		bx = (double*) malloc((N+1) * sizeof(double));
		cx = (double*) malloc( N    * sizeof(double));
		ar = (double*) malloc( N    * sizeof(double));
		br = (double*) malloc((N+1) * sizeof(double));
		cr = (double*) malloc( N    * sizeof(double));	

		// calculate polygonal arc length
		l[0] = 0;
		for (i = 1; i < N+1; i++){
			dx = x[i] - x[i-1];
			dr = r[i] - r[i-1];
			dl = sqrt(dx*dx + dr*dr);
			l[i] = l[i-1] + dl;
		}

		// set integration tolerances
		dl = l[N] - l[0];
		stol = 0.0000001*dl;
		Atol = 0.0000001*dl*dl;
		Vtol = 0.0000001*dl*dl*dl;
		
		// calculate cubic spline coefficients for interpolation
		spline(N, l, x, 0.,  0., ax, bx, cx);
		spline(N, l, r, 1., -1., ar, br, cr);

		/* calculate meridional arc length, area, and volume 
		 * using the extended trapezoidal rule */
		s[0] = 0.;
		A    = 0.;
		V    = 0.;
		for (i = 0; i < N; i++){
			// update current step
			 xi =  x[i];
			 ri =  r[i];
			axi = ax[i];
			bxi = bx[i];
			cxi = cx[i];
			ari = ar[i];
			bri = br[i];
			cri = cr[i];

			// initialize;
			ssum = 0.;
			Asum = 0.;
			Vsum = 0.;
			dl   = l[i+1] - l[i];
			
			// evaluate integrand at l[i], l[i+1]
			for (j = i; j < i+2; j++){
				rj    =  r[j];
				dxdl  =  cxi;
				drdl  =  cri;
				dsdl  =  sqrt(dxdl*dxdl + drdl*drdl);
				dAdl  =  rj*dsdl;
				dVdl  = -rj*rj*dxdl;
				
				ssum += 0.5*dsdl;
				Asum += 0.5*dAdl;
				Vsum += 0.5*dVdl;
			}

			na = 2; // number of points added
			nt = 2; // total number of points

			for (n = 1; n < MAXIT; n++){
				// add 2^(n-1) additional points
				if (n == 1)
					na  = 1;
				else
					na *= 2;

				nt += na;

				// store integrands before adding new points
				ds = ssum*dl;
				dA = Asum*dl;
				dV = Vsum*dl;

				// refine grid spacing
				dl /= 2.;

				// interpolate integrand at additional points
				for (j = 0; j < nt; j++){
					if (j % 2 != 0){ /* avoid double counting
					                  * previous grid points */
						lj    = j*dl;

						rj    =  ((ari*lj + bri)*lj + cri)*lj + ri;
						dxdl  =  (3.*axi*lj + 2.*bxi)*lj + cxi;
						drdl  =  (3.*ari*lj + 2.*bri)*lj + cri;
						dsdl  =  sqrt(dxdl*dxdl + drdl*drdl);
						dAdl  =  rj*dsdl;
						dVdl  = -rj*rj*dxdl;
						
						ssum += dsdl;
						Asum += dAdl;
						Vsum += dVdl;
					}
				}

				// calculate change in integrand
				ds -= ssum*dl; ds = fabs(ds);
				dA -= Asum*dl; dA = fabs(dA);
				dV -= Vsum*dl; dV = fabs(dV);

				// break loop when integrals converge
				if (ds < stol && dA < Atol && dV < Vtol)
					break;
			}
		  
			// diagnose level of refinement [uncomment when required]
		 // printf("%d stages of refinement and %d total points\n", n, nt);
			
			// increment the integrals
			ds = ssum*dl;
			dA = Asum*dl;
			dV = Vsum*dl;

			s[i+1] = s[i] + ds;
			A     += dA;
			V     += dV;
		}

		A *= 2*M_PI;
		V *= M_PI;
		
		/* calculate principal curvatures, meridional tangent vector,
		 * and outward normal vector at the nodal points */
		cs[0]     =  0.;
		cp[0]     =  0.;
		
		tx[0]     =  0.;
		tr[0]     =  1.;
		
		nx[0]     =  1.;
		nr[0]     =  0.;

		for (i = 1; i < N; i++){
			ri      =  r[i];
			bxi     =  bx[i];
			bri     =  br[i];
			cxi     =  cx[i];
			cri     =  cr[i];
			dxdl    =  cxi;
			d2xdl2  =  2*bxi;
			drdl    =  cri;
			d2rdl2  =  2*bri;
			dsdl    =  sqrt(cxi*cxi + cri*cri);

			cs[i]   = -(dxdl*d2rdl2 - d2xdl2*drdl)/pow(dsdl, 3);
			cp[i]   =  dxdl/(ri*dsdl);
			
			tx[i]   =  dxdl/dsdl;
			tr[i]   =  drdl/dsdl;
			
			nx[i]   =  tr[i];
			nr[i]   = -tx[i];
		}
	
		fc        = (s[0]   - s[2])/(s[1] - s[2]);
		cs[0]     = cs[1]   + fc*(cs[1]   - cs[2]);
		cp[0]     = cp[1]   + fc*(cp[1]   - cp[2]);

		fc        = (s[N] - s[N-2])/(s[N-1] - s[N-2]);
		cs[N]     = cs[N-1] + fc*(cs[N-1] - cs[N-2]);
		cp[N]     = cp[N-1] + fc*(cp[N-1] - cp[N-2]);
		
		tx[N]     =  0.;
		tr[N]     = -1.;
		
		nx[N]     = -1.;
		nr[N]     =  0.;
		
		// release memory
		free(l);
		free(ax);
		free(bx);
		free(cx);
		free(ar);
		free(br);
		free(cr);

	}
	
};

#endif
