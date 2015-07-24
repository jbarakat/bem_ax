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
 *  x,r [input]			nodal coordinates
 *  N		[input]			number of boundary elements
 *  s   [output]		meridional arc length
 *  V   [output]		volume
 *  A   [output]		area
 *  g   [output]		metric tensor
 *  h   [output]		curvature tensor
 *  H   [output]		mean curvature
 *  K   [output]		Gaussian curvature
 *  t   [output]		meridional tangent vector
 *  n   [output] 		normal vector
 */

#ifndef GEOM_H
#define GEOM_H

/* HEADER FILES */
#include <math.h>
#include "interp.h"

class geom {
public:
	/* PROTOTYPES */
	
	
	/* IMPLEMENTATIONS */
//	void geom(const int N, double *x, double *r,
//						double *s, double &V, double &A,
//						double *gss, double *gpp,
//						double *hss, double *hpp,
//	          double *tx, double *tr,
//	          double *nx, double *nr){
	void calcParams(const int N, double *x, double *r,
						double *s, double &V, double &A){
		// declare variables
		int i, j, n;
		int na, nt;
		double step;
		double *l, l0, l1, lj;
		double *ax, *bx, *cx;
		double *ar, *br, *cr;
		double axi, bxi, cxi;
		double ari, bri, cri;
		double dx, dr, ds, dl;
		double dxdl, drdl, dsdl, dAdl, dVdl;
		const int MAXIT = 40;
		const double TOL = 0.00001;

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

		// calculate cubic spline coefficients
		spline(N, l, x, 0.,  0., ax, bx, cx);
		spline(N, l, r, 1., -1., ar, br, cr);

		/* calculate the meridional arc length, area, and volume 
		 * using the extended trapezoidal rule */
		s[0] = 0.;
		for (i = 0; i < N; i++){
			// update current step
			l0 = l[i];
			l1 = l[i+1];
			axi = ax[i];
			bxi = bx[i];
			cxi = cx[i];
			ari = ar[i];
			bri = br[i];
			cri = cr[i];

			// initialize;
			dsdl = 0.;
			dl  = l1 - l0;
			
			// evaluate integrand at l[i], l[i+1]
			for (j = i; j < i+2; j++){
				lj = l[i];
				dxdl  = (3.*axi*lj + 2.*bxi)*lj + cxi;
				drdl  = (3.*ari*lj + 2.*bri)*lj + cri;
				step   = sqrt(dxdl*dxdl + drdl*drdl);
				dsdl += 0.5*step;
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
				ds = dsdl*dl;

				// refine grid spacing
				dl /= 2.;

				// evaluate integrand at additional points
				for (j = 0; j < nt; j++){
					if (j % 2 != 0){ /* avoid double counting
					                  * previous grid points */
						lj = j*dl;
						dxdl  = (3.*axi*lj + 2.*bxi)*lj + cxi;
						drdl  = (3.*ari*lj + 2.*bri)*lj + cri;
						step   = sqrt(dxdl*dxdl + drdl*drdl);
						dsdl += step;
					}
				}
				
				// calculate change in integrand
				ds -= dsdl*dl; ds = fabs(ds);
				
				// break loop when integral converges
				if (ds < TOL)
					break;

			}
		  // diagnose level of refinement [uncomment when required]
//		  printf("%d stages of refinement and %d total points\n", n, nt);

			ds = dsdl*dl;
			s[i+1] = s[i] + ds;

		}

		A = 0.;
		V = 0.;

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
