/* INTERPOLATION
 *  Generate interpolants using cubic splines or Lagrange polynomials.
 *
 * REFERENCES
 *  Moin, Cambridge University Press (2010) (Ch. 1)
 *  Pozrikidis, Chapman & Hall/CRC (2002) (Ch. 3)
 *  
 * PARAMETERS
 *  x,y   [input]		set of N+1 grid points
 *  xi,yi [output]  interpolated grid point
 *  N     [input]   number of segments
 *  a     [output]	spline coefficient of 3rd derivative
 *  b     [output]	spline coefficient of 2nd derivative
 *  c     [output]	spline coefficient of 1st derivative
 */

#ifndef INTERP_H
#define INTERP_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <lapacke.h>
#include <gsl/gsl_sf_trig.h>

/* PROTOTYPES */
void spline(const int, double*, double*, double, double,
            double, double&);
void spline(const int, double*, double*, double, double,
            double*, double*, double*);
void lagrange(const int, double*, double*, double, double&);

/* IMPLEMENTATIONS */
/***********************************************************************/

/* Interpolate to point (xi,yi) based on a set of N+1 nodes (x,y) using
 * cubic splines.
 */
void spline(const int N, double *x, double *y, double slope1, double slope2,
            double xi, double &yi){
	// declare variables
  int i, j, n;
	double *a, *b, *c;
	double dx;

	// allocate memory
	a = (double*) malloc( N    * sizeof(double));
	b = (double*) malloc((N+1) * sizeof(double));
	c = (double*) malloc( N    * sizeof(double));

	// compute cubic spline coefficients
	spline(N, x, y, slope1, slope2, a, b, c);

	/* find where xi falls in the domain 
	 *(assume x increases monotonically along the boundary) */
	n = -1;
	for (i = 0; i < N; i++){
		if (xi > x[i]){
			n = i;
			break;
		}
	}
	
	if (n < 0){
		printf("Error: xi not in x.");
		return;
	}

	// compute yi
	dx = xi - x[n];
	yi = x[n] + ((a[n]*dx + b[n])*dx + c[n])*dx;
	
	// release memory
	free(a);
	free(b);
	free(c);
	
}

/* Generate coefficients a, b, c for cubic spline interpolation of
 * N+1 points (x,y) by solving a tridiagonal system of equations
 * with clamped boundary conditions. Note that
 *  a,c  contain N elements, and
 *  b    contains N+1 elements.
 */
void spline(const int N, double *x, double *y, double slope1, double slope2,
            double *a, double *b, double *c){
	// declare variables
	int i, j, info;
	double *h, *rhs, *sln;
	double *dl, *d, *du;

	// allocate memory
	h   = (double*) malloc( N    * sizeof(double));
	dl  = (double*) malloc((N-2) * sizeof(double));
	d   = (double*) malloc((N-1) * sizeof(double));
	du  = (double*) malloc((N-2) * sizeof(double));
	rhs = (double*) malloc((N-1) * sizeof(double));
	sln = (double*) malloc((N-1) * sizeof(double));
	
	// compute intervals h
	for (i = 0; i < N; i++){
		h[i] = x[i+1] - x[i];
	}

	// generate tridiagonal matrix for N-1 equations
	d [0] = 2.*(h[0] + h[1]) - 0.5*h[0];
	dl[0] = h[1];
	du[0] = h[1];

	for (i = 1; i < N-2; i++){
		d [i] = 2.*(h[i] + h[i+1]);
		dl[i] = h[i+1];
		du[i] = h[i+1];
	}

	d[N-2] = 2.*(h[N-2] + h[N-1]) - 0.5*h[N-1];

	for (i = 0; i < N-1; i++){
		rhs[i] = 3.*((y[i+2] - y[i+1])/h[i+1] - (y[i+1] - y[i])/h[i]);
	}

	rhs[0]   -= 1.5*((y[1] - y[0])  /h[0]   - slope1);
	rhs[N-2] += 1.5*((y[N] - y[N-1])/h[N-1] - slope2);

	// solve the tridiagonal system
	info = LAPACKE_dgtsv(LAPACK_COL_MAJOR, N-1, 1, dl, d, du, rhs, N-1);
	
	for (i = 0; i < N-1; i++){
		sln[i] = rhs[i];
	}

	// compute the coefficients a, b, c
	for (i = 0; i < N-1; i++){
		b[i+1] = sln[i];
	}

	b[0] = -0.5*b[1]   + 1.5*((y[1] - y[0]  )/h[0]   - slope1)/h[0];
	b[N] = -0.5*b[N-1] - 1.5*((y[N] - y[N-1])/h[N-1] - slope2)/h[N-1];
	
	for (i = 0; i < N; i++){
		a[i] = (b[i+1] - b[i])/(3.*h[i]);
		c[i] = (y[i+1] - y[i])/h[i] - h[i]*(b[i+1] + 2.*b[i])/3.;
	}
	
	// release memory
	free(h);
	free(dl);
	free(d);
	free(du);
	free(rhs);
	free(sln);
}


/* Interpolate to point (xi,yi) based on a set of N+1 nodes (x,y) using
 * Lagrange polynomials.
 */
void lagrange(const int N, double *x, double *y, double xi, double &yi){
  // declare variables
  int i, j;
  double *rho, psi, L, P, dx; 
  double idx;
  const double TOL = 1e-12;
  
  // check if x contains xi
  idx = 1.; 
  for (i = 0; i < N+1; i++){
    dx = fabs(xi - x[i]);
    idx = fmin(idx, dx);
    if (idx < TOL){
      yi = y[i];
      return;
    }   
  }

  // allocate memory
  rho = (double*) malloc((N+1) * sizeof(double));

  // construction step
  for (i = 0; i < N+1; i++){
    rho[i] = 1.; 
    for (j = 0; j < N+1; j++){
      dx = x[i] - x[j];
      if (j != i)
        rho[i] *= dx; 
    }   
  }

  // evaluation step
  psi = 1.; 
  L = 0;
  P = 0;

  for (i = 0; i < N+1; i++){
    dx = xi - x[i];
    psi *= dx;
    L = 1/(dx*rho[i]);
    P += y[i]*L;
  }

	P *= psi;
  yi = P;

  // free memory
  free(rho);
}



/***********************************************************************/

#endif
