/* QUADRATURE
 *  Integrate the Green's function over a spline element using
 *  Gauss-Legendre quadrature.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [Ch. 6]
 *  
 * PARAMETERS
 *  n      [input]    nodal index
 *  NGL    [input]		number of quadrature points
 *  G      [input]		geometry of the boundary
 */

#ifndef BEDISC_H
#define BEDISC_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>
#include "geom.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Evaluate single layer potential over the boundary surface */
void singleLayer(int NGL, geom Geom){
	// declare variables
	int i, j, nnode, nelem;
	
	double *x,  *r;
	double *ax, *bx, *cx;
  double *ar, *br, *cr;
  double *s, A, V;
  double *ks, *kp;
  double *tx, *tr;
  double *nx, *nr;

	double  x0,  r0;
	double ax0, bx0, cx0;
	double ar0, br0, cr0;
	double  s0;
	double ks0, kp0;
	double tx0, tr0;
	double nx0, nr0;
	
	double  x1,  r1;
	double  s1;
	double ks1, kp1;
	double tx1, tr1;
	double nx1, nr1;
	
	// get number of boundary elements
	nelem = Geom.getNElem();
	nnode = nelem + 1;

	// allocate memory
  x     = (double*) malloc( nnode    * sizeof(double));
  r     = (double*) malloc( nnode    * sizeof(double));
  ax    = (double*) malloc((nnode-1) * sizeof(double));
  bx    = (double*) malloc( nnode    * sizeof(double));
  cx    = (double*) malloc((nnode-1) * sizeof(double));
  ar    = (double*) malloc((nnode-1) * sizeof(double));
  br    = (double*) malloc( nnode    * sizeof(double));
  cr    = (double*) malloc((nnode-1) * sizeof(double));
  s     = (double*) malloc( nnode    * sizeof(double));
  ks    = (double*) malloc( nnode    * sizeof(double));
  kp    = (double*) malloc( nnode    * sizeof(double));
  tx    = (double*) malloc( nnode    * sizeof(double));
  tr    = (double*) malloc( nnode    * sizeof(double));
  nx    = (double*) malloc( nnode    * sizeof(double));
  nr    = (double*) malloc( nnode    * sizeof(double));

	// get geometric parameters	
	Geom.getAll(nelem, x , r ,
              ax,    bx, cx,
              ar,    br, cr,
              s,     A ,  V, 
              ks,    kp,
              tx,    tr,
              nx,    nr);
	
	//

	for (i = 0; i < nelem; i++){
		// update parameters for the ith boundary element
		 x0 = x [i];
		 r0 = r [i];
		ax0 = ax[i];
		bx0 = bx[i];
		cx0 = cx[i];
		ar0 = ar[i];
		br0 = br[i];
		cr0 = cr[i];
		 s0 = s [i];
		ks0 = ks[i];
		kp0 = kp[i];
		tx0 = tx[i];
		tr0 = tr[i];
		nx0 = nx[i];
		nr0 = nr[i];
		
		 x1 = x [i+1];
		 r1 = r [i+1];
	 	 s1 = s [i+1];
		ks1 = ks[i+1];
		kp1 = kp[i+1];
		tx1 = tx[i+1];
		tr1 = tr[i+1];
		nx1 = nx[i+1];
		nr1 = nr[i+1];
		
		// 
	}
}

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
