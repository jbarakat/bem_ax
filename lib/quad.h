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
#include "stokes.h"
#include "gauleg.h"

/* PROTOTYPES */

/* IMPLEMENTATIONS */

/* Evaluate single layer potential over the boundary surface 
 *  IGF = 0  free-space axisymmetric Green's function
 *  IGF = 1  tube-bounded Green's function
 */
void singleLayer(const int IGF, int ngl, stokes Stokes){
	if (IGF != 0 && IGF != 1){
		printf("Error: IGF can only take values of 0 or 1");
		return;
	}
	
	// declare variables
	int i, j, k, nnode, nelem;
	double A, V;
	
	double  ax,  bx, cx;
	double  ar,  br, cr;
	
	double  x ,  r ;
	double  xp,  rp;
	double  Dl,  l , dl;
	double  w;

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
	
	// get number of boundary elements
	nelem = Stokes.getNElem();
	nnode = nelem + 1;

	// allocate memory
  zgl   = (double*) malloc( ngl * sizeof(double));
  wgl   = (double*) malloc( ngl * sizeof(double));

	// get area and volume
	A     = Stokes.getArea();
	V     = Stokes.getVlme();

	// get Gauss-Legendre abscissas and weights
	gauleg(ngl, zgl, wgl);
	
	/* loop over boundary nodes and evaluate integrals
	 * at field points (xp,rp) */
	for (i = 0; i < nnode; i++){
		


		/* loop over boundary elements and carry out
		 * quadrature over source points (x,r) */
		for (j = 0; j < nelem; j++){
			// update parameters for the jth boundary element
			Stokes.getNode(j,    x0,  r0);
			Stokes.getNode(j+1,  x1,  r1);
			Stokes.getPoly(j,    l0);
			Stokes.getPoly(j+1,  l1);
	
			Stokes.getSpln(j,   ax , bx , cx ,
			                    ar , br , cr );
	
			Stokes.getNrml(j,   nx0, nr0);
			Stokes.getNrml(j+1, nx1, nr1);
			
			// calculate length of polygonal line segment
			Dl = l1 - l0;
	
			// loop over Gauss-Legendre grid points
			for (k = 0; k < ngl; k++){
				// get weight for kth node
				w  = wgl[k];
				
				// map grid point onto polygonal arc length
				l  = l0 + 0.5*Dl*(zgl[k] + 1);
				dl = l  - l0;
	
				// interpolate x and r (source point)
				x  = ((ax*dl + bx)*dl + cx)*dl + x0;
				r  = ((ar*dl + br)*dl + cr)*dl + r0;

				// calculate Green's functions
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
	
				
			}
	
			// NEED WAY TO GET BOUNDARY TRACTION AND BOUNDARY VELOCITY 
	
			// map l domain onto [-1, 1] domain
			//
			//  xi = -1 + (2/(l[i+1] - l[i]))*(l - l[i])
			//
			//  or map xi onto l domain
			//
			//  l = l[i] + 0.5*(l[i+1] - l[i])*(1 + xi) 
			
		} // end of boundary elements (integration over source points)
	} // end of boundary nodes (evaluation at field points)
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
