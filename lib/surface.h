/* SURFACE MECHANICS
 *  Boundary class that inherits the members of STOKES.H with type = 1
 *  (fluid-fluid interface) and, in addition, prescribes a surface
 *  constitutive model for the in-plane tensions (shear and dilatation),
 *  transverse shear tensions, and in-plane moments (bending and twisting).
 *
 * REFERENCES
 *  Mollman, John Wiley & Sons (1982)
 *  
 * PARAMETERS
 *  tau   [output]		in-plane tensions
 *  q     [output]		transverse shear tension
 *  m			[output]		moments
 */

#ifndef SURFACE_H
#define SURFACE_H

/* HEADER FILES */
#include "stokes.h"

class surface: public stokes {
private:
  // indicator for constitive model
  int model;
  /*  = 0  (drop w/constant surfactant)
   *  = 1  (drop w/varying surface tension)
   *  = 2  (vesicle w/fixed surface area - Helfrich model)
   *  = 3  (vesicle w/variable surface area)
   *  = 4  (red blood cell - Skalak model) */
  
  // in-plane tensions
  double *tenss, *tensp;

  // transverse shear tensions
  double *tensn;

  // moments
  double *mmnts, *mmntp;
 
public:
  /* PROTOTYPES */
  
  
  /* IMPLEMENTATIONS */
  surface() : stokes() {
  }

  surface(int id, int N, int M,
          double gamm,
          double ES, double ED, double EB, double ET,
          double *x, double *r) : stokes(1, N, M, x, r) {
    // declare variables
    int i, j;

    // set constitutive model
    model = id;

    /* allocate memory for pointer arrays
     * and initialize to zero */
    tenss = (double*) calloc(nglob, sizeof(double));
    tensp = (double*) calloc(nglob, sizeof(double));
    tensn = (double*) calloc(nglob, sizeof(double));
    mmnts = (double*) calloc(nglob, sizeof(double));
    mmntp = (double*) calloc(nglob, sizeof(double));

    if (id == 0){ // drop w/constant surface tension
      for (i = 0; i < nglob; i++){
        tenss[i] = gamm;
        tensp[i] = gamm;
        tensn[i] = 0;
        mmnts[i] = 0;
        mmnts[i] = 0;
      }
    }
    
    if (id == 1){ // drop w/varying surface tension
      for (i = 0; i < nglob; i++){
         
        // NEED TO FINISH THIS

      }
    }
    
    if (id == 2){ // vesicle w/constant surface area
      for (i = 0; i < nglob; i++){
        
        // NEED TO FINISH THIS

      }
    }

    // ADD IF STATMENTS FOR id == 3 (vesicle w/non-constant area)
    // and id == 4 (red blood cell)

  }
  
  /*- Set functions ---*/

  /*- Get functions ---*/
  
  // get tension components at the (iglob)th global basis node
  void getTens(int iglob, double &taus, double &taup, double &q){
    if (iglob >= nglob){
      printf("Error: index out of bounds.\n");
      return;
    }

    taus = tenss[iglob];
    taup = tensp[iglob];
    q    = tensn[iglob];
  }

  // get tension components at all global basis nodes
  void getTens(double *taus, double *taup, double *q){
    int iglob;

    if (taus == NULL || taup == NULL || q == NULL){
      printf("Error: no memory allocated for taus, taup, q.\n");
    }

    for (iglob = 0; iglob < nglob; iglob++){
      taus[iglob] = tenss[iglob];
      taup[iglob] = tensp[iglob];
      q   [iglob] = tensn[iglob];
    }
  }

  /* Function to calculate the force resultants normal and tangential
   * (in the meridional direction) to a free surface with constitutive
   * tensions and moments. The force resultants are evaluated at the
	 * basis nodes (collocation points) of the surface. */
  void calcForce(surface Surface,
                 double *fn, double *fs){
		// declare variables
		int i, j, k;
		double *zlocl, *L, *dLdx;
		double  cf;

		// allocate memory
		zlocl = (double*) malloc(nlocl * sizeof(double));

		/* calculate Gauss-Lobatto points on the
		 * interval [-1,1] */
		// REDUNDANT CALCULATION
		cf = M_PI/(nlocl-1);
		for (i = 0; i < nlocl; i++){
			zlocl[nlocl - i - 1] = gsl_sf_cos(cf*i);
		}

		lagrange(nlocl-1, nquad, zlocl, zquad);

		/* NOTE: Calculating the Lagrange polynomials here
		 * is redundant; that is, they are calculated else-
		 * where for the same interval. Could improve per-
		 * formance by eliminating redundancy, but since
		 * this is a small calculation, we'll just keep it
		 * here for now */

		for (i = 0; i < nelem; i++){		// loop over boundary elements
			for (j = 0; j < nlocl; j++){	// loop over local basis nodes

			}
		}


    // NEED LAGRANGE INTERPOLANT AAAAAND THEIR DERIVATIVES

  }

};

#endif
