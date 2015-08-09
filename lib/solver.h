/* BOUNDARY INTEGRAL SOLVER
 *  Solve the boundary integral equations and time-evolve the
 *  quantities associated with the geometric nodes (e.g., position,
 *  concentration).
 *
 * REFERENCES
 *  
 * PARAMETERS
 *  x	[input/output]	description
 *  y [input/output]	description
 */


/* FOR NOW, FOCUS ON SINGLE LAYER IMPLEMENTATION AND
 * THE IMPLEMENTATION OF THE FREE BOUNDARY PROBLEM
 * (I.E., IMPLEMENTING THE KINEMATIC CONDITION) 
 * JUST USE A FORWARD EULER SCHEME FOR NOW */




#ifndef SOLVER_H
#define SOLVER_H

/* HEADER FILES */
#include "quad.h"


/* PROTOTYPES */


/* IMPLEMENTATIONS */

void timeInt(int nstep, int nquad, double dt, surface Surface){

  /*-----------------------------------------------------*/
  /*----------------------- SETUP -----------------------*/
  /*-----------------------------------------------------*/

	// declare variables
	int    istep;
	int    i, j, k, m, n;
	int    nelem, ngeom;
	int    nlocl, nglob;

	double area, vlme;
	
	double * v ;
	double * vx, * vr, * vn;
	double *Dfx, *Dfr;

	double * x, * r;
	double *nx, *nr;

	int    IGF;
	int    ISURF;
	
  /* get number of boundary elements,
   * geometric nodes, and basis nodes */
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// allocate memory
	v   = (double*) malloc( 2*nglob * sizeof(double));
	vx  = (double*) calloc(   ngeom , sizeof(double));
	vr  = (double*) calloc(   ngeom , sizeof(double));
	vn  = (double*) calloc(   ngeom , sizeof(double));
	x   = (double*) malloc(   ngeom * sizeof(double));
	r   = (double*) malloc(   ngeom * sizeof(double));
	nx  = (double*) malloc(   ngeom * sizeof(double));
	nr  = (double*) malloc(   ngeom * sizeof(double));


  /*-----------------------------------------------------*/
  /*-------------------- INITIALIZE ---------------------*/
  /*-----------------------------------------------------*/
	
	IGF = 1;   // tube Green's function
	IGF = 0;   // free-space Green's function

	ISURF = 2; // vesicle
	ISURF = 0; // drop

	// set surface velocity
	for (i = 0; i < nglob; i++){
		Surface.setVel(i, 0., 0.);
	}
	
	// get initial node coordinates and normal
	Surface.getNode(x , r );
	Surface.getNrml(nx, nr);
	
	
  /*-----------------------------------------------------*/
  /*----------------- TIME INTEGRATION ------------------*/
  /*-----------------------------------------------------*/
	
	for (istep = 0; istep < nstep; istep++){
	
		/*-- Step 1: Solve the boundary integral equation. --*/
		singleLayer(IGF, nquad, Surface, v);

		// NOTE: FOR SINGLE LAYER POTENTIAL, THERE IS NO NEED
		// TO CALCULATE A MATRIX INVERSE


		// write to file before evolving system
		// ADD WRITE FUNCTION

	
		/*-- Step 2: Update surface velocity and boundary ---*
		 *---------- shape using the kinematic condition. ---*/
		
		for (i = 0; i < ngeom; i++){
			// get global index
			n = i*(nlocl - 1);
	
			// calculate surface velocity
			vx[i]  = v[2*n  ];
			vr[i]  = v[2*n+1];
			vn[i]  = vx[i]*nx[i] + vr[i]*nr[i];
			
			// advect geometric nodes using forward Euler scheme
			x [i] += nx[i]*vn[i]*dt;
			r [i] += nr[i]*vn[i]*dt;
			
			// NOTE: SHOULD ALSO USE BACKWARD EULER SCHEME TO CHECK
			// ERROR.

		}
	
		/*-- Step 3: Update surface fields. -----------------*/
		Surface.setGeomParams(ngeom, x, r);
		Surface.getNrml(nx, nr);
	
		// NOTE: FOR NOW, DON'T RECALCULATE TENSION
		// AND MOMENTS, BECAUSE WE'RE USING A DROP.
		
		// NOTE #2: ALSO, WOULD HAVE TO UPDATE SURFACTANT
		// CONCENTRATION FIELD HERE.
		
		
	// free memory	
	free(v );
	free(vx);
	free(vr);
	free(vn);
	free(x );
	free(r );
	free(nx);
	free(nr);
	}
}




#endif
