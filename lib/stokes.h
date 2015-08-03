/* STOKES FLOW
 *  Boundary class that inherits the members of GEOM.H and, in addition,
 *  prescribes the material properties (e.g., viscosity) as well as the
 *  kinematic (e.g., velocity) and dynamic (e.g., stress) quantities 
 *  associated with the boundary in a Stokes flow.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  x,r    [input]			nodal coordinates
 *  N      [input]			number of (global) elements
 *  M      [input]			number of (local) subelements
 *  lamb   [input]			viscosity ratio
 *  u      [output]			displacement
 *  v      [output]			velocity
 *  f      [output]			traction
 *  c		   [output]			concentration (or any scalar field)
 */

#ifndef STOKES_H
#define STOKES_H

/* HEADER FILES */
#include "geom.h"

const int RIGID = 0;
const int FLUID = 1;

class stokes: public geom {
friend class surface;
private:
	// indicator for type of boundary
	int     type ;
	/*       = RIGID  (rigid boundary)
	 *       = FLUID  (fluid-fluid interface) */

	/* local and global number of basis nodes 
	 * for density function interpolation */
	int     nlocl,  nglob;
	
	// displacement
	double *dispx, *dispr;

	// velocity
	double *velx , *velr ;

	// traction
	double *trctx, *trctr;

	// concentration
	double *conc ;

	// viscosity ratio
	double  visc ;
	/* NOTE: By convention, the denominator
	 *       corresponds to the phase into
	 *       which the normal vector points. */

public:

	/* PROTOTYPES */
	
	
	/* IMPLEMENTATIONS */

	// Constructors
	stokes() : geom() {
	}

	stokes(int id, int N, int M, 
	       double *x, double *r) : geom(N, x, r) {
		// set viscosity ratio = 1 if not given
		stokes(id, N, M, 1., x, r);
	}

	stokes(int id, int N, int M, 
	       double lamb, double *x, double *r) : geom(N, x, r) {
		if (id != 0 && id != 1){
			printf("Error: id can only take values of 0 or 1.\n");
			return;
		}
		
		if (M < 1){
			printf("Error: cannot have fewer than 1 subelements.\n");
			return;
		}

		int i, n;

		// set boundary type
		type = id;

		// set local and global number of basis nodes
		nlocl = M + 1;
		nglob = N*M + 1;

		/* allocate memory for pointer arrays
		 * and initialize to zero */
		dispx = (double*) calloc(nglob, sizeof(double));
		dispr = (double*) calloc(nglob, sizeof(double));
		velx  = (double*) calloc(nglob, sizeof(double));
		velr  = (double*) calloc(nglob, sizeof(double));
		trctx = (double*) calloc(nglob, sizeof(double));
		trctr = (double*) calloc(nglob, sizeof(double));
		conc  = (double*) calloc(nglob, sizeof(double));

		// initialize viscosity ratio
		visc = lamb;
	}

	// Destructor

	/*- Set functions ---*/

	/*- Get functions ---*/
	// get local number of basis nodes
	int getNLocl(){
		int n = nlocl;
		return(n);
	}
	
	void getNLocl(int n){
		n = nlocl;
	}

	// get global number of basis nodes
	int getNGlob(){
		int n = nglob;
		return(n);
	}
	
	void getNGlob(int n){
		n = nglob;
	}

	/* get traction of the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void getTrct(int ielem, int ilocl, double &fx, double &fr){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds.\n");
			return;
		}

		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;

		fx = trctx[iglob];
		fr = trctr[iglob];
	}

	// get traction at the (iglob)th global basis node
	void getTrct(int iglob, double &fx, double &fr){
		if (iglob >= nglob){
			printf("Error: index out of bounds.\n");
			return;
		}
		
		fx = trctx[iglob];
		fr = trctr[iglob];
	}

	// get traction at all global basis nodes
	void getTrct(double *fx, double *fr){
		int iglob;

		if (fx == NULL || fr == NULL){
			printf("Error: no memory allocated for fx, fr.\n");
		}

		for (iglob = 0; iglob < nglob; iglob++){
			fx[iglob] = trctx[iglob];
			fr[iglob] = trctr[iglob];
		}
	}

	// get viscosity ratio
	double getVisc(){
		double lamb;
		lamb = visc;
		return(lamb);
	}
	
	void getVisc(double &lamb){
		lamb = visc;
	}


};

#endif
