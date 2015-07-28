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
 *  u      [output]			displacement
 *  v      [output]			velocity
 *  f      [output]			traction
 *  c		   [output]			concentration (or any scalar field)
 *  lamb   [output]			viscosity ratio
 */

#ifndef STOKES_H
#define STOKES_H

/* HEADER FILES */
#include "geom.h"

const int RIGID = 0;
const int FLUID = 1;

class stokes: public geom {
private:
	// indicator for type of boundary
	int     IBE;
		/*     = RIGID  (rigid boundary)
		 *     = FLUID  (fluid-fluid interface) */
	
	// displacement
	double *dispx, *dispr;
	double *disps, *dispp;

	// velocity
	double *velx , *velr ;
	double *vels , *velp ;

	// traction
	double *trctx, *trctr;
	double *trcts, *trctp;

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

	stokes(int ibe, int N, double *x, double *r) : geom(N, x, r) {
		if (ibe != 0 && ibe != 1){
			printf("Error: ibe can only take values of 0 or 1.");
			return;
		}

		int i, n;

		// set indicator function
		IBE = ibe;
		
		// get number of boundary nodes
		n = getNNode();

		/* allocate memory for pointer arrays
		 * and initialize to zero */
		dispx = (double*) calloc(n, sizeof(double));
		dispr = (double*) calloc(n, sizeof(double));
		disps = (double*) calloc(n, sizeof(double));
		dispp = (double*) calloc(n, sizeof(double));
		velx  = (double*) calloc(n, sizeof(double));
		velr  = (double*) calloc(n, sizeof(double));
		vels  = (double*) calloc(n, sizeof(double));
		velp  = (double*) calloc(n, sizeof(double));
		trctx = (double*) calloc(n, sizeof(double));
		trctr = (double*) calloc(n, sizeof(double));
		trcts = (double*) calloc(n, sizeof(double));
		trctp = (double*) calloc(n, sizeof(double));
		conc  = (double*) calloc(n, sizeof(double));

		// initialize viscosity ratio
		visc = 1.;
	}

	stokes(int ibe, int N, double lamb, double *x, double *r) : geom(N, x, r) {
		if (ibe != 0 && ibe != 1){
			printf("Error: ibe can only take values of 0 or 1.");
			return;
		}

		int i, n;
		
		// set indicator function
		IBE = ibe;
		
		// get number of boundary nodes
		n = getNNode();

		/* allocate memory for pointer arrays
		 * and initialize to zero */
		dispx = (double*) calloc(n, sizeof(double));
		dispr = (double*) calloc(n, sizeof(double));
		disps = (double*) calloc(n, sizeof(double));
		dispp = (double*) calloc(n, sizeof(double));
		velx  = (double*) calloc(n, sizeof(double));
		velr  = (double*) calloc(n, sizeof(double));
		vels  = (double*) calloc(n, sizeof(double));
		velp  = (double*) calloc(n, sizeof(double));
		trctx = (double*) calloc(n, sizeof(double));
		trctr = (double*) calloc(n, sizeof(double));
		trcts = (double*) calloc(n, sizeof(double));
		trctp = (double*) calloc(n, sizeof(double));
		conc  = (double*) calloc(n, sizeof(double));

		// initialize viscosity ratio
		visc = lamb;
	}

	// Destructor

	// Set functions

	// Get functions
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
