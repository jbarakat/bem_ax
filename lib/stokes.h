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
		velx  = (double*) calloc(n, sizeof(double));
		velr  = (double*) calloc(n, sizeof(double));
		trctx = (double*) calloc(n, sizeof(double));
		trctr = (double*) calloc(n, sizeof(double));
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
		velx  = (double*) calloc(n, sizeof(double));
		velr  = (double*) calloc(n, sizeof(double));
		trctx = (double*) calloc(n, sizeof(double));
		trctr = (double*) calloc(n, sizeof(double));
		conc  = (double*) calloc(n, sizeof(double));

		// initialize viscosity ratio
		visc = lamb;
	}

	// Destructor

	// Set functions

	// Get functions
	void getTrct(int i, double &fx, double &fr){
		int n = getNNode();

		if (i >= n){
			printf("Error: index out of bounds.\n");
			return;
		}
		
		fx = trctx[i];
		fr = trctr[i];
	}

	void getTrct(double *fx, double *fr){
		int i;
		int n = getNNode();

		if (fx == NULL || fr == NULL){
			printf("Error: no memory allocated for fx, fr\n");
		}

		for (i = 0; i < n-1; i++){
			fx[i] = trctx[i];
			fr[i] = trctr[i];
		}
		
	}

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
