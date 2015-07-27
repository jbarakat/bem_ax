/* STOKES FLOW
 *  Boundary class that inherits the members of GEOM.H and, in addition,
 *  prescribes the kinematic (e.g., velocity) and dynamic (e.g., stress)
 *  quantities associated with a Stokes flow.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  u      [output]			displacement
 *  v      [output]			velocity
 *  S      [output]			Cauchy stress
 *  f      [output]			traction
 */

#ifndef STOKES_H
#define STOKES_H

/* HEADER FILES */
#include "geom.h"

class stokes: public geom {
private:
	double *dispx, *dispr;
	double *disps, *dispp;
	double *velx , *velr ;
	double *vels , *velp ;
	double *trctx, *trctr;
	double *trcts, *trctp;
public:

	/* PROTOTYPES */
	
	
	/* IMPLEMENTATIONS */

	// Constructors

	// Destructor

	// Set functions

	// Get functions

};

#endif
