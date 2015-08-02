/* SURFACE MECHANICS
 *  Produce the local in-plane tensions, transverse shear tension,
 *  and bending/twisting moments for an axisymmetric free surface.
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

#include "stokes.h"

/* HEADER FILES */


/* PROTOTYPES */


/* IMPLEMENTATIONS */

void drop(
     stokes Stokes,
     double taus, double taup,
	 	 double q   , 
	 	 double ms  , double mp){
	// declare variables
}


void helfrich(
     stokes Stokes,
     double taus, double taup,
	 	 double q   , 
	 	 double ms  , double mp){
	// declare variables
}

void skalak(
     stokes Stokes,
     double taus, double taup,
	 	 double q   , 
	 	 double ms  , double mp){
	// declare variables
}



#endif
