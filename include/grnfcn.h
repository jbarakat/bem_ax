/* GREEN'S FUNCTIONS
 *  Evaluate axisymmetric Green's functions for Stokes flow.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [pp. 89-91]
 *  Tozeren, Inter. J. Num. Meth. Fluids 4, 159-170 (1984) 
 *  
 * PARAMETERS
 *  x,r   [input]   field point
 *  x0,r0 [input]   source point
 *  rc    [input]   cylindrical tube radius
 *  f     [input]   force density
 *  u     [output]  velocity
 *  M     [output]  Green's function
 */

#ifndef GRNFCN_H
#define GRNFCN_H


/* HEADER FILES */
#include "bessel.h"
#include "ellint.h"
#include <math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>

/* PROTOTYPES */
void gf_axR(double, double, double, double,
            double&, double&, double&, double&);
void gf_axR(double, double, double, double,
            double&, double&, double&, double&,
            double&, double&, double&, double&,
            double&, double&, double&, double&);
void gf_axR_vel(double, double, double, double,
                double, double, double&, double&);
void gf_axT(double, double, double, double, double,
            double&, double&, double&, double&);
void gf_axT(double, double, double, double, double,
            double&, double&, double&, double&,
            const int, const double);
void gf_axT_ker(double, double, double, double, double, double,
            double &, double &, double&, double &);
void gf_axT_vel(double, double, double, double, double,
                double, double, double&, double&);


#endif
