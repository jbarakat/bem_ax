/* GAUSS-LOGARITHMIC QUADRATURE
 *  Table of quadrature points that is exact for integrands of the
 *  form f(x) + g(x)*log(x), where f and g are polynomials
 *
 * REFERENCES
 *  Crow, Mathematics of Computation 60-201 (1993)
 * 
 * PARAMETERS
 *  n	[input]		number of points
 *  X	[output]	abcsissas
 *  W	[output]	weights
 */

#ifndef GAULOG_H
#define GAULOG_H

/* HEADER FILES */
#include <stdlib.h>
#include <stdio.h>

/* PROTOTYPES */
void gaulog(int n, double * X, double * W);

#endif

