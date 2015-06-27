/* GAUSS-LEGENDRE QUADRATURE
 *  Generate the abscissas and weigts for Gauss-Legendre quadrature by
 *  solving an eigenvalue proglem for a symmetric, tridiagonal matrix.
 *
 * REFERENCES
 *  Golub and Welsch, Mathematics of Computation 23-106 (1969)
 * 
 * PARAMETERS
 *  n	[input]		number of points
 *  X	[output]	abcsissas
 *  W	[output]	weights
 */

#ifndef GAUSSLEG_H
#define GAUSSLEG_H

/* HEADER FILES */
#include <lapacke.h>
#include <cblas.h>
#include <math.h>

/* PROTOTYPES */
void gauleg(int, double*, double*);

#endif
