/* NEWTON'S METHOD
 *  Newton-Raphson method with bisection to guide the solution.
 *
 * REFERENCES
 *  Press et al, Cambridge University Press (2007)
 *  
 * PARAMETERS
 */

#ifndef NEWTON_H
#define NEWTON_H

/* HEADER FILES */

/* PROTOTYPES */


/* IMPLEMENTATIONS */
/***********************************************************************/
template <class T>
double newt(T &f, T &dfdx, const double XMIN, const double XMAX){
	const int MAXIT = 100;
	const double XACC = 1e-14;


}


/***********************************************************************/

#endif
