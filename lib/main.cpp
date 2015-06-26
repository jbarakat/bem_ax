/* MAIN PROGRAM
 * Execute library functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include "gauleg.h"

int main(){
	// declare variables
	int n;
	double *X, *W;
	
	printf("N = ");
	scanf("%u",&n);

	// note: have to make sure the range of integration is from -1 to 1, see Numerical Recipes p. 207 (p. 184 on pdf)

	// allocate memory
	X = (double*) calloc(n,sizeof(double));
	W = (double*) calloc(n,sizeof(double));

	gauleg(n, X, W);

	// print abscissas and weights
	int i;
	printf("\nX = ");
	for (i = 0; i < n; i++)
		printf("%.6f ",X[i]);
	printf("\n");

	printf("\nW = ");
	for (i = 0; i < n; i++)
		printf("%.6f ",W[i]);
	printf("\n");

	return(0);
}
