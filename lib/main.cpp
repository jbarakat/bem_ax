/* MAIN PROGRAM
 * Execute library functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include "gauleg.h"
#include "bessel.h"

void testGauleg();
void testBessel();

int main(){
	testBessel();
	//testGauleg();

	return(0);
}

void testBessel(){
	int n, nmin, nmax, *narray, nsize;
	double x, nu;
	double Jn, Yn, In, Kn;
	double *Jnarray, *Ynarray, *Inarray, *Knarray;

	x = 4.0;
	nsize = 10;
	
	n = 0;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	n = 1;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	n = 2;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	
	nu = 0.5;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	narray = (int*) calloc(nsize,sizeof(int));
	Jnarray = (double*) calloc(nsize,sizeof(double));
	Ynarray = (double*) calloc(nsize,sizeof(double));
	Inarray = (double*) calloc(nsize,sizeof(double));
	Knarray = (double*) calloc(nsize,sizeof(double));
	
	for (int i = 0; i < nsize; i++){
		narray[i] = i; 
	}
	nmin = narray[0];
	nmax = narray[nsize-1];
	
	Jnarray = besselJArray(nmin, nmax, x);
	Ynarray = besselYArray(nmin, nmax, x);
	Inarray = besselIArray(nmin, nmax, x);
	Knarray = besselKArray(nmin, nmax, x);
	printf("n  Jn      Yn      In      Kn\n");
	for (int i = 0; i < nsize; i++){
		printf("%d  %.4f  %.4f  %.4f  %.4f\n", narray[i], Jnarray[i], Ynarray[i], Inarray[i], Knarray[i]);
	}
}

void testGauleg(){
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
}
