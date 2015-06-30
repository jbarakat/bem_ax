/* MAIN PROGRAM
 * Execute library functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "gauleg.h"
#include "bessel.h"
//#include <boost/math/special_functions>
//#include "grnfcn.h"

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

void testGauleg();
void testBessel();
void testGrnfcn();

int main(){
	//testGrnfcn();
	//testBessel();
	//testGauleg();

	double complex z;
	z = lapack_make_complex_double(2.4, 1.1);
	printf("%.4f\n",creal(z));

	return(0);
}

void testGrnfcn(){
	double Dk, dDkds;
	int k;
	double complex z;
	z = lapack_make_complex_double(2.4, 1.1);
	
	// get k
	printf("k = ");
	scanf("%d",&k);

	double s = 4;

//	calcDk(k, s, Dk, dDkds);
	printf("%.16f\n%.16f\n", Dk, dDkds);

	//int n = 1;
	//double an, bn, cn;
	//lapack_complex_double *xn, *yn;
	//calcDkRoots(n, k, an, bn, cn, xn, yn);
	//printf("an = %.16f\nbn = %.16f\ncn = %.16f\n", an, bn, cn);
}

void testBessel(){
	int n, nmin, nmax, *narray, nsize;
	double x, nu;
	double Jn, Yn, In, Kn;
	double *Jnarray, *Ynarray, *Inarray, *Knarray;

	// allocate memory
	narray = (int*) calloc(nsize,sizeof(int));
	Jnarray = (double*) calloc(nsize,sizeof(double));
	Ynarray = (double*) calloc(nsize,sizeof(double));
	Inarray = (double*) calloc(nsize,sizeof(double));
	Knarray = (double*) calloc(nsize,sizeof(double));
	
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
	printf("n = %.4f, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", nu, x, Jn, Yn, In, Kn);

	for (int i = 0; i < nsize; i++){
		narray[i] = i; 
	}
	nmin = narray[0];
	nmax = narray[nsize-1];
	
	besselJArray(nmin, nmax, x, Jnarray);
	besselYArray(nmin, nmax, x, Ynarray);
	besselIArray(nmin, nmax, x, Inarray);
	besselKArray(nmin, nmax, x, Knarray);
	printf("n  Jn      Yn      In      Kn\n");
	for (int i = 0; i < nsize; i++){
		printf("%d  %.4f  %.4f  %.4f  %.4f\n", narray[i], Jnarray[i], Ynarray[i], Inarray[i], Knarray[i]);
	}
	
	free(narray);
	free(Jnarray);
	free(Ynarray);
	free(Inarray);
	free(Knarray);
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
