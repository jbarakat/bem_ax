/* MAIN PROGRAM
 * Execute library functions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "gauleg.h"
#include "bessel.h"
#include "ellint.h"
//#include <boost/math/special_functions>
#include "grnfcn.h"
#include "interp.h"
#include "quad.h"
#include "surface.h"
#include <math.h>
#include <vector>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif


void testSingleLayer();
void testQuadrature();
void testLagrange();
void testGeom();
void testInterp();

int main(){
	testSingleLayer();
//	testQuadrature();
//	testLagrange();
//	testGeom();
//	testInterp();

	return(0);
}

void testSingleLayer(){
	int i, j, k;
	int IGF = 0;
	int N = 2;
	int M = 1;
	int nquad;
	double lamb = 1;
	double *x, *r, *thet;
	double a, b;
	double *A, *df, *v;

	// interface parameters
	double gamm = 0.2;
	double ES = 0.;
	double ED = 0.;
	double EB = 0.;
	double ET = 0.;

	// constitutive model
	int model = 0;

	// get number of collocation points
	int ncoll = N*M + 1;

	// get number of quadrature points
	printf("nquad = ");
	scanf("%u", &nquad);

	// allocate memory
	x     = (double*) malloc((N+1)         * sizeof(double));
	r     = (double*) malloc((N+1)         * sizeof(double));
	thet  = (double*) malloc((N+1)         * sizeof(double));
	A     = (double*) malloc(4*ncoll*ncoll * sizeof(double));
	df    = (double*) malloc(2*ncoll       * sizeof(double));
	v     = (double*) malloc(2*ncoll       * sizeof(double));

	// define coordinates on an ellipse
	a = 1.;
	b = 1.;
	for (i = 0; i < N+1; i++){
		thet[i] = (N-i)*M_PI/N;
		x[i] = a*gsl_sf_cos(thet[i]);
		r[i] = b*gsl_sf_sin(thet[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	surface spheroid(model, N, M, lamb, gamm,
	                 ES, ED, EB, ET,
									 x, r);

	singleLayer(IGF, nquad, spheroid, A, df, v);

	printf("A = \n");
	for (i = 0; i < 2*ncoll; i++){
		for (j = 0; j < 2*ncoll; j++){
			if (A[2*ncoll*i + j] >= 0 && A[2*ncoll*i + j] < 10)
				printf(" ");
			printf("%.4f ", A[2*ncoll*i + j]);
		}
		printf("\n");
	}
	printf("\n");
	
	printf("df = \n");
	for (i = 0; i < ncoll; i++){
		for (j = 0; j < 2; j++){
			if (df[2*i+j] >= 0 && df[2*i+j] < 10)
				printf(" ");
			printf("%.4f", df[2*i+j]);
			printf(" ");
		}
		printf("\n");
	}
	printf("\n");
	
	printf("v = \n");
	for (i = 0; i < ncoll; i++){
		for (j = 0; j < 2; j++){
			if (v[2*i+j] >= 0 && v[2*i+j] < 10)
				printf(" ");
			printf("%.4f", v[2*i+j]);
			printf(" ");
		}
		printf("\n");
	}

	double *ax, *bx, *cx;
	double *ar, *br, *cr;
	ax = (double*) malloc( N    * sizeof(double));
	bx = (double*) malloc((N+1) * sizeof(double));
	cx = (double*) malloc( N    * sizeof(double));
	ar = (double*) malloc( N    * sizeof(double));
	br = (double*) malloc((N+1) * sizeof(double));
	cr = (double*) malloc( N    * sizeof(double));

	spline(N, x, r, 0, 0, 1, -1, ax, bx, cx, ar, br, cr);
	printf("a      b      c");
	printf("\n");
	for (i = 0; i < N+1; i++){
		printf("%.4f ", ar[i]);
		printf("%.4f ", br[i]);
		printf("%.4f ", cr[i]);
		printf("\n");
	}
}

void testQuadrature(){
	int i, j, k;
	double *x, *r, *thet;
	double *P;
	double *zq, *wq;
	double z, w;
	double h;
	double integral;
	int nquad;

	// get number of quadrature points
	printf("nquad = ");
	scanf("%u", &nquad);

	// allocate memory
	x     = (double*) malloc(nquad * sizeof(double));
	r     = (double*) malloc(nquad * sizeof(double));
	P     = (double*) malloc(nquad * sizeof(double));
	thet  = (double*) malloc(nquad * sizeof(double));
	zq    = (double*) malloc(nquad * sizeof(double));
	wq    = (double*) malloc(nquad * sizeof(double));
	
	// get quadrature points
	gauleg(nquad, zq, wq);

	// define coordinates on a circle
	for (i = 0; i < nquad; i++){
		thet[i] = M_PI/4*(1 + zq[i]);
		x[i] = gsl_sf_cos(thet[i]);
		r[i] = gsl_sf_sin(thet[i]);
	}

	// evaluate quadrature
	integral = 0;
	for (i = 0; i < nquad; i++){
		integral += wq[i]*x[i];
	}
	integral *= M_PI/4;
	
	printf("Exact answer = 1\n");
	printf("Approximate answer = %.4f\n", integral);


}

void testLagrange(){
	int i, j, k;
	int N = 10; // order of Lagrange polynomial
	int M = 10; // number of points to interpolate to
	double lamb = 1;
	double *x, *r, *thet;
	double *xi;
	double *theti;
	double a, b;
	double *L, *L1;
	double *dLdx, *dL1dx;

	// allocate memory
	x     = (double*) malloc((N+1)   * sizeof(double));
	xi    = (double*) malloc(      M * sizeof(double));
	r     = (double*) malloc((N+1)   * sizeof(double));
	thet  = (double*) malloc((N+1)   * sizeof(double));
	theti = (double*) malloc((N+1)   * sizeof(double));
	L     = (double*) malloc((N+1)*M * sizeof(double));
	L1    = (double*) malloc((N+1)   * sizeof(double));
	dLdx  = (double*) malloc((N+1)*M * sizeof(double));
	dL1dx = (double*) malloc((N+1)   * sizeof(double));

	// define coordinates on an ellipse
	a = 1.;
	b = 1.;
	printf("x = \n");
	for (i = 0; i < N+1; i++){
		thet[i] = i*M_PI/N;
		x[i] = a*gsl_sf_cos(thet[i]);
		r[i] = b*gsl_sf_sin(thet[i]);
		printf("%.4f\n", x[i]);
	}
	printf("\n");

	printf("theti = \n");
	for (i = 0; i < M; i++){
		xi   [i] = x   [0] + i*(x   [N] - x   [0])/(M-1);
		theti[i] = thet[0] + i*(thet[N] - thet[0])/(M-1);
	//	xi[i] *= 0.9;
	//	printf("%.4f\n", xi[i]);
		printf("%.4f\n", theti[i]);
	}
	printf("\n");

//	lagrange(N, M, x, xi, L);
//	
//	printf("L = \n");
//	for (i = 0; i < N+1; i++){		// rows = polynomial
//		for (j = 0; j < M; j++){		// columns = grid point
//			if (L[i*M+j] >= 0)
//				printf(" ");
//			printf("%.4f ", L[i*M + j]);
//
//	//		if (L[i*M+j] < 10)
//	//			printf(" ");
//	//		printf("%d ",i*M + j);
//		}
//		printf("\n");
//	}
//	
//	lagrange(N, x, xi[1], L1);
//	printf("L = \n");
//	for (i = 0; i < N+1; i++){		// polynomial
//		if (L1[i] >= 0)
//			printf(" ");
//		printf("%.4f", L1[i]);
//		printf("\n");
//	}

	lagrange(N, M, x, xi, L, dLdx);
	
	printf("L = \n");
	for (i = 0; i < N+1; i++){		// rows = polynomial
		for (j = 0; j < M; j++){		// columns = grid point
			if (L[i*M+j] >= 0)
				printf(" ");
			printf("%.4f ", L[i*M + j]);
		}
		printf("\n");
	}
	
	printf("dLdx = \n");
	for (i = 0; i < N+1; i++){		// rows = polynomial
		for (j = 0; j < M; j++){		// columns = grid point
			if (dLdx[i*M+j] >= 0)
				printf(" ");
			printf("%.4f ", dLdx[i*M + j]);
		}
		printf("\n");
	}
	
	for (k = 0; k < M; k++){
		lagrange(N, x, xi[k], L1, dL1dx);
		for (i = 0; i < N+1; i++){		// polynomial
			printf(" L%d(xi[%d]) = ",i, k);
			if (L1[i] >= 0)
				printf(" ");
			printf("%.4f    ", L1[i]);
		}
		printf("\n");
		
		for (i = 0; i < N+1; i++){		// polynomial
			printf("dL%d(xi[%d]) = ", i, k);
			if (dL1dx[i] >= 0)
				printf(" ");
			printf("%.4f    ", dL1dx[i]);
		}
		printf("\n");
	}
	
	lagrange(N, M, thet, theti, L, dLdx);
	
	double yi, dyi;
	for (k = 0; k < M; k++){
		lagrange(N, thet, theti[k], L1, dL1dx);
		yi  = 0;
		dyi = 0;
		for (i = 0; i < N+1; i++){
			yi  += L1   [i]*r[i];
			dyi += dL1dx[i]*r[i];
		}
		printf("yint = %.4f, yexact = %.4f   ", yi, gsl_sf_sin(theti[k]));
		printf("dyint = %.4f, dyexact = %.4f", dyi, gsl_sf_cos(theti[k]));
		printf("\n");
		
	}
}

void testGeom(){
	int i, j;
	double *x, *r, *thet;
	double *l, *s;
	double V, A;
	double V0, A0;
	double relerrA, relerrV;
	double *ks, *kp, *nx, *nr, *tx, *tr;
	double *ax, *bx, *cx, *ar, *br, *cr;
	double a, b;
	int N = 40;
	
	// allocate memory
	ax    = (double*) malloc( N    * sizeof(double));
	bx    = (double*) malloc((N+1) * sizeof(double));
	cx    = (double*) malloc( N    * sizeof(double));
	ar    = (double*) malloc( N    * sizeof(double));
	br    = (double*) malloc((N+1) * sizeof(double));
	cr    = (double*) malloc( N    * sizeof(double));
	x     = (double*) malloc((N+1) * sizeof(double));
	r     = (double*) malloc((N+1) * sizeof(double));
	thet  = (double*) malloc((N+1) * sizeof(double));
	l     = (double*) malloc((N+1) * sizeof(double));
	s     = (double*) malloc((N+1) * sizeof(double));
	ks    = (double*) malloc((N+1) * sizeof(double));
	kp    = (double*) malloc((N+1) * sizeof(double));
	tx    = (double*) malloc((N+1) * sizeof(double));
	tr    = (double*) malloc((N+1) * sizeof(double));
	nx    = (double*) malloc((N+1) * sizeof(double));
	nr    = (double*) malloc((N+1) * sizeof(double));
	
	// define coordinates on an ellipse
	a = 10.;
	b = 1.;
	for (i = 0; i < N+1; i++){
		thet[i] = i*M_PI/N;
		x[i] = a*gsl_sf_cos(thet[i]);
		r[i] = b*gsl_sf_sin(thet[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	// calculate analytical area and volume
	double e2, e1;
	double asine, atanhe;
	A0 = 1.;
	if (a < b){ // oblate
		e2 = 1. - a*a/(b*b);
		e1 = sqrt(e2);
		atanhe = gsl_complex_abs(gsl_complex_arctanh_real(e1));
		A0 += ((1 - e2)/e1)*atanhe;
		A0 *= 2*M_PI*b*b;
	}

	if (a > b){ // prolate
		e2 = 1. - b*b/(a*a);
		e1 = sqrt(e2);
		asine = gsl_complex_abs(gsl_complex_arcsin_real(e1));
		A0 += (a/(b*e1))*asine;
		A0 *= 2*M_PI*b*b;

	}
	
	if (a == b){
		A0 = 4*M_PI*a*a;
	}

	V0 = 4*M_PI*a*b*b/3;

	// test constructor, set and get functions
	stokes spheroid(0, N, 1, x, r);
//	spheroid.getArcl(s);
//	spheroid.getArea(A); //	OR A = spheroid.getArea();
//	spheroid.getVlme(V); //	OR V = spheroid.getVlme();
//	spheroid.getCurv(ks, kp);
//	spheroid.getTang(tx, tr);
//	spheroid.getNrml(nx, nr);
	spheroid.getAll(N, x, r, ax, bx, cx, ar, br, cr, l, s, A, V, ks, kp, tx, tr, nx, nr);

	relerrA = (A-A0)/A0*100;
	relerrV = (V-V0)/V0*100;

	printf("theta  l       s       ks      kp      tx      tr      nx      nr\n");
	for (i = 0; i < N+1; i++)
		printf("%.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f\n", thet[i], l[i], s[i], ks[i], kp[i], tx[i], tr[i], nx[i], nr[i]);

	printf("\n a = %.1f, b = %.1f\n", a, b);
	printf("\n A = %.4f, A0 = %.4f, relerr = %.4f%\n", A, A0, relerrA);
	printf("\n V = %.4f, V0 = %.4f, relerr = %.4f%\n", V, V0, relerrV);
	
}

void testInterp(){
	const int N = 10;


	int i, j, k;
	double *X, *Y, *T;
	double x, y;
	double dx;
	double a, b;
	double *L;
	double *A, *B, *C;
	double *Ax, *Bx, *Cx;
	double *Ay, *By, *Cy;
	
	// allocate memory
	X  = (double*) malloc((N+1) * sizeof(double));
	Y  = (double*) malloc((N+1) * sizeof(double));
	T  = (double*) malloc((N+1) * sizeof(double));
	L  = (double*) malloc((N+1) * sizeof(double));
	A  = (double*) malloc( N    * sizeof(double));
	B  = (double*) malloc((N+1) * sizeof(double));
	C  = (double*) malloc( N    * sizeof(double));
	Ax = (double*) malloc( N    * sizeof(double));
	Bx = (double*) malloc((N+1) * sizeof(double));
	Cx = (double*) malloc( N    * sizeof(double));
	Ay = (double*) malloc( N    * sizeof(double));
	By = (double*) malloc((N+1) * sizeof(double));
	Cy = (double*) malloc( N    * sizeof(double));
	
	// define coordinates on an ellipse
	a = 10.;
	b = 1.;
	for (i = 0; i < N+1; i++){
		T[i] = i*M_PI/N;
		X[i] = a*gsl_sf_cos(T[i]);
		Y[i] = b*gsl_sf_sin(T[i]);
		//printf("%.4f %.4f %.4f\n", T[i], X[i], Y[i]);
	}

	// interpolate (Lagrange)
	x = 0.6;
	lagrange(N, X, Y, x, y);
	printf("y = %.4f\n", y);

	lagrange(N, X, x, L);
	printf("L = \n");
	for (i = 0; i < N+1; i++){
		printf("%.4f\n", L[i]);
	}

	// interpolate (cubic spline)
//	spline(N, X, Y, slope1, slope2, x, y);

	double slopex1 = 0;
	double slopex2 = 0;
	double slopey1 = 1.;
	double slopey2 = -1.;
	spline(N, X, Y, slopex1, slopex2, slopey1, slopey2, Ax, Bx, Cx, Ay, By, Cy);
	
	printf("Ax        Bx        Cx\n");
	for (i = 0; i < N+1; i++){
		printf("%.4f %.4f %.4f\n", Ax[i], Bx[i], Cx[i]);
	}
	printf("Ay        By        Cy\n");
	for (i = 0; i < N+1; i++){
		printf("%.4f %.4f %.4f\n", Ay[i], By[i], Cy[i]);
	}
	
	free(X);
	free(Y);
	free(T);
	free(A);
	free(B);
	free(C);
	free(Ax);
	free(Bx);
	free(Cx);
	free(Ay);
	free(By);
	free(Cy);
}
