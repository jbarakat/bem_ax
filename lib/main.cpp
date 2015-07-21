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

// not sure if these declarations are neceessary...
#ifndef lapack_complex_float
#define lapack_complex_float float complex
#endif

#ifndef lapack_complex_double
#define lapack_complex_double double complex
#endif

void testGauleg();
void testBessel();
void testBesselComplex();
void testBesselNegativeOrder();
void testBesselComplexNegativeOrder();
void testEllint();
void testGrnfcnR();
void testGrnfcnT();

int main(){
//	testGrnfcnR();
	testGrnfcnT();
//	testBessel();
//	testBesselComplex();
//	testBesselNegativeOrder();
//	testBesselComplexNegativeOrder();
//	testEllint();
//	testGauleg();

	return(0);
}

void testGrnfcnT(){
	// declare variables
	int i, j, k;
	double x, x0, r, r0, rc;
	double xmin, xmax, rmin, rmax;
	int Nx, Nr;
	double *M, *u, *f;
	double Mxx, Mxr, Mrx, Mrr;
	double fx, fr, ux, ur;
	
	double xring = 0.;
	double yring = 1.;
	int istr, istep;
	double *xstream, *ystream;
	double *xstr, *ystr;
	double xvel, yvel;
	double velx, velx1, vely, vely1;
	int Nstr, Nstep;
	double Dr, Dt;

	// assign variables
	x = 0.;
	x0 = 0.1;
	r = 1.;
	r0 = 0.1;
	rc = 1.;

	fx = 1.; // axial point force
	fr = 0.; // radial point force
	
	/* complementary Green's function */
	// allocate memory
	M = (double*) malloc(2 * 2 * sizeof(double));
	u = (double*) malloc(2 * sizeof(double));

	// calculate Green's function
	gf_axT(x, r, x0, r0, rc, Mxx, Mxr, Mrx, Mrr);
	M[0] = Mxx/(8*M_PI);
	M[1] = Mxr/(8*M_PI);
	M[2] = Mrx/(8*M_PI);
	M[3] = Mrr/(8*M_PI);
	printf("M = \n");
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			printf("%.16f ", M[i*2 + j]);
		}
		printf("\n");
	}
	
	// calculate the velocity at x, r
	gf_axT_vel(x, r, x0, r0, rc, fx, fr, ux, ur);
	u[0] = ux;
	u[1] = ur;
	printf("\n");
	printf("u = \n");
	for (i = 0; i < 2; i++){
		printf("%.16f ", u[i]);
		printf("\n");
	}

	free(M);
	free(u);


	/* calculate velocity field and plot streamlines
	 * (based on the sgf_ax_fs_str.f Fortran code from Pozrikidis) */
	FILE *pFile;
	pFile =  fopen("str.dat", "w+");
	fprintf(pFile, "xstr ystr\n");
	if (fx == 1. && fr == 0. || fx == 0. && fr == 1.){
		if (fx == 1. && fr == 0.){
			Nstr = 20;
			Nstep = 2*128;
			Dr = 0.02;

			xstream = (double*) malloc(Nstr * sizeof(double));
			ystream = (double*) malloc(Nstr * sizeof(double));

			for (istr = 0; istr < Nstr; istr++){
				xstream[istr] = 0.01;
				ystream[istr] = 0.10*(istr+1);
			}
			ystream[19] = 1.95;
			
		}
		else if (fx == 0. && fr == 1.){
			Nstr = 13;
			Nstep = 2*128;
			Dr = -0.02;

			xstream = (double*) malloc(Nstr * sizeof(double));
			ystream = (double*) malloc(Nstr * sizeof(double));

			for (istr = 0; istr < Nstr; istr++){
				xstream[istr] = -0.1 + 0.1*istr;
				ystream[istr] = 1.99;
			}
			xstream[0] = 0.01;
			xstream[1] = 0.05;
			
		}

		for (istr = 0; istr < Nstr; istr++){
			printf("%d\n",istr);

			// allocate memory
			xstr = (double*) malloc(1 * sizeof(double));
			ystr = (double*) malloc(1 * sizeof(double));
		
			// initialize streamline
			xstr[0] = xstream[istr];
			ystr[0] = ystream[istr];

			for (istep = 0; istep < Nstep; istep++){
				fprintf(pFile, "%.8f %.8f\n", xstr[istep], ystr[istep]);

				// update field point
				xvel = xstr[istep];
				yvel = ystr[istep];

				// reallocate memory
				xstr = (double*) realloc(xstr, (istep + 2) * sizeof(double));
				ystr = (double*) realloc(ystr, (istep + 2) * sizeof(double));
				
				// calculate the velocity induced by a ring of point forces
				gf_axT_vel(xvel, yvel, xring, yring, 3., fx, fr, velx, vely);
				Dt = Dr/sqrt(velx*velx + vely*vely);

				//  update field point
				xvel = xstr[istep] + velx*Dt;
				yvel = ystr[istep] + vely*Dt;

				// calculate the induced velocity at the new field point
				gf_axT_vel(xvel, yvel, xring, yring, 3., fx, fr, velx1, vely1);
				
				// update the streamline using a simple midpoint rule for integration
				xstr[istep + 1] = xstr[istep] + 0.5*(velx + velx1)*Dt;
				ystr[istep + 1] = ystr[istep] + 0.5*(vely + vely1)*Dt;

				// condition for breaking for loop
				if (xstr[istep + 1] > 2)
					break;
				if (xstr[istep + 1] < -2)
					break;
				if (ystr[istep + 1] > 2)
					break;
				if (ystr[istep + 1] < -2)
					break;
			}
			fprintf(pFile, "\n");

			free(xstr);
			free(ystr);
		}
	}

	fclose(pFile);

}


void testGrnfcnR(){
	// declare variables
	int i, j, k;
	double x, x0, r, r0;
	double xmin, xmax, rmin, rmax;
	int Nx, Nr;
	double *M, *Q, *f, *u;
	double Mxx, Mxr, Mrx, Mrr;
	double Qxxx, Qxxr, Qxrx, Qxrr, Qrxx, Qrxr, Qrrx, Qrrr;
	double fx, fr, ux, ur;
	
	double xring = 0.;
	double yring = 1.;
	int istr, istep;
	double *xstream, *ystream;
	double *xstr, *ystr;
	double xvel, yvel;
	double velx, velx1, vely, vely1;
	int Nstr, Nstep;
	double Dr, Dt;

	// assign variables
	x = 2.;
	x0 = 1.;
	r = 2.;
	r0 = 1.;

	fx = 1.; // axial point force
	fr = 0.; // radial point force

	/* Green's function for a ring of point forces in free space */
	// allocate memory
	M = (double*) malloc(2 * 2 * sizeof(double));
	Q = (double*) malloc(2 * 2 * 2 * sizeof(double));
	u = (double*) malloc(2 * sizeof(double));

	// calculate stokeslet
	gf_axR(x, r, x0, r0, Mxx, Mxr, Mrx, Mrr);
	M[0] = Mxx;
	M[1] = Mxr;
	M[2] = Mrx;
	M[3] = Mrr;
	printf("M = \n");
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			printf("%.16f ", M[i*2 + j]);
		}
		printf("\n");
	}

	// calculate stokeslet and stresslet
	gf_axR(x, r, x0, r0, Mxx, Mxr, Mrx, Mrr, Qxxx, Qxxr, Qxrx, Qxrr, Qrxx, Qrxr, Qrrx, Qrrr);
	M[0] = Mxx;
	M[1] = Mxr;
	M[2] = Mrx;
	M[3] = Mrr;

	Q[0] = Qxxx;
	Q[1] = Qxxr;
	Q[2] = Qxrx;
	Q[3] = Qxrr;
	Q[4] = Qrxx;
	Q[5] = Qrxr;
	Q[6] = Qrrx;
	Q[7] = Qrrr;

	printf("M = \n");
	for (i = 0; i < 2; i++){
		for (j = 0; j < 2; j++){
			printf("%.16f ", M[i*2 + j]);
		}
		printf("\n");
	}

	printf("Q = \n");
	for (i = 0; i < 8; i++){
		printf("%.16f ", Q[i]);
		printf("\n");
	}

	// calculate the velocity at x, r
	gf_axR_vel(x, r, x0, r0, fx, fr, ux, ur);
	u[0] = ux;
	u[1] = ur;
	printf("\n");
	printf("u = \n");
	for (i = 0; i < 2; i++){
		printf("%.16f ", u[i]);
		printf("\n");
	}
	
	free(M);
	free(Q);
	free(u);
	

	/* calculate velocity field and plot streamlines
	 * (based on the sgf_ax_fs_str.f Fortran code from Pozrikidis) */
	FILE *pFile;
	pFile =  fopen("str.dat", "w+");
	fprintf(pFile, "xstr ystr\n");
	if (fx == 1. && fr == 0. || fx == 0. && fr == 1.){
		if (fx == 1. && fr == 0.){
			Nstr = 20;
			Nstep = 2*128;
			Dr = 0.02;

			xstream = (double*) malloc(Nstr * sizeof(double));
			ystream = (double*) malloc(Nstr * sizeof(double));

			for (istr = 0; istr < Nstr; istr++){
				xstream[istr] = 0.01;
				ystream[istr] = 0.10*(istr+1);
			}
			ystream[19] = 1.95;
			
		}
		else if (fx == 0. && fr == 1.){
			Nstr = 13;
			Nstep = 2*128;
			Dr = -0.02;

			xstream = (double*) malloc(Nstr * sizeof(double));
			ystream = (double*) malloc(Nstr * sizeof(double));

			for (istr = 0; istr < Nstr; istr++){
				xstream[istr] = -0.1 + 0.1*istr;
				ystream[istr] = 1.99;
			}
			xstream[0] = 0.01;
			xstream[1] = 0.05;
			
		}

		for (istr = 0; istr < Nstr; istr++){
			printf("%d\n",istr);

			// allocate memory
			xstr = (double*) malloc(1 * sizeof(double));
			ystr = (double*) malloc(1 * sizeof(double));
		
			// initialize streamline
			xstr[0] = xstream[istr];
			ystr[0] = ystream[istr];

			for (istep = 0; istep < Nstep; istep++){
				fprintf(pFile, "%.8f %.8f\n", xstr[istep], ystr[istep]);

				// update field point
				xvel = xstr[istep];
				yvel = ystr[istep];

				// reallocate memory
				xstr = (double*) realloc(xstr, (istep + 2) * sizeof(double));
				ystr = (double*) realloc(ystr, (istep + 2) * sizeof(double));
				
				// calculate the velocity induced by a ring of point forces
				gf_axR_vel(xvel, yvel, xring, yring, fx, fr, velx, vely);
				Dt = Dr/sqrt(velx*velx + vely*vely);

				//  update field point
				xvel = xstr[istep] + velx*Dt;
				yvel = ystr[istep] + vely*Dt;

				// calculate the induced velocity at the new field point
				gf_axR_vel(xvel, yvel, xring, yring, fx, fr, velx1, vely1);
				
				// update the streamline using a simple midpoint rule for integration
				xstr[istep + 1] = xstr[istep] + 0.5*(velx + velx1)*Dt;
				ystr[istep + 1] = ystr[istep] + 0.5*(vely + vely1)*Dt;

				// condition for breaking for loop
				if (xstr[istep + 1] > 2)
					break;
				if (xstr[istep + 1] < -2)
					break;
				if (ystr[istep + 1] > 2)
					break;
				if (ystr[istep + 1] < -2)
					break;
			}
			fprintf(pFile, "\n");

			free(xstr);
			free(ystr);
		}
	}

	fclose(pFile);

	
//	int k;
//	double complex z, s;
//	double complex Dk, dDkds;
//	z = 2.4 + 1.1*I;
//	
//	// get k
//	printf("k = ");
//	scanf("%d",&k);
//
//	s = 4. + 1.*I;
//
//	calcDk(k, s, Dk, dDkds);
//	printf("Dk = %.4f + %.4fi\n", creal(Dk), cimag(Dk));
//	printf("dDk/ds = %.4f + %.4fi\n", creal(dDkds), cimag(dDkds));	
//	
//	calcDk1Q(k, s, Dk, dDkds);
//	printf("Dk = %.4f + %.4fi\n", creal(Dk), cimag(Dk));
//	printf("dDk/ds = %.4f + %.4fi\n", creal(dDkds), cimag(dDkds));
	
//	int n = 1;
//	double an, bn, cn;
//	double complex xn[1], yn[1];
//	calcDkRoots(n, k, an, bn, cn, xn, yn);
//	printf("an = %.16f\nbn = %.16f\ncn = %.16f\n", an, bn, cn);
}

void testEllint(){
	double phi, k, n;
	double F, E, Pi, Kcomp, Ecomp, Picomp;
	
	phi = 1.2;
	k = 0.8;
	n = 0.3;

	F = ellintF(phi, k);
	E = ellintE(phi, k);
//	Pi = ellintPi(phi, k, n);
	Kcomp = ellintK(k);
	Ecomp = ellintE(k);
//	Picomp = ellintPi(k, n);

	printf("F = %.6f\n", F);
	printf("E = %.6f\n", E);
//	printf("Pi = %.6f\n", Pi);
	printf("Kcomp = %.6f\n", Kcomp);
	printf("Ecomp = %.6f\n", Ecomp);
//	printf("Picomp = %.6f\n", Picomp);
}

void testBesselNegativeOrder(){
	int n, nmin, nmax, *narray, nsize;
	double x, nu;
	double Jn, Yn, In, Kn;
	double *Jnarray, *Ynarray, *Inarray, *Knarray;

	x = 4.0;
	nsize = 10;
	
	// allocate memory
	narray = (int*) calloc(nsize,sizeof(int));
	Jnarray = (double*) calloc(nsize,sizeof(double));
	Ynarray = (double*) calloc(nsize,sizeof(double));
	Inarray = (double*) calloc(nsize,sizeof(double));
	Knarray = (double*) calloc(nsize,sizeof(double));
	
	n = -2;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	n = n + 1;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	n = n + 1;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	
	n = n + 1;
	Jn = besselJ(n,x);
	Yn = besselY(n,x);
	In = besselI(n,x);
	Kn = besselK(n,x);
	printf("n = %d, x = %.4f, Jn = %.4f, Yn = %.4f, In = %.4f, Kn = %.4f\n", n, x, Jn, Yn, In, Kn);

	for (int i = 0; i < nsize; i++){
		narray[i] = -2 + i; 
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
	
//	for (int i = 0; i < nsize; i++)
//		printf("%d  %.4f  %.4f  %.4f  %.4f\n", narray[i], Jnarray[i], Jnarray[i], Jnarray[i], Jnarray[i]);
	free(narray);
	free(Jnarray);
	free(Ynarray);
	free(Inarray);
	free(Knarray);
	
}

void testBesselComplexNegativeOrder(){
	double complex z;
	double complex Jn;
	double complex Yn;
	double complex In;
	double complex Kn;

	int n, nmin, nmax, nsize;
	double complex *Jnarray;
	double complex *Ynarray;
	double complex *Inarray;
	double complex *Knarray;

	n = -4;
	nmin = n;
	nmax = n + 7;
	nsize = nmax - nmin + 1;

	Jnarray = (double complex *) calloc(nsize, sizeof(double complex));
	Ynarray = (double complex *) calloc(nsize, sizeof(double complex));
	Inarray = (double complex *) calloc(nsize, sizeof(double complex));
	Knarray = (double complex *) calloc(nsize, sizeof(double complex));

	z = 2.4 + 1.1*I;
	printf("n = %d\n",n);
	printf("z = %.4f + %.4fi\n",creal(z), cimag(z));
	
	Jn = besselJ(n, z);
	printf("Jn = %.4f + %.4fi\n",creal(Jn), cimag(Jn));

	Yn = besselY(n, z);
	printf("Yn = %.4f + %.4fi\n",creal(Yn), cimag(Yn));
	
	In = besselI(n, z);
	printf("In = %.4f + %.4fi\n",creal(In), cimag(In));
	
	Kn = besselK(n, z);
	printf("Kn = %.4f + %.4fi\n",creal(Kn), cimag(Kn));
	
	
	besselJArray(nmin, nmax, z, Jnarray);
	besselYArray(nmin, nmax, z, Ynarray);
	besselIArray(nmin, nmax, z, Inarray);
	besselKArray(nmin, nmax, z, Knarray);
	
	printf("\n");
	printf("z = %.4f + %.4fi\n",creal(z), cimag(z));
	printf("n  Jn                Yn                 In                 Kn\n");
	for (int i = 0; i < nsize; i++){
		printf("%d", nmin + i);
		printf("  ");
		printf("%.4f + %.4fi",creal(Jnarray[i]), cimag(Jnarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Ynarray[i]), cimag(Ynarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Inarray[i]), cimag(Inarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Knarray[i]), cimag(Knarray[i]));
		printf("\n");
	//	
	//	printf("Jn = %.4f + %.4fi\n",creal(Jnarray[i]), cimag(Jnarray[i]));
	//	printf("Yn = %.4f + %.4fi\n",creal(Ynarray[i]), cimag(Ynarray[i]));
	//	printf("In = %.4f + %.4fi\n",creal(Inarray[i]), cimag(Inarray[i]));
	//	printf("Kn = %.4f + %.4fi\n",creal(Knarray[i]), cimag(Knarray[i]));
	}
	
	free(Jnarray);
	free(Ynarray);
	free(Inarray);
	free(Knarray);

}

void testBesselComplex(){
	double complex z;
	double complex Jn;
	double complex Yn;
	double complex In;
	double complex Kn;

	int n, nmin, nmax, nsize;
	double complex *Jnarray;
	double complex *Ynarray;
	double complex *Inarray;
	double complex *Knarray;

	n = 2;
	nmin = n;
	nmax = n + 4;
	nsize = nmax - nmin + 1;

	Jnarray = (double complex *) calloc(nsize, sizeof(double complex));
	Ynarray = (double complex *) calloc(nsize, sizeof(double complex));
	Inarray = (double complex *) calloc(nsize, sizeof(double complex));
	Knarray = (double complex *) calloc(nsize, sizeof(double complex));

	z = 2.4 + 1.1*I;
	printf("n = %d\n",n);
	printf("z = %.4f + %.4fi\n",creal(z), cimag(z));
	
	Jn = besselJ(n, z);
	printf("Jn = %.4f + %.4fi\n",creal(Jn), cimag(Jn));

	Yn = besselY(n, z);
	printf("Yn = %.4f + %.4fi\n",creal(Yn), cimag(Yn));
	
	In = besselI(n, z);
	printf("In = %.4f + %.4fi\n",creal(In), cimag(In));
	
	Kn = besselK(n, z);
	printf("Kn = %.4f + %.4fi\n",creal(Kn), cimag(Kn));
	
	
	besselJArray(nmin, nmax, z, Jnarray);
	besselYArray(nmin, nmax, z, Ynarray);
	besselIArray(nmin, nmax, z, Inarray);
	besselKArray(nmin, nmax, z, Knarray);
	
	printf("\n");
	printf("z = %.4f + %.4fi\n",creal(z), cimag(z));
	printf("n  Jn                Yn                 In                 Kn\n");
	for (int i = 0; i < nsize; i++){
		printf("%d", nmin + i);
		printf("  ");
		printf("%.4f + %.4fi",creal(Jnarray[i]), cimag(Jnarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Ynarray[i]), cimag(Ynarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Inarray[i]), cimag(Inarray[i]));
		printf("  ");
		printf("%.4f + %.4fi",creal(Knarray[i]), cimag(Knarray[i]));
		printf("\n");
	//	
	//	printf("Jn = %.4f + %.4fi\n",creal(Jnarray[i]), cimag(Jnarray[i]));
	//	printf("Yn = %.4f + %.4fi\n",creal(Ynarray[i]), cimag(Ynarray[i]));
	//	printf("In = %.4f + %.4fi\n",creal(Inarray[i]), cimag(Inarray[i]));
	//	printf("Kn = %.4f + %.4fi\n",creal(Knarray[i]), cimag(Knarray[i]));
	}
	
	free(Jnarray);
	free(Ynarray);
	free(Inarray);
	free(Knarray);
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
	Jn = besselJ(nu,x);
	Yn = besselY(nu,x);
	In = besselI(nu,x);
	Kn = besselK(nu,x);
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
