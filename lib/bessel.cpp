/* BESSEL FUNCTIONS
 *  Evaluate Bessel functions and modified Bessel functions.
 *
 * REFERENCES
 *  Abramowitz and Stegun, Dover Publications (1964)
 *
 * PARAMETERS
 *  x	[input]		field point
 *  n,m	[input]		integer order
 *  nu	[input]		fractional order
 *  J	[output]        Bessel function of the first kind
 *  Y	[output]        Bessel function of the second kind
 *  I	[output]        modified Bessel function of the first kind
 *  K	[output]        modified Bessel function of the second kind
 *
 * NOTES
 *  GSL provides a basic implementation of Bessel functions, but
 *  does not support Bessel functions with complex arguments. The
 *  latter is handled by SLATEC subroutines.
 *
 *  The current implementation cannot handle negative non-integer
 *  orders. For negative integer orders, the following rules for
 *  Bessel functions are used:
 *   Jn(x) = -J(-n)(x)	if n is odd	[ditto for Yn(x)]
 *   Jn(x) = J(-n)(x)	if n is even	[ditto for Yn(x)]
 *   In(x) = I(-n)(x)	for all n	[ditto for Kn(x)]
 */

/* HEADER FILES */
#include "bessel.h"

/* IMPLEMENTATIONS */
// Bessel functions of the first kind
double besselJ(int n, double x){
	double Jn;
	
	if (n < 0){
		if (n % 2)
			Jn = -gsl_sf_bessel_Jn(-n, x);
		else
			Jn = gsl_sf_bessel_Jn(-n, x);
	}
	else if (n == 0)
		Jn = gsl_sf_bessel_J0(x);
	else if (n == 1)
		Jn = gsl_sf_bessel_J1(x);
	else
		Jn = gsl_sf_bessel_Jn(n, x);
	
	return(Jn);
}

double besselJ(double nu, double x){
	double Jnu;
	
	Jnu = gsl_sf_bessel_Jnu(nu, x);

	return(Jnu);
}

double complex besselJ(int m, double complex z){
	double nu = double(m);
	double complex Jm;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	if (m < 0){
		if (m % 2) {
			nu = -nu;
			zbesj_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
			Jm = -(cyr[0] + I*cyi[0]);
		}
		else {
			nu = -nu;
			zbesj_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
			Jm = cyr[0] + I*cyi[0];
		}
	}
	else {
		zbesj_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
		Jm = cyr[0] + I*cyi[0];
	}
	
	if (ierr != 0)
		printf("Error! Zbesj did not return successfully.\n");

	return(Jm);
}

double complex besselJ(double nu, double complex z){
	double complex Jnu;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	zbesj_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);

	if (ierr != 0)
		printf("Error! Zbesj did not return successfully.\n");
	
	Jnu = cyr[0] + I*cyi[0];

	return(Jnu);
}

void besselJArray(int nmin, int nmax, double x, double *Jn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int info, i;
	int nsize = nmax - nmin + 1;

	if (nmin < 0 && nmax <= 0) {
		double* tmp;
		tmp = (double*) calloc(nsize, sizeof(double));

		info = gsl_sf_bessel_Jn_array(-nmax, -nmin, x, tmp);

		if ((-nmin) % 2) {
			for (i = 0; i < nsize; i++){
				Jn[i] = pow(-1, i + 1)*tmp[nsize - i - 1];
			}
		}
		else {
			for (i = 0; i < nsize; i++){
				Jn[i] = pow(-1, i)*tmp[nsize - i - 1];
			}
		}

		free(tmp);
	}
	else if (nmin < 0 && nmax > 0) {
		int info2;
		double *tmp, *tmp2;
		tmp = (double*) calloc(-nmin, sizeof(double));
		tmp2 = (double*) calloc(nmax + 1, sizeof(double));

		info = gsl_sf_bessel_Jn_array(1, -nmin, x, tmp);
		info2 = gsl_sf_bessel_Jn_array(0, nmax, x, tmp2);
		
		if ((-nmin) % 2) {
			for (i = 0; i < (-nmin); i++){
				Jn[i] = pow(-1, i + 1)*tmp[-nmin - i - 1];
			}
		}
		else {
			for (i = 0; i < (-nmin); i++){
				Jn[i] = pow(-1, i)*tmp[-nmin - i - 1];
			}
		}
		
		for (i = 0; i < (nmax + 1); i++){
			Jn[-nmin + i] = tmp2[i];
		}

		free(tmp);
		free(tmp2);
	}
	else
		info = gsl_sf_bessel_Jn_array(nmin, nmax, x, Jn);
}

void besselJArray(int nmin, int nmax, double complex z, double complex *Jn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int i;
	int nsize = nmax - nmin + 1;
	double nu;
	int kode = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double *cyr, *cyi;

	cyr = (double*) calloc(nsize, sizeof(double));
	cyi = (double*) calloc(nsize, sizeof(double));
	
	if (nmin < 0 && nmax <= 0) {
		nu = double(-nmax);

		zbesj_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		
		if ((-nmin) % 2) {
			for (i = 0; i < nsize; i++){
				Jn[i] = pow(-1, i + 1)*(cyr[nsize - i - 1] + I*cyi[nsize - i - 1]);
			}
		}
		else {
			for (i = 0; i < nsize; i++){
				Jn[i] = pow(-1, i)*(cyr[nsize - i - 1] + I*cyi[nsize - i - 1]);
			}
		}
	}
	else if (nmin < 0 && nmax > 0) {
		int nneg = -nmin;
		int npos = nmax + 1;
	
		nu = 1.;
		zbesj_(&zr, &zi, &nu, &kode, &nneg, cyr, cyi, &nz, &ierr);

		if (nneg % 2) {
			for (i = 0; i < nneg; i++){
				Jn[i] = pow(-1, i + 1)*(cyr[nneg - i - 1] + I*cyi[nneg - i - 1]);
			}
		}
		else {
			for (i = 0; i < nneg; i++){
				Jn[i] = pow(-1, i)*(cyr[nneg - i - 1] + I*cyi[nneg - i - 1]);
			}
		}

		nu = 0.;
		zbesj_(&zr, &zi, &nu, &kode, &npos, cyr, cyi, &nz, &ierr);
		for (i = 0; i < npos; i++){
			Jn[nneg + i] = cyr[i] + I*cyi[i];
		}

	}
	else {
		nu = double(nmin);
		zbesj_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			Jn[i] = cyr[i] + I*cyi[i];
		}
	}

	if (ierr != 0)
		printf("Error! Zbesj did not return successfully.\n");

	free(cyr);
	free(cyi);
}

// Bessel functions of the second kind
double besselY(int n, double x){
	double Yn;
	
	if (n < 0){
		if (n % 2)
			Yn = -gsl_sf_bessel_Yn(-n, x);
		else
			Yn = gsl_sf_bessel_Yn(-n, x);
	}
	else if (n == 0)
		Yn = gsl_sf_bessel_Y0(x);
	else if (n == 1)
		Yn = gsl_sf_bessel_Y1(x);
	else
		Yn = gsl_sf_bessel_Yn(n, x);
	
	return(Yn);
}

double besselY(double nu, double x){
	double Ynu;
	
	Ynu = gsl_sf_bessel_Ynu(nu, x);

	return(Ynu);
}

double complex besselY(int m, double complex z){
	double nu = double(m);
	double complex Ym;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	double cwrkr[1], cwrki[1];
	
	if (m < 0){
		if (m % 2) {
			nu = -nu;
			zbesy_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
			Ym = -(cyr[0] + I*cyi[0]);
		}
		else {
			nu = -nu;
			zbesy_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
			Ym = cyr[0] + I*cyi[0];
		}
	}
	else {
		zbesy_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
		Ym = cyr[0] + I*cyi[0];
	}

	if (ierr != 0)
		printf("Error! Zbesy did not return successfully.\n");

	return(Ym);
}

double complex besselY(double nu, double complex z){
	double complex Ynu;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	double cwrkr[1], cwrki[1];
	
	zbesy_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
	
	if (ierr != 0)
		printf("Error! Zbesy did not return successfully.\n");

	Ynu = cyr[0] + I*cyi[0];

	return(Ynu);
}

void besselYArray(int nmin, int nmax, double x, double *Yn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int info, i;
	int nsize = nmax - nmin + 1;

	if (nmin < 0 && nmax <= 0) {
		double* tmp;
		tmp = (double*) calloc(nsize, sizeof(double));

		info = gsl_sf_bessel_Yn_array(-nmax, -nmin, x, tmp);

		if ((-nmin) % 2) {
			for (i = 0; i < nsize; i++){
				Yn[i] = pow(-1, i + 1)*tmp[nsize - i - 1];
			}
		}
		else {
			for (i = 0; i < nsize; i++){
				Yn[i] = pow(-1, i)*tmp[nsize - i - 1];
			}
		}

		free(tmp);
	}
	else if (nmin < 0 && nmax > 0) {
		int info2;
		double *tmp, *tmp2;
		tmp = (double*) calloc(-nmin, sizeof(double));
		tmp2 = (double*) calloc(nmax + 1, sizeof(double));

		info = gsl_sf_bessel_Yn_array(1, -nmin, x, tmp);
		info2 = gsl_sf_bessel_Yn_array(0, nmax, x, tmp2);
		
		if ((-nmin) % 2) {
			for (i = 0; i < (-nmin); i++){
				Yn[i] = pow(-1, i + 1)*tmp[-nmin - i - 1];
			}
		}
		else {
			for (i = 0; i < (-nmin); i++){
				Yn[i] = pow(-1, i)*tmp[-nmin - i - 1];
			}
		}
		
		for (i = 0; i < (nmax + 1); i++){
			Yn[-nmin + i] = tmp2[i];
		}

		free(tmp);
		free(tmp2);
	}
	else
		info = gsl_sf_bessel_Yn_array(nmin, nmax, x, Yn);
}

void besselYArray(int nmin, int nmax, double complex z, double complex *Yn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int i;
	int nsize = nmax - nmin + 1;
	double nu;
	int kode = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double *cyr, *cyi;
	double *cwrkr, *cwrki;

	cyr = (double*) calloc(nsize, sizeof(double));
	cyi = (double*) calloc(nsize, sizeof(double));
	cwrkr = (double*) calloc(nsize,sizeof(double));
	cwrki = (double*) calloc(nsize,sizeof(double));
	
	if (nmin < 0 && nmax <= 0) {
		nu = double(-nmax);

		zbesy_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
		
		if ((-nmin) % 2) {
			for (i = 0; i < nsize; i++){
				Yn[i] = pow(-1, i + 1)*(cyr[nsize - i - 1] + I*cyi[nsize - i - 1]);
			}
		}
		else {
			for (i = 0; i < nsize; i++){
				Yn[i] = pow(-1, i)*(cyr[nsize - i - 1] + I*cyi[nsize - i - 1]);
			}
		}
	}
	else if (nmin < 0 && nmax > 0) {
		int nneg = -nmin;
		int npos = nmax + 1;
	
		nu = 1.;
		zbesy_(&zr, &zi, &nu, &kode, &nneg, cyr, cyi, &nz, cwrkr, cwrki, &ierr);

		if (nneg % 2) {
			for (i = 0; i < nneg; i++){
				Yn[i] = pow(-1, i + 1)*(cyr[nneg - i - 1] + I*cyi[nneg - i - 1]);
			}
		}
		else {
			for (i = 0; i < nneg; i++){
				Yn[i] = pow(-1, i)*(cyr[nneg - i - 1] + I*cyi[nneg - i - 1]);
			}
		}

		nu = 0.;
		zbesy_(&zr, &zi, &nu, &kode, &npos, cyr, cyi, &nz, cwrkr, cwrki, &ierr);
		for (i = 0; i < npos; i++){
			Yn[nneg + i] = cyr[i] + I*cyi[i];
		}

	}
	else {
		nu = double(nmin);
		zbesj_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			Yn[i] = cyr[i] + I*cyi[i];
		}
	}
	
	if (ierr != 0)
		printf("Error! Zbesy did not return successfully.\n");

	free(cyr);
	free(cyi);
	free(cwrkr);
	free(cwrki);
}

// Modified Bessel functions of the first kind
double besselI(int n, double x){
	double In;
	
	if (n < 0)
		In = gsl_sf_bessel_In(-n, x);
	else if (n == 0)
		In = gsl_sf_bessel_I0(x);
	else if (n == 1)
		In = gsl_sf_bessel_I1(x);
	else
		In = gsl_sf_bessel_In(n, x);
	
	return(In);
}

double besselI(double nu, double x){
	double Inu;
	
	Inu = gsl_sf_bessel_Inu(nu, x);

	return(Inu);
}

double complex besselI(int m, double complex z){
	double nu = double(m);
	double complex Im;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	if (m < 0){
		nu = -nu;
		zbesi_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
		Im = cyr[0] + I*cyi[0];
	}
	else {
		zbesi_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
		Im = cyr[0] + I*cyi[0];
	}

	if (ierr != 0)
		printf("Error! Zbesi did not return successfully.\n");

	return(Im);
}

double complex besselI(double nu, double complex z){
	double complex Inu;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	zbesi_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
	
	if (ierr != 0)
		printf("Error! Zbesi did not return successfully.\n");

	Inu = cyr[0] + I*cyi[0];

	return(Inu);
}

void besselIArray(int nmin, int nmax, double x, double *In){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int info, i;
	int nsize = nmax - nmin + 1;

	if (nmin < 0 && nmax <= 0) {
		double* tmp;
		tmp = (double*) calloc(nsize, sizeof(double));

		info = gsl_sf_bessel_In_array(-nmax, -nmin, x, tmp);

		for (i = 0; i < nsize; i++){
			In[i] = tmp[nsize - i - 1];
		}

		free(tmp);
	}
	else if (nmin < 0 && nmax > 0) {
		int info2;
		double *tmp, *tmp2;
		tmp = (double*) calloc(-nmin, sizeof(double));
		tmp2 = (double*) calloc(nmax + 1, sizeof(double));

		info = gsl_sf_bessel_In_array(1, -nmin, x, tmp);
		info2 = gsl_sf_bessel_In_array(0, nmax, x, tmp2);
		
		for (i = 0; i < (-nmin); i++){
			In[i] = tmp[-nmin - i - 1];
		}
		
		for (i = 0; i < (nmax + 1); i++){
			In[-nmin + i] = tmp2[i];
		}

		free(tmp);
		free(tmp2);
	}
	else
		info = gsl_sf_bessel_In_array(nmin, nmax, x, In);
}

void besselIArray(int nmin, int nmax, double complex z, double complex *In){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int i;
	int nsize = nmax - nmin + 1;
	double nu = double(nmin);
	int kode = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double *cyr, *cyi;

	cyr = (double*) calloc(nsize, sizeof(double));
	cyi = (double*) calloc(nsize, sizeof(double));
	
	if (nmin < 0 && nmax <= 0) {
		nu = double(-nmax);

		zbesi_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			In[i] = cyr[nsize - i - 1] + I*cyi[nsize - i - 1];
		}
	}
	else if (nmin < 0 && nmax > 0) {
		int nneg = -nmin;
		int npos = nmax + 1;
	
		nu = 1.;
		zbesi_(&zr, &zi, &nu, &kode, &nneg, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nneg; i++){
			In[i] = cyr[nneg - i - 1] + I*cyi[nneg - i - 1];
		}

		nu = 0.;
		zbesi_(&zr, &zi, &nu, &kode, &npos, cyr, cyi, &nz, &ierr);
		for (i = 0; i < npos; i++){
			In[nneg + i] = cyr[i] + I*cyi[i];
		}

	}
	else {
		nu = double(nmin);
		zbesi_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			In[i] = cyr[i] + I*cyi[i];
		}
	}

	if (ierr != 0)
		printf("Error! Zbesi did not return successfully.\n");
	
	free(cyr);
	free(cyi);
}

// Modified bessel functions of the second kind
double besselK(int n, double x){
	double Kn;
	
	if (n < 0)
		Kn = gsl_sf_bessel_Kn(-n, x);
	if (n == 0)
		Kn = gsl_sf_bessel_K0(x);
	else if (n == 1)
		Kn = gsl_sf_bessel_K1(x);
	else
		Kn = gsl_sf_bessel_Kn(n, x);
	
	return(Kn);
}

double besselK(double nu, double x){
	double Knu;
	
	Knu = gsl_sf_bessel_Knu(nu, x);

	return(Knu);
}

double complex besselK(int m, double complex z){
	double nu = double(m);
	double complex Km;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	if (m < 0){
		nu = -nu;
		zbesk_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
		Km = cyr[0] + I*cyi[0];
	}
	else {
		zbesk_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
		Km = cyr[0] + I*cyi[0];
	}
	
	if (ierr != 0)
		printf("Error! Zbesk did not return successfully.\n");

	return(Km);
}

double complex besselK(double nu, double complex z){
	double complex Knu;
	int kode = 1;
	int n = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double cyr[1], cyi[1];
	
	zbesk_(&zr, &zi, &nu, &kode, &n, cyr, cyi, &nz, &ierr);
	
	if (ierr != 0)
		printf("Error! Zbesk did not return successfully.\n");

	Knu = cyr[0] + I*cyi[0];

	return(Knu);
}
void besselKArray(int nmin, int nmax, double x, double *Kn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int info, i;
	int nsize = nmax - nmin + 1;

	if (nmin < 0 && nmax <= 0) {
		double* tmp;
		tmp = (double*) calloc(nsize, sizeof(double));

		info = gsl_sf_bessel_Kn_array(-nmax, -nmin, x, tmp);

		for (i = 0; i < nsize; i++){
			Kn[i] = tmp[nsize - i - 1];
		}

		free(tmp);
	}
	else if (nmin < 0 && nmax > 0) {
		int info2;
		double *tmp, *tmp2;
		tmp = (double*) calloc(-nmin, sizeof(double));
		tmp2 = (double*) calloc(nmax + 1, sizeof(double));

		info = gsl_sf_bessel_Kn_array(1, -nmin, x, tmp);
		info2 = gsl_sf_bessel_Kn_array(0, nmax, x, tmp2);
		
		for (i = 0; i < (-nmin); i++){
			Kn[i] = tmp[-nmin - i - 1];
		}
		
		for (i = 0; i < (nmax + 1); i++){
			Kn[-nmin + i] = tmp2[i];
		}

		free(tmp);
		free(tmp2);
	}
	else
		info = gsl_sf_bessel_Kn_array(nmin, nmax, x, Kn);
}

void besselKArray(int nmin, int nmax, double complex z, double complex *Kn){
	if (nmin > nmax)
		printf("Error! nmin > nmax.\n");
	
	int i;
	int nsize = nmax - nmin + 1;
	double nu = double(nmin);
	int kode = 1;
	double zr = creal(z);
	double zi = cimag(z);
	int nz, ierr;
	double *cyr, *cyi;

	cyr = (double*) calloc(nsize, sizeof(double));
	cyi = (double*) calloc(nsize, sizeof(double));
	
	if (nmin < 0 && nmax <= 0) {
		nu = double(-nmax);

		zbesk_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			Kn[i] = cyr[nsize - i - 1] + I*cyi[nsize - i - 1];
		}
	}
	else if (nmin < 0 && nmax > 0) {
		int nneg = -nmin;
		int npos = nmax + 1;
	
		nu = 1.;
		zbesk_(&zr, &zi, &nu, &kode, &nneg, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nneg; i++){
			Kn[i] = cyr[nneg - i - 1] + I*cyi[nneg - i - 1];
		}

		nu = 0.;
		zbesk_(&zr, &zi, &nu, &kode, &npos, cyr, cyi, &nz, &ierr);
		for (i = 0; i < npos; i++){
			Kn[nneg + i] = cyr[i] + I*cyi[i];
		}

	}
	else {
		nu = double(nmin);
		zbesk_(&zr, &zi, &nu, &kode, &nsize, cyr, cyi, &nz, &ierr);
		for (i = 0; i < nsize; i++){
			Kn[i] = cyr[i] + I*cyi[i];
		}
	}

	if (ierr != 0)
		printf("Error! Zbesk did not return successfully.\n");
	
	free(cyr);
	free(cyi);
}
