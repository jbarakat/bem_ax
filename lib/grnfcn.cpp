/* GREEN'S FUNCTIONS
 *  Evaluate axisymmetric Green's functions for Stokes flow.
 *
 * REFERENCES
 *  Pozrikidis, Cambridge University Press (1992) [pp. 51, 89-91]
 *  Tozeren, Inter. J. Num. Meth. Fluids 4, 159-170 (1984) 
 *  
 * PARAMETERS
 *  x,r   [input]   field point
 *  x0,r0 [input]   source point
 *  rc    [input]   cylindrical tube radius
 *  f     [input]   force density
 *  u     [output]  velocity
 *  M 		[output]  Green's function
 */

/* HEADER FILES */
#include "grnfcn.h"


/* NOTE: SHOULD ALSO PROGRAM THE GREEN'S FUNCTIONS FOR FLOW BOUNDED
 * INTERNALLY BY A SOLID SPHERE, REWRITTEN IN AXISYMMETRIC COORDINATES
 * IN THE SAME WAY AS THE FREE-SPACE GREEN'S FUNCTION. THIS WILL
 * REQUIRE QUITE A BIT OF ALGEBRA, SO DO THIS LATER
 * see p. 51, 99 in Pozrikidis
 *
 * NOTE: I have not removed the singularity at r = 0 in the free-space
 * Green's function w/the stresslet... need to address this at some
 * point.
 *
 */



/* IMPLEMENTATIONS */
/* Green's function M evaluated at (x,r) due to a ring of point forces
 * located at (x0,r0) in an unbound domain (free space).
 */
void gf_axR(double x, double r, double x0, double r0,
            double &Mxx, double &Mxr, double &Mrx, double &Mrr){
	// declare variables
	double K  , E  ;
	double x2 , x02;
	double r2 , r02;
	double X  , X2 ;
	double R  , R2 ;
	double k  , k2 , k4, kp2;
	double fc , fc3;
	double fcM;
  double I10, I11, I30, I31, I32;

	x2  = x*x;
	x02 = x0*x0;
	r2  = r*r;
	r02 = r0*r0;
	X   = x - x0; 
	X2  = X*X;
	R2  = X2 + (r - r0)*(r - r0);
	R   = pow(R2, 0.5);
	k2  = 4*r*r0/(X2 + (r + r0)*(r + r0));
	k   = pow(k2, 0.5);
	k4  = k2*k2;
	kp2 = 1 - k2;
	
	if (r < 1e-8){
		Mxx =  2*M_PI*(r0/R)*(1 + X2/R2);
		Mxr = -2*M_PI*(X*r02)/(R*R2); 
		Mrx =  0;
		Mrr =  0;
	}
	else {
	  fc  = 0.5*k/pow(r*r0, 0.5);
	  fc3 = fc*fc*fc;
		fcM = r0;
		
		fc  *= 4;
		fc3 *= (4/kp2);
	
		// evaluate complete elliptic integrals
		K = ellintK(k);
		E = ellintE(k);
		
	  // reduce integrals in terms of complete elliptic integrals
	  I10 = fc  *                          K;
		I11 = fc  * ((  2 -    k2        ) * K
						    -   2                  * E) / k2;
		I30 = fc3 *                          E; 
		I31 = fc3 * ((- 2 +  2*k2        ) * K  
		            +(  2 -    k2        ) * E) / k2;
		I32 = fc3 * ((- 8 + 12*k2 - 4*k4 ) * K  
		            +(  8 -  8*k2 +   k4 ) * E) / k4; 
		
		if (r < 0.001)
			printf("%.16f\n", I31);
	
	  // calculate components of the Stokeslet M
	  Mxx =  fcM*  (     I10 + X2        *I30);
	  Mxr = -fcM*X*(r0  *I30 - r         *I31);
	  Mrx = -fcM*X*(r0  *I31 - r         *I30);
	  Mrr =  fcM*  (     I11 + (r02 + r2)*I31 - r*r0*(I30 + I32));
	}

}

/* Green's functions M and Q evaluated at (x,r) due to a ring of point
 * forces located at (x0,r0) in an unbound domain (free space).
 */

/* NOTE: I have not removed the singularity at r = 0 in this function yet! */
/* FIX THIS LATER (SEE FUNCTION ABOVE) */

void gf_axR(double x, double r, double x0, double r0, 
            double &Mxx, double &Mxr, double &Mrx, double &Mrr,
            double &Qxxx, double &Qxxr, double &Qxrx, double &Qxrr,
            double &Qrxx, double &Qrxr, double &Qrrx, double &Qrrr){
  // declare variables
  double K, E;
  double x2 = x*x;
  double x02 = x0*x0;
  double r2 = r*r;
  double r02 = r0*r0;
  double X = x - x0; 
  double X2 = X*X;
  double R2 = X2 + (r - r0)*(r - r0);
  double R = pow(R2, 0.5);
  double k2 = 4*r*r0/(X2 + (r + r0)*(r + r0));
  double k = pow(k2, 0.5);
  double k4 = k2*k2;
  double k6 = k4*k2;
  double k8 = k6*k2;
	double kp2 = 1 - k2;
	double kp4 = kp2*kp2;
  double fc = 0.5*k/pow(r*r0, 0.5);
  double fc3 = fc*fc*fc;
	double fc5 = fc3*fc*fc;
	double fcM = r0;
	double fcQ = 6*r0;
  double I10, I11, I30, I31, I32, I50, I51, I52, I53;

	fc  *= 4;
	fc3 *= (4/kp2);
	fc5 *= (4/(3*kp4));

  // evaluate complete elliptic integrals
  K = ellintK(k);
  E = ellintE(k);

  // reduce integrals in terms of complete elliptic integrals
  I10 = fc  *                                           K;
	I11 = fc  * ((  2 -     k2                        ) * K
					    -   2                                   * E) / k2;
	I30 = fc3 *                                           E; 
	I31 = fc3 * ((- 2 +   2*k2                        ) * K  
	            +(  2 -     k2                        ) * E) / k2;
	I32 = fc3 * ((- 8 +  12*k2 -   4*k4               ) * K  
	            +(  8 -   8*k2 +     k4               ) * E) / k4; 
	I50 = fc5 * ((- 1 +     k2                        ) * K
	            +(  4 -   2*k2                        ) * E);
	I51 = fc5 * ((- 2 +   3*k2 -     k4               ) * K
	            +(  2 -   2*k2 +   2*k4               ) * E) / k2;
	I52 = fc5 * ((  8 -  16*k2 +   7*k4 +    k6       ) * K
	            +(- 8 +  12*k2          -  2*k6       ) * E) / k4; 
	I53 = fc5 * (( 64 - 160*k2 + 126*k4 - 29*k6 -   k8) * K
	            +(-64 + 128*k2 -  66*k4 +  2*k6 + 2*k8) * E) / k6;

  // calculate components of the Stokeslet M
  Mxx =  fcM*  (     I10 + X2        *I30);
  Mxr = -fcM*X*(r0  *I30 - r         *I31);
  Mrx = -fcM*X*(r0  *I31 - r         *I30);
  Mrr =  fcM*  (     I11 + (r02 + r2)*I31 - r*r0*(I30 + I32));

	// calculate components of the stresslet Q
	Qxxx =  fcQ*X *  X2       *I50;
	Qxxr = -fcQ*X *( r0       *I50 - r    * I51);
	Qxrr =  fcQ*X *( r2       *I52 + r02  * I50 - 2*r0*r*I51);
	Qrxx = -fcQ*X2*( r0       *I51 - r    * I50);
	Qrxr =  fcQ*X *((r02 + r2)*I51 - r0 *r*(I50 +   I52));
	Qrrr = -fcQ*   ( r02*r0   *I51 - r02*r*(I50 + 2*I52) + r0*r2*(I53 + 2*I51) - r*r2*I52);
	Qxrx = Qxxr;
	Qrrx = Qrxr;

}


/* velocity u evaluated at (x,r) due to a ring of point forces 
 * of strength f located at (x0,r0) in an unbound domain (free space).
 */
void gf_axR_vel(double x, double r, double x0, double r0,
                double fx, double fr, double &ux, double &ur){
	// declare variables
	double Mxx, Mxr, Mrx, Mrr;
	
	// calculate components of the Green's function
	gf_axR(x, r, x0, r0, Mxx, Mxr, Mrx, Mrr);

	// calculate velocity components
	ux = (Mxx*fx + Mxr*fr)/(8*M_PI);
	ur = (Mrx*fx + Mrr*fr)/(8*M_PI);

}

/* Complementary Green's function MC evaluated at (x,r) due to a ring of 
 * point forces located at (x0,r0), bounded externally by a cylindrical
 * tube of radius rc.
 */
void gf_axC(double x, double r, double x0, double r0, double rc,
            double &MCxx, double &MCxr, double &MCrx, double &MCrr){
	/* set maximum number of refinement stages
	 * and tolerance for Fourier integrals */
	const int MAXIT = 20;
	const double TOL = 0.00001;

	// evaluate Green's function
	gf_axC(x, r, x0, r0, rc, MCxx, MCxr, MCrx, MCrr, MAXIT, TOL);
}

void gf_axC(double x, double r, double x0, double r0, double rc,
            double &MCxx, double &MCxr, double &MCrx, double &MCrr,
						const int MAXIT, const double TOL){
	// declare variables
	int i, j, k, n;
	double mCxx, mCxr, mCrx, mCrr;
	double dmCxx, dmCxr, dmCrx, dmCrr;
	double fc, dev;
	
	int na, nt;
	double s, ds;
	double smin = 0.;
	double smax = 1.;
	double Ds = smax - smin;
	
	// initialize
	mCxx = 0.;
	mCxr = 0.;
	mCrx = 0.;
	mCrr = 0.;

	dmCxx = 0.;
	dmCxr = 0.;
	dmCrx = 0.;
	dmCrr = 0.;
	
	s = 0.5*Ds;
	ds = Ds;
	na = 1;  // number of points added
	nt = 1;  // total number of points

	/* evaluate Fourier integrals using the extended midpoint rule,
	 * tripling the number of integration points at each stage of
	 * refinement (incremented by n).
	 */
	for (n = 0; n < MAXIT; n++){
		if (n == 0){
			// evaluate kernels at midpoint of the domain
			gf_axC_ker(x, r, x0, r0, rc, s, mCxx, mCxr, mCrx, mCrr);

			
			// calculate change in kernels
			dmCxx -= mCxx*ds; dmCxx = fabs(dmCxx);
			dmCxr -= mCxr*ds; dmCxr = fabs(dmCxr);
			dmCrx -= mCrx*ds; dmCrx = fabs(dmCrx);
			dmCrr -= mCrr*ds; dmCrr = fabs(dmCrr);

			// prepare for next stage of refinement
			na = 2;
			nt = 3;
		}
		else {
			// add 2 * 3^(n-1) additional points
			if (n != 1){
				na *= 3;
				nt *= 3;
			}
			
			// store kernel before adding new points
			dmCxx = mCxx*ds;
			dmCxr = mCxr*ds;
			dmCrx = mCrx*ds;
			dmCrr = mCrr*ds;

			// refine grid spacing
			ds /= 3.0;
			
			// evaluate kernels at additional grid points
			for (k = 0; k < nt; k++){
				if ((2*k + 1) % 3 != 0){ /* avoid double counting
																	  previous grid points */
					s = smin + (double(k) + 0.5)*ds;
					gf_axC_ker(x, r, x0, r0, rc, s, mCxx, mCxr, mCrx, mCrr);
				}
			}

			// calculate change in kernels
			dmCxx -= mCxx*ds; dmCxx = fabs(dmCxx);
			dmCxr -= mCxr*ds; dmCxr = fabs(dmCxr);
			dmCrx -= mCrx*ds; dmCrx = fabs(dmCrx);
			dmCrr -= mCrr*ds; dmCrr = fabs(dmCrr);

		}

		// calculate maximum deviation
		dev = fmax(dmCxx, dmCxr);
		dev = fmax(dev, dmCrx);
		dev = fmax(dev, dmCrr);

		// break loop when integral converges
		if (dev < TOL)
			break;
	}
	
	// diagnose level of refinement [uncomment when required]
//	printf("%d stages of refinement and %d total points\n", n, nt);
	
	// calculate components of the complementary Green's function
	MCxx = r0*mCxx*ds - 4*M_PI*r0/pow((x - x0)*(x - x0) + (2*rc - r - r0)*(2*rc - r - r0), 0.5);
	MCxr = r0*mCxr*ds;
	MCrx = r0*mCrx*ds;
	MCrr = r0*mCrr*ds;
	
}

/* Total Green's function MT = MR + MC evaluated at (x,r) due to a ring of 
 * point forces located at (x0,r0), bounded externally by a cylindrical
 * tube of radius rc.
 */
void gf_axT(double x, double r, double x0, double r0, double rc,
            double &MTxx, double &MTxr, double &MTrx, double &MTrr){
	/* set maximum number of refinement stages
	 * and tolerance for Fourier integrals */
	const int MAXIT = 20;
	const double TOL = 0.00001;

	// evaluate Green's function
	gf_axT(x, r, x0, r0, rc, MTxx, MTxr, MTrx, MTrr, MAXIT, TOL);
}

void gf_axT(double x, double r, double x0, double r0, double rc,
            double &MTxx, double &MTxr, double &MTrx, double &MTrr,
						const int MAXIT, const double TOL){
	// declare variables
	double MRxx, MRxr, MRrx, MRrr;
	double MCxx, MCxr, MCrx, MCrr;
	
	// calculate components of the free-space Green's function
	gf_axR(x, r, x0, r0, MRxx, MRxr, MRrx, MRrr);

	// calculate components of the complementary Green's function
	gf_axC(x, r, x0, r0, rc, MCxx, MCxr, MCrx, MCrr);
	
	// calculate components of the total Green's function
	MTxx = MRxx + MCxx;
	MTxr = MRxr + MCxr;
	MTrx = MRrx + MCrx;
	MTrr = MRrr + MCrr;

}

/* Evaluate the kernel of the complementary Green's function, MC.
 *
 * NOTE: The lines commented out comprise the code for evaluating the
 * kernel of the free-space Green's function, MR, in terms of Fourier
 * integrals. This method is not effective if r < r0, for which the
 * components of B diverge like exp[(r0 - r)*t] and the integrals are
 * made improper. Instead, the closed-form expression for MR in terms
 * of complete elliptic integrals is used (see gf_axR above).
 */
void gf_axC_ker(double x, double r, double x0, double r0, double rc, double s,
            double &mCxx, double &mCxr, double &mCrx, double &mCrr){
	// declare variables
	double t;
	double X = x - x0;
	double Xt, cosXt, sinXt;
	double omeg, omeg0, omegc, omegn;
	double I0, I1, I00, I10, I0c, I1c, I0n, I1n;
	double K0, K1, K00, K10, K0c, K1c, K0n, K1n;
	double Axx, Axr, Arx, Arr;
//	double Bxx, Bxr, Brx, Brr;
	double Bxxc, Bxrc, Brxc, Brrc;
	double Lxx, Lxr, Lrx, Lrr;
	double Fxx, Fxr, Frx, Frr;
	double detL, fc;

	// change of integration coordinate
	t = -gsl_sf_log(s);

	
	// evaluate trigonometric functions
	Xt = X*t;
	cosXt = gsl_sf_cos(Xt);
	sinXt = gsl_sf_sin(Xt);
	
	// scale radial coordinate
	omeg = t*r;
	omeg0 = t*r0;
	omegc = t*rc;
	omegn = 2*omegc - omeg - omeg0;

	// evaluate modified Bessel functions
	I0  = besselI(0, omeg);
	I1  = besselI(1, omeg);
	I00 = besselI(0, omeg0);
	I10 = besselI(1, omeg0);
	I0c = besselI(0, omegc);
	I1c = besselI(1, omegc);
	I0n = besselI(0, omegn);
	I1n = besselI(1, omegn);
	
	K0  = besselK(0, omeg);
	K1  = besselK(1, omeg);
	K00 = besselK(0, omeg0);
	K10 = besselK(1, omeg0);
	K0c = besselK(0, omegc);
	K1c = besselK(1, omegc);
	K0n = besselK(0, omegn);
	K1n = besselK(1, omegn);

	// calculate components of L matrix in La = 4b
	Lxx =     t*I0c;
	Lxr = omegc*I1c + 2*I0c;
	Lrx =     t*I1c;
	Lrr = omegc*I0c;
	detL = Lxx*Lrr - Lxr*Lrx;
	if (detL == 0)
		printf("Error: matrix L is singular.\n");
	fc = 4/detL; 

	// calculate components of B matrix
	Bxxc = (-2*K0c +  omegc*K1c)*I00 - omeg0*K0c*I10;
	Bxrc =           -omegc*K0c *I00 + omeg0*K1c*I10;
	Brxc =           -omegc*K1c *I10 + omeg0*K0c*I00;
	Brrc = ( 2*K1c +  omegc*K0c)*I10 - omeg0*K1c*I00;
	
	// calculate components of A matrix from inverting La = 4b
	Axx = fc*( Lrr*Bxxc - Lxr*Bxrc);
	Axr = fc*(-Lrx*Bxxc + Lxx*Bxrc);
	Arx = fc*( Lrr*Brxc - Lxr*Brrc);
	Arr = fc*(-Lrx*Brxc + Lxx*Brrc);

	// calculate components of F matrix
	Fxx = t*I0*Axx + (omeg*I1 + 2*I0)*Axr;
	Fxr = t*I0*Arx + (omeg*I1 + 2*I0)*Arr;
	Frx = t*I1*Axx +  omeg*I0        *Axr;
	Frr = t*I1*Arx +  omeg*I0        *Arr;

	// regularize the xx integrand
	Fxx = Fxx + 8*K0n;

	// increment the integrands
//	mRxx +=  Bxx*cosXt;
//	mRxr +=  Brx*sinXt; // note that mRxr is related to Brx
//	mRrx +=  Bxr*sinXt; // and mRrx is related to Bxr
//	mRrr += -Brr*cosXt;

	mCxx +=  Fxx*cosXt/s;
	mCxr +=  Fxr*sinXt/s;
	mCrx +=  Frx*sinXt/s;
	mCrr += -Frr*cosXt/s;

//	MRxx = -4*r0*mRxx*dt;
//	MRxr = -4*r0*mRxr*dt;
//	MRrx = -4*r0*mRrx*dt;
//	MRrr = -4*r0*mRrr*dt;

}

/* velocity u evaluated at (x,r) due to a ring of point forces 
 * of strength f located at (x0,r0), bounded externally by a
 * cylindrical tube of radius rc.
 */
void gf_axT_vel(double x, double r, double x0, double r0, double rc,
                double fx, double fr, double &ux, double &ur){
	// declare variables
	double Mxx, Mxr, Mrx, Mrr;
	
	// calculate components of the Green's function
	gf_axT(x, r, x0, r0, rc, Mxx, Mxr, Mrx, Mrr);

	// calculate velocity components
	ux = (Mxx*fx + Mxr*fr)/(8*M_PI);
	ur = (Mrx*fx + Mrr*fr)/(8*M_PI);

}
