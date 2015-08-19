/* BOUNDARY INTEGRAL SOLVER
 *  Solve the boundary integral equations and time-evolve the
 *  quantities associated with the geometric nodes (e.g., position,
 *  concentration).
 *
 * REFERENCES
 *  
 * PARAMETERS
 *  x	[input/output]	description
 *  y [input/output]	description
 */


/* FOR NOW, FOCUS ON SINGLE LAYER IMPLEMENTATION AND
 * THE IMPLEMENTATION OF THE FREE BOUNDARY PROBLEM
 * (I.E., IMPLEMENTING THE KINEMATIC CONDITION) 
 * JUST USE A FORWARD EULER SCHEME FOR NOW */




#ifndef SOLVER_H
#define SOLVER_H

/* HEADER FILES */
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include "quad.h"

using namespace std;

/* PROTOTYPES */
void checkAngle(surface &, double, int &);
void checkSpacing(surface &, double, double, double);
void writeNode(int, int, double*, double*, string);

/* IMPLEMENTATIONS */

/* Time integration */
void timeInt(int nstep, int nquad, double dt, surface Surface, string opath){

  /*-----------------------------------------------------*/
  /*----------------------- SETUP -----------------------*/
  /*-----------------------------------------------------*/

	// declare variables
	int    i, j, n;
	int    istep;
	int    istop;
	int    nelem, ngeom;
	int    nlocl, nglob;

	double nx  , nr  ;
	double vx  , vr  , vn;
	double area, vlme;
	double thetmax, smax, smin;
	double gamm;
	
	int    IGF;
	int    ISURF;
	
  /* get number of boundary elements,
   * geometric nodes, and basis nodes 
	 * (local and global) */
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// initialize vectors
	vector<double> x  (ngeom), r  (ngeom);
	vector<double> v  (2*nglob);
	vector<double> Df (2*nglob);


  /*-----------------------------------------------------*/
  /*-------------------- INITIALIZE ---------------------*/
  /*-----------------------------------------------------*/

	IGF = 1;   // tube Green's function
	IGF = 0;   // free-space Green's function

	ISURF = 2; // vesicle
	ISURF = 0; // drop

	istop = 0;

	if (ISURF == 0){
		gamm = Surface.getMeanTens();
	}

	// set surface velocity
	for (i = 0; i < nglob; i++){
    Surface.setVel (i, 0.0, 0.0);
		v [2*i  ] = 0.;
		v [2*i+1] = 0.;
    
		Surface.setTrct(i, 0.0, 0.0);
		Df[2*i  ] = 0.;
		Df[2*i+1] = 0.;
	}
	
	// get initial geometry
	Surface.getNode(x .data(), r .data());
	area = Surface.getArea();
	vlme = Surface.getVlme();

	// extrema for node redistribution
	thetmax = M_PI/8;
	smax = 2.0*M_PI/nelem;
	smin = smax/5.0;
	
	
  /*-----------------------------------------------------*/
  /*----------------- TIME INTEGRATION ------------------*/
  /*-----------------------------------------------------*/
	
	for (istep = 0; istep < nstep; istep++){
		cout << "ts = " << istep << ", V = " << vlme << endl;
		
		/*-- Step 1: Solve the boundary integral equation. --*/
		singleLayer(IGF, nquad, Surface, v.data(), Df.data());

		// update kinematic and dynamic parameters
		for (i = 0; i < nglob; i++){
			Surface.setVel (i, v [2*i], v [2*i+1]);
			Surface.setTrct(i, Df[2*i], Df[2*i+1]);

			// SHOULD ALSO SET DISPLACEMENT AND CONCENTRATION FIELDS
		}

		// write to file before evolving system
		writeNode(istep, ngeom, x.data(), r.data(), opath);
		
		// NOTE: FOR SINGLE LAYER POTENTIAL, THERE IS NO NEED
		// TO CALCULATE A MATRIX INVERSE
		
		// NOTE: NEED TO CHECK VOLUME


		/*-- Step 2: Update surface velocity and boundary ---*
		 *---------- shape using the kinematic condition. ---*/
		
		for (i = 0; i < ngeom; i++){
			// get global index
			n = i*(nlocl - 1);

      // get normal at the ith geometric node
      Surface.getNrml(i, nx, nr);

			// calculate surface velocity
			vx  = v[2*n  ];
			vr  = v[2*n+1];
			vn  = vx*nx + vr*nr;
			
			// advect geometric nodes using forward Euler scheme
			x[i] += nx*vn*dt;
			r[i] += nr*vn*dt;
			
			// NOTE: SHOULD ALSO USE BACKWARD EULER SCHEME TO CHECK
			// ERROR.

		}

		/*-- Step 3: Update surface fields. -----------------*/
		Surface.setGeomParams(nelem, x.data(), r.data());

		// NOTE: WOULD HAVE TO UPDATE STOKES.H DISPLACEMENT,
		// VELOCITY, TRACTION, AND SURFACTANT CONCENTRATION FIELDS HERE
		Surface.setSurfParams(nelem, nlocl-1, gamm);

		/*-- Step 4: Check node separation and subtended angle
		 *---------- and redistribute nodes accordingly -----*/
		checkAngle  (Surface, thetmax, istop);
		checkSpacing(Surface, smin, smax, thetmax);
		
		// check to break loop
		if (istop == 1)
			return;
		
		
		/*-- Step 5: Get parameters for next timestep  ------*/
 		nelem = Surface.getNElem();
 		ngeom = nelem + 1;
 		nlocl = Surface.getNLocl();
 		nglob = Surface.getNGlob();
		
		Surface.getArea(area);
		Surface.getVlme(vlme);
	
		x.resize(ngeom);
    r.resize(ngeom);
		Surface.getNode(x.data(), r.data());
		
		if (2*nglob > v.size()){
			v .resize(2*nglob);
			Df.resize(2*nglob);
		}

		// didn't resize tensions/moments... so do that now
		Surface.setSurfParams(nelem, nlocl-1);

		// SHOULD WRITE VELOCITIES AND TRACTIONS AND SEE HOW THESE LOOK
		

	}
}

/*- NODE REDISTRIBUTION ---*/

/* Redistribute nodes on the contour if the angle 
 * subtended by an arc is too large.
 *
 * If i = 0, add another point. If i > 0, remove
 * the middle point and add two evenly spaced
 * points. */
void checkAngle(surface &Surface, double thetmax, int &istop){
	// declare variables
  int     i, i0, i1;
	int     j, j0, j1;
	int     k, k0, k1;
	int     m;
  int     nelem, ngeom;
	int     nlocl, nglob;
	int     counter;
	double  Dthet;

	double  ks , kp ;
	
	double  axM, bxM, cxM;
	double  arM, brM, crM;
	double  ax0, bx0, cx0;
	double  ar0, br0, cr0;
	double  ax1, bx1, cx1;
	double  ar1, br1, cr1;
	
	double  xM, rM, sM, lM, dlM; 
	double  x0, r0, s0, l0, dl0;
	double  x1, r1, s1, l1, dl1;

	double  li, li0, li1, lj;
	
	double  uxp, urp, vxp, vrp, fxp, frp, cp;

	double  cf, *zlocl;
	double  zp, *L;

	// constants
	const double PIH = 0.5*M_PI;

	// initialize
	istop   = 0;	// indicator for whether simulation should terminate
	counter = 0;	// counter for how many nodes added/removed

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// allocate memory
	zlocl = (double*) malloc (nlocl * sizeof(double));
	L     = (double*) malloc (nlocl * sizeof(double));

  // calculate Gauss-Lobatto points on the interval [-1,1]
  cf = M_PI/(nlocl-1);
  for (i = 0; i < nlocl; i++){
    zlocl[nlocl - i - 1] = gsl_sf_cos(cf*i);
  }

	// initialize containers (may be resized)
	vector<double>  x  (ngeom),  r  (ngeom);
	vector<double> ux  (nglob), ur  (nglob);
	vector<double> vx  (nglob), vr  (nglob);
	vector<double> fx  (nglob), fr  (nglob);
	vector<double>  c  (nglob)             ;

//	vector<double> taus(nglob), taup(nglob);
//	vector<double> q   (nglob)             ;
//	vector<double> ms  (nglob), mp  (nglob);

	// get node coordinates and basis functions
  Surface.getNode( x.data(),  r.data());
	Surface.getDisp(ux.data(), ur.data());
	Surface.getVel (vx.data(), vr.data());
	Surface.getTrct(fx.data(), fr.data());
	Surface.getConc( c.data()           );

	// NOTE: AND WHAT ABOUT CONSTITUTIVE FIELDS? (TENSIONS, MOMENTS, etc.)
	// (WORRY ABOUT THIS LATER)

  for (i = 0; i < nelem; i++){ // loop over boundary elements
		i0 = i - 1;
		i1 = i + 1;

  	if (i == 0){
			Surface.getArcl(1, sM);
			Surface.getCurv(0, ks, kp);
			Dthet = 2*sM*ks;
			Dthet = fabs(Dthet);
  	}
  	else {
			Surface.getArcl(i0, s0);
			Surface.getArcl(i1, s1);
			Surface.getCurv(i , ks, kp);
			Dthet = (s1 - s0)*ks;
			Dthet = fabs(Dthet);
  	}
		
		/*--------------------------------------------*/
		/*---------- CHECK EXCESSIVE ANGLE -----------*/
		/*--------------------------------------------*/

		if (Dthet > PIH){
			cout << "checkAngle: an arc is excessive at the " <<
				i1 << "th index" << endl;
			
			istop = 1;

			free(L    );
			free(zlocl);

			return;
		}
		
		/*--------------------------------------------*/
		/*----------- CHECK MAXIMUM ANGLE ------------*/
		/*--------------------------------------------*/
		
		if (Dthet > thetmax){
			counter++;
			
			if (i == 0){	// add a new node at the midpoint
				// temporary containers for new basis functions
				vector<double> uxnew(2*nlocl-3), urnew(2*nlocl-3);
				vector<double> vxnew(2*nlocl-3), vrnew(2*nlocl-3);
				vector<double> fxnew(2*nlocl-3), frnew(2*nlocl-3);
				vector<double>  cnew(2*nlocl-3)                  ;

				// get spline coefficients
				Surface.getSpln(i, axM, bxM, cxM,
				                   arM, brM, crM);

				// get midpoint
				Surface.getPoly(i , li );
				Surface.getPoly(i1, li1);

				 lM = (li + li1)/2.0;
				dlM = lM - li;
				 xM = ((axM*dlM + bxM)*dlM + cxM)*dlM + x[i];
				 rM = ((arM*dlM + brM)*dlM + crM)*dlM + r[i];
				
				// add new geometric node
				x.push_back(0.0);
				r.push_back(0.0);

				for (j = ngeom; j > 0; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[1] = xM;
				r[1] = rM;

				// NOTE: MAYBE HAVE if (nlocl > 2) statement? (might not be necessary...)
				//
				/* interpolate stokes fields to new basis nodes
				 * to the left of and including the midpoint */
				for (j = 1; j < nlocl; j++){   /* loop over new basis nodes to
				                                * the left of the midpoint 
																				* (including the midpoint) */
					j0 = j-1;

					// initialize
					zp = 0.5*(zlocl[j] + 1) - 1;
					lagrange(nlocl-1, zlocl, zp, L);
					uxp = 0; urp = 0;
					vxp = 0; vrp = 0;
					fxp = 0; frp = 0;
					 cp = 0;

					// interpolate functions to zp
					for (m = 0; m < nlocl; m++){
					//	k = i*(nlocl-1) + m;
						k = m;

						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
						 cp +=  c[k]*L[m];
					}

					// add new basis node
					uxnew[j0] = uxp; urnew[j0] = urp;
					vxnew[j0] = vxp; vrnew[j0] = vrp;
					fxnew[j0] = fxp; frnew[j0] = frp;
					 cnew[j0] =  cp;
				}
				
				/* interpolate stokes fields to new basis nodes
				 * to the right of the midpoint */
				for (j = 1; j < nlocl-1; j++){ /* loop over new basis nodes to
				                                * the right of the midpoint
																				* (excluding the midpoint) */
					j0 = nlocl + j - 2;

					// initialize
					zp = 0.5*(zlocl[j] + 1);
					lagrange(nlocl-1, zlocl, zp, L);
					uxp = 0; urp = 0;
					vxp = 0; vrp = 0;
					fxp = 0; frp = 0;
					 cp = 0;

					// interpolate functions to zp
					for (m = 0; m < nlocl; m++){
					//	k = i*(nlocl-1) + m;
						k = m;

						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
						 cp +=  c[k]*L[m];
					}

					// add new basis node
					uxnew[j0] = uxp; urnew[j0] = urp;
					vxnew[j0] = vxp; vrnew[j0] = vrp;
					fxnew[j0] = fxp; frnew[j0] = frp;
					 cnew[j0] =  cp;
				}

				/* delete original basis functions between the
				 * ith and (i+1)st geometric nodes and replace
				 * with the new interpolants */
				for (j = 1; j < nlocl-1; j++){ /* left of midpoint:
				                                * replace old functions */
					j0 = j - 1;
					
					ux[j] = uxnew[j0]; ur[j] = urnew[j0];
					vx[j] = vxnew[j0]; vr[j] = vrnew[j0];
					fx[j] = fxnew[j0]; fr[j] = frnew[j0];
					 c[j] =  cnew[j0];
				}

				for (j = 0; j < nlocl-1; j++){ /* right of midpoint
				                                * (including midpoint):
				                                * increase vector size,
																				* shift vector elements,
				                                * and add new functions */
					j0 = j + nlocl - 2;
					j1 = j + nlocl - 1;

					ux.push_back(0.0); ur.push_back(0.0);
					vx.push_back(0.0); vr.push_back(0.0);
					fx.push_back(0.0); fr.push_back(0.0);
					c .push_back(0.0);

					for (k = nglob+j; k > j0; k--){
						k1 = k+1;

						ux[k1] = ux[k]; ur[k1] = ur[k];
						vx[k1] = vx[k]; vr[k1] = vr[k];
						fx[k1] = fx[k]; fr[k1] = fr[k];
						 c[k1] =  c[k];
					}
					ux[j1] = uxnew[j0]; ur[j1] = urnew[j0];
					vx[j1] = vxnew[j0]; vr[j1] = vrnew[j0];
					fx[j1] = fxnew[j0]; fr[j1] = frnew[j0];
					 c[j1] =  cnew[j0];
				}
			}
			else {				/* remove the middle node and add
			               * two evenly spaced nodes */
				// temporary containers for new basis functions
				vector<double> uxnew(3*nlocl-4), urnew(3*nlocl-4);
				vector<double> vxnew(3*nlocl-4), vrnew(3*nlocl-4);
				vector<double> fxnew(3*nlocl-4), frnew(3*nlocl-4);
				vector<double>  cnew(3*nlocl-4)                  ;

				// get spline coefficients
				Surface.getSpln(i0, ax0, bx0, cx0,
				                    ar0, br0, cr0);
				Surface.getSpln(i , ax1, bx1, cx1,
				                    ar1, br1, cr1);

				// compute new points
				Surface.getPoly(i0, li0);
				Surface.getPoly(i , li );
				Surface.getPoly(i1, li1);

				 l0 = (li0 + 2.0*li)/3.0; 
				dl0 = l0 - li0;
				 x0 = ((ax0*dl0 + bx0)*dl0 + cx0)*dl0 + x[i0];
				 r0 = ((ar0*dl0 + br0)*dl0 + cr0)*dl0 + r[i0];

				 l1 = (2.0*li + li1)/3.0; 
				dl1 = l1 - li;
				 x1 = ((ax1*dl1 + bx1)*dl1 + cx1)*dl1 + x[i ];
				 r1 = ((ar1*dl1 + br1)*dl1 + cr1)*dl1 + r[i ];

				// add new geometric nodes
				x.push_back(0.0);
				r.push_back(0.0);

				for(j = ngeom; j > i-1; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[i]   = x1;
				r[i]   = r1;

				x[i-1] = x0;
				r[i-1] = r0;
			
				/* interpolate stokes fields to new basis nodes
				 * in the first element */
				for (j = 1; j < nlocl; j++){
					j0 = j-1;

					// calculate polygonal line segment
					lj = 0.5*dl0*(zlocl[j] + 1) + li0;

					// initialize
					uxp = 0; urp = 0;
					vxp = 0; vrp = 0;
					fxp = 0; frp = 0;
					 cp = 0;

					if (lj < li){
						// compute Lagrange polynomials
						zp = 2*(lj - li0)/(li - li0) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i0*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}
					else {
						// compute Lagrange polynomials
						zp = 2*(lj - li)/(li - li) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}

					// add new basis node
					uxnew[j0] = uxp; urnew[j0] = urp;
					vxnew[j0] = vxp; vrnew[j0] = vrp;
					fxnew[j0] = fxp; frnew[j0] = frp;
					 cnew[j0] =  cp;
				}

				/* interpolate stokes fields to new basis nodes
				 * in the second element */
				dlM = l1 - l0;
				for (j = 1; j < nlocl; j++){
					j0 = nlocl + j - 2;
					
					// calculate polygonal line segment
					lj = 0.5*dlM*(zlocl[j] + 1) + l0;
					
					// initialize
					uxp = 0; urp = 0;
					vxp = 0; vrp = 0;
					fxp = 0; frp = 0;
					 cp = 0;

					if (lj < li){
						// compute Lagrange polynomials
						zp = 2*(lj - li0)/(li - li0) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i0*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}
					else {
						// compute Lagrange polynomials
						zp = 2*(lj - li)/(li - li) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}

					// add new basis node
					uxnew[j0] = uxp; urnew[j0] = urp;
					vxnew[j0] = vxp; vrnew[j0] = vrp;
					fxnew[j0] = fxp; frnew[j0] = frp;
					 cnew[j0] =  cp;
				}
				
				/* interpolate stokes fields to new basis nodes
				 * in the third element (exclude z = 1) */
				for (j = 1; j < nlocl-1; j++){
					j0 = 2*nlocl + j - 3;
					
					// calculate polygonal line segment
					lj = 0.5*dl1*(zlocl[j] + 1) + li1;
					
					// initialize
					uxp = 0; urp = 0;
					vxp = 0; vrp = 0;
					fxp = 0; frp = 0;
					 cp = 0;

					if (lj < li){
						// compute Lagrange polynomials
						zp = 2*(lj - li0)/(li - li0) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i0*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}
					else {
						// compute Lagrange polynomials
						zp = 2*(lj - li)/(li - li) - 1;
						lagrange(nlocl-1, zlocl, zp, L);
	
						// interpolate functions to zp
						for (m = 0; m < nlocl; m++){
							k = i*(nlocl-1) + m;
	
							uxp += ux[k]*L[m]; urp += ur[k]*L[m];
							vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
							fxp += fx[k]*L[m]; frp += fr[k]*L[m];
							 cp +=  c[k]*L[m];
						}
					}

					// add new basis node
					uxnew[j0] = uxp; urnew[j0] = urp;
					vxnew[j0] = vxp; vrnew[j0] = vrp;
					fxnew[j0] = fxp; frnew[j0] = frp;
					 cnew[j0] =  cp;
				}

				/* delete original basis functions between the
				 * (i-1)st and ith geometric nodes and replace
				 * with interpolants in the first new element */
				for (j = 1; j < nlocl; j++){
					j0 = j-1;
					k  = i0*(nlocl-1) + j;
					
					ux[k] = uxnew[j0]; ur[k] = urnew[j0];
					vx[k] = vxnew[j0]; vr[k] = vrnew[j0];
					fx[k] = fxnew[j0]; fr[k] = frnew[j0];
					 c[k] =  cnew[j0];
				}

				/* delete original basis functions between the
				 * ith and (i+1)st geometric nodes and replace
				 * with interpolants in the second new element */
				for (j = 1; j < nlocl-1; j++){
					j0 = nlocl + j - 2;
					k  = i*(nlocl-1) + j;
					
					ux[k] = uxnew[j0]; ur[k] = urnew[j0];
					vx[k] = vxnew[j0]; vr[k] = vrnew[j0];
					fx[k] = fxnew[j0]; fr[k] = frnew[j0];
					 c[k] =  cnew[j0];
				}

				/* add interpolants in the third element */
				for (j = 0; j < nlocl-1; j++){
					j0 = 2*nlocl + j - 3;
					k0 = i1*(nlocl-1) + j - 1;
					
					ux.push_back(0.0); ur.push_back(0.0);
					vx.push_back(0.0); vr.push_back(0.0);
					fx.push_back(0.0); fr.push_back(0.0);
					c .push_back(0.0);

					for (k = nglob + j; k > k0; k--){
						k1 = k+1;

						ux[k1] = ux[k]; ur[k1] = ur[k];
						vx[k1] = vx[k]; vr[k1] = vr[k];
						fx[k1] = fx[k]; fr[k1] = fr[k];
						 c[k1] =  c[k];
					}
					ux[k0+1] = uxnew[j0]; ur[k0+1] = urnew[j0];
					vx[k0+1] = vxnew[j0]; vr[k0+1] = vrnew[j0];
					fx[k0+1] = fxnew[j0]; fr[k0+1] = frnew[j0];
					 c[k0+1] =  cnew[j0];
				}
			}

			// update geometric nodes
			nelem++;
			ngeom++;
			Surface.geomPushBack();
			Surface.setNNode(ngeom);
			Surface.setNElem(nelem);
			Surface.setGeomParams(nelem, x.data(), r.data());
	
			// update stokes fields
			nglob = nelem*(nlocl-1) + 1;
			Surface.setNGlob(nglob);
			for (j = 0; j < nlocl-1; j++){
				Surface.stksPushBack();
			}
			Surface.setStksParams(nelem, nlocl-1,
			                      ux.data(), ur.data(),
			                      vx.data(), vr.data(),
			                      fx.data(), fr.data(),
			                      c .data());

			cout << "1 node added. Total of "
				<< nelem << " boundary elements." << endl;
		
			free(zlocl);
			free(L    );
		
			return;
		}
	}

//	if (counter > 0){
//		cout << counter << " node(s) added. Total of "
//			<< nelem << " boundary elements." << endl;
//	}

	free(zlocl);
	free(L    );
}

/* Redistribute nodes on the contour if two nodes 
 * are too close or too far apart. 
 * 
 * If a segment is too long, add a point in the middle.
 * If a segment is too short, remove the ith and (i+1)st
 * points and replace with a midpoint. */
void checkSpacing(surface &Surface, double smin, double smax, double thetmax){
	// declare variables
  int     i, im, i0, i1, i2, i3;
	int     j, j0, j1;
	int     k, k0, k1;
	int     m;
  int     nelem, ngeom;
	int     nlocl, nglob;
	int     counter;
	double  Dthet;
	double  Ds;

	double  ks , kp ;
	double  ksM, kpM;
	double  ks0, kp0;
	double  ks1, kp1;
	double  ks2, kp2;
	
	double  axM, bxM, cxM;
	double  arM, brM, crM;
	double  ax0, bx0, cx0;
	double  ar0, br0, cr0;
	double  ax1, bx1, cx1;
	double  ar1, br1, cr1;
	
	double  xM, rM, sM, lM, dlM; 
	double  dl0, dl1;

	double  li, li0, li1, li2, lj;
	double  sim, si, si0, si1, si2, si3;
	
	double  uxp, urp, vxp, vrp, fxp, frp, cp;

	double  cf, *zlocl;
	double  zp, *L;

	// constants
	const double PIH = 0.5*M_PI;

	// initialize
	counter = 0;	// counter for how many nodes added/removed

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// allocate memory
	zlocl = (double*) malloc (nlocl * sizeof(double));
	L     = (double*) malloc (nlocl * sizeof(double));

  // calculate Gauss-Lobatto points on the interval [-1,1]
  cf = M_PI/(nlocl-1);
  for (i = 0; i < nlocl; i++){
    zlocl[nlocl - i - 1] = gsl_sf_cos(cf*i);
  }

	// initialize containers (may be resized)
	vector<double>  x  (ngeom),  r  (ngeom);
	vector<double> ux  (nglob), ur  (nglob);
	vector<double> vx  (nglob), vr  (nglob);
	vector<double> fx  (nglob), fr  (nglob);
	vector<double>  c  (nglob)             ;

	// get node coordinates and basis functions
  Surface.getNode( x.data(),  r.data());
	Surface.getDisp(ux.data(), ur.data());
	Surface.getVel (vx.data(), vr.data());
	Surface.getTrct(fx.data(), fr.data());
	Surface.getConc( c.data()           );

	// NOTE: AND WHAT ABOUT CONSTITUTIVE FIELDS? (TENSIONS, MOMENTS, etc.)
	// (WORRY ABOUT THIS LATER)

  for (i = 0; i < nelem; i++){ // loop over boundary elements
		im = i - 2;
		i0 = i - 1;
		i1 = i + 1;
		i2 = i + 2;
		i3 = i + 3;

		// get separation
		Surface.getArcl(i , si );
		Surface.getArcl(i1, si1);

		Ds = si1 - si;

		/*--------------------------------------------*/
		/*--------- CHECK MAXIMUM SEPARATION ---------*/
		/*--------------------------------------------*/

		if (Ds > smax){						 /* if segment length exceeds
															  * maximum separation, add
																* a point in the middle */
			counter++;

			// temporary containers for new basis functions
			vector<double> uxnew(2*nlocl-3), urnew(2*nlocl-3);
			vector<double> vxnew(2*nlocl-3), vrnew(2*nlocl-3);
			vector<double> fxnew(2*nlocl-3), frnew(2*nlocl-3);
			vector<double>  cnew(2*nlocl-3)                  ;

			// get midpoint
			Surface.getSpln(i, axM, bxM, cxM,
			                   arM, brM, crM);
			Surface.getPoly(i , li );
			Surface.getPoly(i1, li1);
			lM = 0.5*(li1 + li);
			dlM = lM - li;
			xM = ((axM*dlM + bxM)*dlM + cxM)*dlM + x[i];
			rM = ((arM*dlM + brM)*dlM + crM)*dlM + r[i];

			// add new geometric node
			x.push_back(0.0);
			r.push_back(0.0);

			for (j = ngeom; j >i-1; j--){
				j1 = j+1;

				x[j1] = x[j];
				r[j1] = r[j];
			}
			x[i] = xM;
			r[i] = rM;

			/* interpolate stokes fields to new basis nodes
		   * to the left of and including the midpoint */
			for (j = 1; j < nlocl; j++){	/* loop over new basis nodes to
			                               * the left of the midpoint 
																		 * (including the midpoint */
				j0 = j-1;

				// initialize
				zp = 0.5*(zlocl[j] + 1) - 1;
				lagrange(nlocl-1, zlocl, zp, L);
				uxp = 0; urp = 0;
				vxp = 0; vrp = 0;
				fxp = 0; frp = 0;
				 cp = 0;

				// interpolate functions to zp
				for (m = 0; m < nlocl; m++){
					k = i*(nlocl - 1) + m;

					uxp += ux[k]*L[m]; urp += ur[k]*L[m];
					vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
					fxp += fx[k]*L[m]; frp += fr[k]*L[m];
					 cp +=  c[k]*L[m];
				}

				// add new basis node
				uxnew[j0] = uxp; urnew[j0] = urp;
				vxnew[j0] = vxp; vrnew[j0] = vrp;
				fxnew[j0] = fxp; frnew[j0] = frp;
				 cnew[j0] =  cp;
			}

			/* interpolate stokes fields to new basis nodes
			 * to the right of the midpoint */
			for (j = 1; j < nlocl-1; j++){ /* loop over new basis nodes to
			                                * the right of the midpoint
																			* (excluding the midpoint */
				j0 = nlocl + j - 2;

				// initialize
				zp = 0.5*(zlocl[j] + 1);
				lagrange(nlocl-1, zlocl, zp, L);
				uxp = 0; urp = 0;
				vxp = 0; vrp = 0;
				fxp = 0; frp = 0;
				 cp = 0;

				// interpolate functions to zp
				for (m = 0; m < nlocl; m++){
					k = i*(nlocl - 1) + m;

					uxp += ux[k]*L[m]; urp += ur[k]*L[m];
					vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
					fxp += fx[k]*L[m]; frp += fr[k]*L[m];
					 cp +=  c[k]*L[m];
				}

				// add new basis node
				uxnew[j0] = uxp; urnew[j0] = urp;
				vxnew[j0] = vxp; vrnew[j0] = vrp;
				fxnew[j0] = fxp; frnew[j0] = frp;
				 cnew[j0] =  cp;
			}

			/* delete original basis functions between the
			 * ith and (i+1)st geometric nodes and replace
			 * with the new interpolants */
			for (j = 1; j < nlocl-1; j++){ /* left of midpoint:
			                                * replace old functions */
				j0 = j - 1;
				k  = i*(nlocl-1) + j;

				ux[k] = uxnew[j0]; ur[k] = urnew[j0];
				vx[k] = vxnew[j0]; vr[k] = vrnew[j0];
				fx[k] = fxnew[j0]; fr[k] = frnew[j0];
				 c[k] =  cnew[j0];
			}

			for (j = 0; j < nlocl-1; j++){ /* right of midpoint
																			* (including midpoint):
			                                * increase vector size,
																			* shift vector elements,
																			* and add new functions */
				j0 = j + nlocl - 2;
				j1 = j + nlocl - 1;
				k0 = i1*(nlocl-1) + j - 1;

				ux.push_back(0.0); ur.push_back(0.0);
				vx.push_back(0.0); vr.push_back(0.0);
				fx.push_back(0.0); fr.push_back(0.0);
				 c.push_back(0.0);

				for (k = j+nglob; k > k0; k--){
					k1 = k + 1;

					ux[k1] = ux[k]; ur[k1] = ur[k];
					vx[k1] = vx[k]; vr[k1] = vr[k];
					fx[k1] = fx[k]; fr[k1] = fr[k];
					 c[k1] =  c[k];
				}
				ux[k0+1] = uxnew[j0]; ur[k0+1] = ur[j0];
				vx[k0+1] = vxnew[j0]; vr[k0+1] = vr[j0];
				fx[k0+1] = fxnew[j0]; fr[k0+1] = fr[j0];
				 c[k0+1] =  cnew[j0];
			}

			// update geometric nodes
			nelem++;
			ngeom++;
			Surface.geomPushBack();
			Surface.setNNode(ngeom);
			Surface.setNElem(nelem);
			Surface.setGeomParams(nelem, x.data(), r.data());
	
			// update stokes fields
			nglob = nelem*(nlocl-1) + 1;
			Surface.setNGlob(nglob);
			for (j = 0; j < nlocl-1; j++){
				Surface.stksPushBack();
			}
			Surface.setStksParams(nelem, nlocl-1,
			                      ux.data(), ur.data(),
			                      vx.data(), vr.data(),
			                      fx.data(), fr.data(),
			                      c .data());

			cout << "1 node added. Total of "
				<< nelem << " boundary elements." << endl;
		
			free(zlocl);
			free(L    );
		
			return;
		} // if (Ds > smax)
		
		/*--------------------------------------------*/
		/*--------- CHECK MINIMUM SEPARATION ---------*/
		/*--------------------------------------------*/

		if (Ds < smin){						 /* if segment length falls
                    					  * below minimum separation,
                    						* remove both endpoints and
																* replace with the midpoint */
			/* NOTE: A short segment is eliminated ONLY if this
			 * action does not violate the pre-established
			 * requirements (maximum angle, maximum separation). */
		
			if (i > 0 && i < nelem-1){

			counter--;
			
			double sep1, sep2;
			double totang0, totang1, totang2;

		//	if (i == 0){
		//		i0 = i1;
		//		im = i2;
		//	}
			if (i == 1){
				im = i ;
			}
			if (i == nelem-2){
				i3 = i1;
			}
		//	if (i == nelem-1){
		//		i2 = i ;
		//		i3 = i0;
		//	}
			
			// get additional geometric parameters
			Surface.getCurv(i , ks , kp );
			Surface.getCurv(i0, ks0, kp0);
			Surface.getCurv(i1, ks1, kp1);
			Surface.getCurv(i2, ks2, kp2);
			Surface.getArcl(im, sim);
			Surface.getArcl(i0, si0);
			Surface.getArcl(i2, si2);
			Surface.getArcl(i3, si3);

		//	if (i == 0){
		//		si0 = -si0;
		//		sim = -sim;
		//	}
			if (i == 1){
				sim = -sim;
			}
			if (i == nelem-2){
				si3 = 2*si1 - si3;
			}
		//	if (i == nelem-1){
		//		si2 = 2*si1 - si2;
		//		si3 = 2*si1 - si3;
		//	}
			
			sM  = 0.5*(si + si1);
			ksM = 0.5*(ks + ks1);
			
			// get new separations
			sep1 = sM  - si0;
			sep2 = si2 - sM ;
			
			// get new angles
			totang0 = ks0*(sM  - sim);
			totang1 = ksM*(si2 - si0);
			totang2 = ks2*(si3 - sM );

			if (totang0 < thetmax &&
			    totang1 < thetmax &&
					totang2 < thetmax &&
					sep1    < smax    &&
					sep2    < smax    ){
				// temporary containers for new basis functions
				vector<double> uxnew(2*nlocl-3), urnew(2*nlocl-3);
				vector<double> vxnew(2*nlocl-3), vrnew(2*nlocl-3);
				vector<double> fxnew(2*nlocl-3), frnew(2*nlocl-3);
				vector<double>  cnew(2*nlocl-3)                  ;
				
				// get midpoint
				Surface.getSpln(i, axM, bxM, cxM,
				                   arM, brM, crM);
				Surface.getPoly(i , li );
				Surface.getPoly(i1, li1);
				lM = 0.5*(li1 + li);
				dlM = lM - li;
				xM = ((axM*dlM + bxM)*dlM + cxM)*dlM + x[i];
				rM = ((arM*dlM + brM)*dlM + crM)*dlM + r[i];
				
				Surface.getPoly(i0, li0);
				Surface.getPoly(i2, li2);
				dl0 = lM  - li0;
				dl1 = li1 - lM ;

				/* remove geometric nodes and 
				 * replace with midpoint */
				for (j = i1; j < nelem; j++){
					j1 = j + 1;
					
					x[j] = x[j1];
					r[j] = r[j1];
				}
				x[i] = xM;
				r[i] = rM;

	//			/* interpolate stokes fields to new basis nodes
	//		   * to the left of and including the midpoint */
	//			for (j = 1; j < nlocl; j++){	/* loop over new basis nodes to
	//			                               * the left of the midpoint 
	//																		 * (including the midpoint */
	//				j0 = j-1;

	//				// get polygonal arclength
	//				lj = 0.5*dl0*(zlocl[j] + 1) + li0;

	//				// interpolate in the appropriate element
	//				if (lj < li){
	//					zp = 2*(lj - li0)/(li  - li0) - 1;
	//					lagrange(nlocl-1, zlocl, zp, L);
	//					uxp = 0; urp = 0;
	//					vxp = 0; vrp = 0;
	//					fxp = 0; frp = 0;
	//					 cp = 0;

	//					for (m = 0; m < nlocl; m++){
	//						k = i0*(nlocl - 1) + m;

	//						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
	//						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
	//						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
	//						 cp +=  c[k]*L[m];
	//					}
	//				}
	//				else {
	//					zp = 2*(lj - li )/(li1 - li ) - 1;
	//					lagrange(nlocl-1, zlocl, zp, L);
	//					uxp = 0; urp = 0;
	//					vxp = 0; vrp = 0;
	//					fxp = 0; frp = 0;
	//					 cp = 0;
	//					
	//					for (m = 0; m < nlocl; m++){
	//						k = i *(nlocl - 1) + m;

	//						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
	//						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
	//						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
	//						 cp +=  c[k]*L[m];
	//					}
	//				}

	//				// add new basis node
	//				uxnew[j0] = uxp; urnew[j0] = urp;
	//				vxnew[j0] = vxp; vrnew[j0] = vrp;
	//				fxnew[j0] = fxp; frnew[j0] = frp;
	//				 cnew[j0] =  cp;
	//			}

	//			/* interpolate stokes fields to new basis nodes
	//			 * to the right of the midpoint */
	//			for (j = 1; j < nlocl; j++){	/* loop over new basis nodes to
	//			                               * the right of the midpoint 
	//																		 * (excluding the midpoint */
	//				j0 = nlocl + j - 2;

	//				// get polygonal arclength
	//				lj = 0.5*dl1*(zlocl[j] + 1) + li0;

	//				// interpolate in the appropriate element
	//				if (lj < li1){
	//					zp = 2*(lj - li )/(li1 - li ) - 1;
	//					lagrange(nlocl-1, zlocl, zp, L);
	//					uxp = 0; urp = 0;
	//					vxp = 0; vrp = 0;
	//					fxp = 0; frp = 0;
	//					 cp = 0;
	//					
	//					for (m = 0; m < nlocl; m++){
	//						k = i *(nlocl - 1) + m;

	//						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
	//						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
	//						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
	//						 cp +=  c[k]*L[m];
	//					}
	//				}
	//				else {
	//					zp = 2*(lj - li1)/(li2 - li1) - 1;
	//					lagrange(nlocl-1, zlocl, zp, L);
	//					uxp = 0; urp = 0;
	//					vxp = 0; vrp = 0;
	//					fxp = 0; frp = 0;
	//					 cp = 0;
	//					
	//					for (m = 0; m < nlocl; m++){
	//						k = i1*(nlocl - 1) + m;

	//						uxp += ux[k]*L[m]; urp += ur[k]*L[m];
	//						vxp += vx[k]*L[m]; vrp += vr[k]*L[m];
	//						fxp += fx[k]*L[m]; frp += fr[k]*L[m];
	//						 cp +=  c[k]*L[m];
	//					}
	//				}

	//				// add new basis node
	//				uxnew[j0] = uxp; urnew[j0] = urp;
	//				vxnew[j0] = vxp; vrnew[j0] = vrp;
	//				fxnew[j0] = fxp; frnew[j0] = frp;
	//				 cnew[j0] =  cp;
	//			}

	//			/* delete original basis functions between the
	//			 * (i-1)st and (i+2)nd geometric nodes and replace
	//			 * with the new interpolants */
	//			for (j = 1; j < nlocl-1; j++){ /* left of midpoint:
	//			                                * replace old functions */
	//				j0 = j - 1;
	//				k  = i0*(nlocl-1) + j;

	//				ux[k] = uxnew[j0]; ur[k] = urnew[j0];
	//				vx[k] = vxnew[j0]; vr[k] = vrnew[j0];
	//				fx[k] = fxnew[j0]; fr[k] = frnew[j0];
	//				 c[k] =  cnew[j0];
	//			}

	//			for (j = 0; j < nlocl-1; j++){ /* right of midpoint
	//																			* (including midpoint):
	//			                                * replace old functions */
	//				j0 = j + nlocl - 2;
	//				k  = i *(nlocl-1) + j;

	//				ux[k] = uxnew[j0]; ur[k] = ur[j0];
	//				vx[k] = vxnew[j0]; vr[k] = vr[j0];
	//				fx[k] = fxnew[j0]; fr[k] = fr[j0];
	//				 c[k] =  cnew[j0];
	//			}

	//			for (j = 0; j < nlocl-1; j++){ /* delete local basis nodes 
	//																			* between (i+1)st and (i+2)nd
	//			                                * geometric nodes */
	//				j0 = j + nlocl - 2;
	//				k0 = i1*(nlocl-1);

	//				for (k = k0; k < nglob-j; k++) {
	//					k1 = k + 1;

	//					ux[k] = ux[k1]; ur[k] = ur[k1];
	//					vx[k] = vx[k1]; vr[k] = vr[k1];
	//					fx[k] = fx[k1]; fr[k] = fr[k1];
	//					 c[k] =  c[k1];
	//				}

	//				ux.pop_back(); ur.pop_back();
	//				vx.pop_back(); vr.pop_back();
	//				fx.pop_back(); fr.pop_back();
	//				 c.pop_back();
	//			}

				// update geometric nodes
				nelem--;
				ngeom--;
				Surface.geomPopBack();
				Surface.setNNode(ngeom);
				Surface.setNElem(nelem);
				Surface.setGeomParams(nelem, x.data(), r.data());
				
				// update stokes fields
				nglob = nelem*(nlocl-1) + 1;
				Surface.setNGlob(nglob);
				for (j = 0; j < nlocl-1; j++){
					Surface.stksPopBack();
				}
				Surface.setStksParams(nelem, nlocl-1,
				                      ux.data(), ur.data(),
				                      vx.data(), vr.data(),
				                      fx.data(), fr.data(),
				                      c .data());

				cout << "1 node removed. Total of "
					<< nelem << " boundary elements." << endl;
				
				free(zlocl);
				free(L    );

				return;
			} // if (totang and sep meet criteria)
			} // if (i > 0 && i < nelem-1)
		} // if (Ds < smin)
	} // end of boundary elements


//	if (counter > 0){
//		cout << counter << " node(s) added/removed. Total of "
//			<< nelem << " boundary elements." << endl;
//	}

	free(zlocl);
	free(L    );
}



/*- WRITE FUNCTIONS -------*/

/* Write geometric coordinates to file */
void writeNode(int istep, int ngeom, double *x, double *r, string path){
	// declare variables
	int i;
	ostringstream num;
	ofstream file;
	string id, fn, line;

	// get timestep
	num << istep;
	if (istep < 10)
		id = "0000" + num.str();
	else if (istep < 100)
		id = "000"  + num.str();
	else if (istep < 1000)
		id = "00"   + num.str();
	else if (istep < 10000)
		id = "0"    + num.str();
	
	// define filename
	fn = path + "nodexr" + id + ".txt";
	
	// write to file
	file.open(fn.c_str());
	for (i = 0; i < ngeom; i++){
		ostringstream xs, rs;
		xs << fixed << setprecision(8) << x[i];
		rs << fixed << setprecision(8) << r[i];
		line = xs.str() + " " + rs.str();
		
		if (x[i] > 0)
			line = " " + line;

		file << line << endl;
	}
	file.close();
}



#endif
