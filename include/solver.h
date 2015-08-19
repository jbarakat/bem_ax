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
	
	// extrema for node redistribution
	thetmax = M_PI/8;
	smax = 2.0*M_PI/nelem;
	smin = smax/4.0;
	
	
  /*-----------------------------------------------------*/
  /*----------------- TIME INTEGRATION ------------------*/
  /*-----------------------------------------------------*/
	
	for (istep = 0; istep < nstep; istep++){
		/*-- Step 0: Initialize -----------------------------*/
		vector<double> x (ngeom), r (ngeom);
		vector<double> v (2*nglob, 0.0);
		vector<double> Df(2*nglob, 0.0);
	
		// get geometry
		Surface.getNode(x.data(), r.data());
		area = Surface.getArea();
		vlme = Surface.getVlme();
		
		cout << "ts = " << istep << ", V = " << vlme << endl;
		
		/*-- Step 1: Solve the boundary integral equation. --*/
		singleLayer(IGF, nquad, Surface, v.data(), Df.data());

	//	// update kinematic and dynamic parameters
	//	for (i = 0; i < nglob; i++){
	//		Surface.setVel (i, v [2*i], v [2*i+1]);
	//		Surface.setTrct(i, Df[2*i], Df[2*i+1]);

	//		// SHOULD ALSO SET DISPLACEMENT AND CONCENTRATION FIELDS
	//	}

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
	//	Surface.setSurfParams(nelem, nlocl-1, gamm);

		/*-- Step 4: Check node separation and subtended angle
		 *---------- and redistribute nodes accordingly -----*/
		checkAngle  (Surface, thetmax, istop);
		
		// check to break loop
		if (istop == 1)
			return;
		
		checkSpacing(Surface, smin, smax, thetmax);
		
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
		
		//Surface.setSurfParams(nelem, nlocl-1);
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

	double  li, li0, li1;

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

	// initialize containers (may be resized)
	vector<double>  x  (ngeom),  r  (ngeom);

	// get node coordinates
  Surface.getNode( x.data(),  r.data());

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

			return;
		}
		
		/*--------------------------------------------*/
		/*----------- CHECK MAXIMUM ANGLE ------------*/
		/*--------------------------------------------*/
		
		if (Dthet > thetmax){
			counter++;
			
			if (i == 0){	// add a new node at the midpoint
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
			}
			else {				/* remove the middle node and add
			               * two evenly spaced nodes */
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

				for(j = ngeom; j > i; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[i1] = x1;
				r[i1] = r1;

				x[i ] = x0;
				r[i ] = r0;
			}

			// update geometry
			nelem++;
			ngeom++;
			nglob = nelem*(nlocl-1) + 1;
			Surface.geomPushBack();
			Surface.setNNode(ngeom);
			Surface.setNElem(nelem);
			Surface.setNGlob(nglob);
			Surface.setGeomParams(nelem, x.data(), r.data());
	
			cout << "checkAngle: 1 node added. Total of "
				<< nelem << " boundary elements." << endl;
		
			return;
		}
	}

//	if (counter > 0){
//		cout << counter << " node(s) added. Total of "
//			<< nelem << " boundary elements." << endl;
//	}
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
			
	double sep1, sep2;
	double totang0, totang1, totang2;

	// constants
	const double PIH = 0.5*M_PI;

	// initialize
	counter = 0;	// counter for how many nodes added/removed

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// initialize containers (may be resized)
	vector<double>  x  (ngeom),  r  (ngeom);

	// get node coordinates
  Surface.getNode( x.data(),  r.data());

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

			for (j = ngeom; j > i; j--){
				j1 = j+1;

				x[j1] = x[j];
				r[j1] = r[j];
			}
			x[i1] = xM;
			r[i1] = rM;

			// update geometric nodes
			nelem++;
			ngeom++;
			nglob = nelem*(nlocl-1) + 1;
			Surface.geomPushBack();
			Surface.setNNode(ngeom);
			Surface.setNElem(nelem);
			Surface.setNGlob(nglob);
			Surface.setGeomParams(nelem, x.data(), r.data());

			cout << "checkSpacing: 1 node added. Total of "
				<< nelem << " boundary elements." << endl;
		
			return;
		} // if (Ds > smax)
		
		/*--------------------------------------------*/
		/*--------- CHECK MINIMUM SEPARATION ---------*/
		/*--------------------------------------------*/

		// START HERE

		if (Ds < smin && Ds == 0){						 /* if segment length falls
                    					  * below minimum separation,
                    						* remove both endpoints and
																* replace with the midpoint */

			/* NOTE: A short segment is eliminated ONLY if this
			 * action does not violate the pre-established
			 * requirements (maximum angle, maximum separation). */
		
			if (i > 0 && i < nelem-1){

			counter--;

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

				// update geometric nodes
				nelem--;
				ngeom--;
				nglob = nelem*(nlocl-1) + 1;
				Surface.geomPopBack();
				Surface.setNNode(ngeom);
				Surface.setNElem(nelem);
				Surface.setNGlob(nglob);
				Surface.setGeomParams(nelem, x.data(), r.data());
				
				cout << "checkSpacing: 1 node removed. Total of "
					<< nelem << " boundary elements." << endl;

				return;
			} // if (totang and sep meet criteria)
			} // if (i > 0 && i < nelem-1)
		} // if (Ds < smin)
	} // end of boundary elements


//	if (counter > 0){
//		cout << counter << " node(s) added/removed. Total of "
//			<< nelem << " boundary elements." << endl;
//	}
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
