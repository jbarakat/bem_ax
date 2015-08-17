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
#include "quad.h"
#include <vector>

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
	double thetmax, lmax, lmin;
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
	lmax = 2.0*M_PI/nelem;
	lmin = lmax/2.0;
	
	
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
	//	checkAngle  (Surface, thetmax, istop);
	//	checkSpacing(Surface, lmin, lmax, thetmax);
		
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
  int     i, j, i0, i1, j1;
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

	double  cf, *zlocl;

	// constants
	const double PIH = 0.5*M_PI;

	// initialize
	istop = 0;
	counter = 0;

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// allocate memory
	zlocl = (double*) malloc (nelem * sizeof(double));

  // calculate Gauss-Lobatto points on the interval [-1,1]
  cf = M_PI/(nlocl-1);
  for (i = 0; i < nlocl; i++){
    zlocl[nlocl - i - 1] = gsl_sf_cos(cf*i);
  }

	
	/* initialize containers (may be resized) */
	vector<double> x   (ngeom), r   (ngeom);

	vector<double> ux  (nglob), ur  (nglob);
	vector<double> vx  (nglob), vr  (nglob);
	vector<double> fx  (nglob), fr  (nglob);
	vector<double> c   (nglob)             ;

	vector<double> taus(nglob), taup(nglob);
	vector<double> q   (nglob)             ;
	vector<double> ms  (nglob), mp  (nglob);

	// get node coordinates and basis functions
  Surface.getNode(x .data(), r .data());
	Surface.getVel (vx.data(), vx.data());
	Surface.getTrct(fx.data(), fx.data());

	// WILL ALSO NEED TO UPDATE DISPLACEMENT AND
	// CONCENTRATION

	// loop over boundary elements
  for (i = 0; i < nelem; i++){
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

		if (Dthet > PIH){
			cout << "checkAngle: an arc is excessive at the " <<
				(i+1) << "th index" << endl;
			
			istop = 1;
			return;
		}
		
		if (Dthet > thetmax){
			counter++;

			// resize containers
			x.push_back(0.0);
			r.push_back(0.0);

			for (j = 0; j < nlocl-1; j++){
				ux  .push_back(0.0);
				ur  .push_back(0.0);
				vx  .push_back(0.0);
				vr  .push_back(0.0);
				fx  .push_back(0.0);
				fr  .push_back(0.0);
				c   .push_back(0.0);

				taus.push_back(0.0);
				taup.push_back(0.0);
				q   .push_back(0.0);
				ms  .push_back(0.0);
				mp  .push_back(0.0);
			}

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

				// interpolate stokes fields to midpoint

				// ALSO NEED TO INTERPOLATE TO NEW LOCAL POINTS CREATED WITHIN EACH ELEMENT
				// START FROM HERE




				// add new point
				for (j = ngeom; j > 0; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[1] = xM;
				r[1] = rM;
	
				// update geometric nodes
				nelem++;
				ngeom++;
				nglob = nelem*(nlocl-1) + 1;
				Surface.geomPushBack();
				Surface.setNNode(ngeom);
				Surface.setNElem(nelem);
				Surface.setGeomParams(nelem, x.data(), r.data());

				// interpolate stokes fields
				


				// update stokes
				for (j = 0; j < nlocl-1; j++){
					Surface.stksPushBack();
					Surface.surfPushBack();
				}






				// calculate surface tensions and moments


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

				 l0 = (li0 + 2.0*li)/3.0; 
				dl0 = l0 - li;
				 x0 = ((ax0*dl0 + bx0)*dl0 + cx0)*dl0 + x[i0];
				 r0 = ((ar0*dl0 + br0)*dl0 + cr0)*dl0 + r[i0];

				Surface.getPoly(i , li );
				Surface.getPoly(i1, li1);

				 l1 = (2.0*li + li1)/3.0; 
				dl1 = l1 - li;
				 x1 = ((ax1*dl1 + bx1)*dl1 + cx1)*dl1 + x[i ];
				 r1 = ((ar1*dl1 + br1)*dl1 + cr1)*dl1 + r[i ];

				// add new points
				for(j = ngeom; j > i-1; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[i]   = x1;
				r[i]   = r1;

				x[i-1] = x0;
				r[i-1] = r0;
				
				// update geometric nodes
				nelem++;
				ngeom++;
				Surface.geomPushBack();
				Surface.setGeomParams(nelem, x.data(), r.data());
			}

	//		Surface.updateStorage(nelem, nlocl-1);
		}
	}

	if (counter > 0){
		// LAUNCH PUSHBACK FUNCTIONS FOR AS MANY TIMES AS NEEDED
		// THEN UPDATE SURFACE GEOMETRY HERE
		
		cout << counter << " node(s) added. Total of "
			<< nelem << " boundary elements." << endl;
	}
}

/* Redistribute nodes on the contour if two nodes 
 * are too close or too far apart. 
 * 
 * If a segment is too long, add a point in the middle.
 * If a segment is too short, remove the ith and (i+1)th
 * points and replace with a midpoint. */

 // WORK ON THIS AFTER YOU'RE DONE WITH checkAngle

void checkSpacing(surface &Surface, double lmin, double lmax, double thetmax){
	// declare variables
  int     i, j, i1, j1;
  int     nelem, ngeom, nlocl;
	int     counter;

	double  ax, bx, cx;
	double  ar, br, cr;

	double  xM, rM, lM; 

	double  li, li1;
	double  Dl, dl;


	// initialize
	counter = 0;

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
	
	/* initialize containers for node coordinates
	 * (these may be resized) */
	vector<double> x(ngeom), r(ngeom);

	// get geometric parameters
  Surface.getNode(x.data(), r.data());

	// loop over boundary elements
  for (i = 0; i < nelem; i++){
		i1 = i + 1;
	
		// get separation
		Surface.getArcl(i , li );
		Surface.getArcl(i1, li1);

		Dl = li1 - li;
		
		// check if segment exceeds maximum separation
		if (Dl > lmax){
			counter++;

			// resize containers
			nelem++;
			ngeom++;
			x.resize(ngeom);
			r.resize(ngeom);

			// get midpoint
			Surface.getSpln(i, ax, bx, cx,
			                   ar, br, cr);
			lM = 0.5*(li1 + li);
			dl = lM - li;
			xM = ((ax*dl + bx)*dl + cx)*dl + x[i];
			rM = ((ar*dl + br)*dl + cr)*dl + r[i];
		
			// add new point
			for(j = ngeom-1; j > i-1; j--){
				j1    = j+1;
				x[j1] = x[j];
				r[j1] = r[j];
			}
			x[i]   = xM;
			r[i]   = rM;
			
			// update geometry
//			Surface.updateStorage(nelem, nlocl-1);
			Surface.setGeomParams(nelem, x.data(), r.data());
		}
		
//		// check if segment falls below minimum separation
//		if (Dl < lmin){
//			counter--;
//
//			/* don't resize containers yet, but do 
//			 * update number of elements and nodes */
//			nelem--;
//			ngeom--;
//
//			// get midpoint
//			Surface.getSpln(i, ax, bx, cx,
//			                   ar, br, cr);
//			lM = 0.5*(li1 + li);
//			dl = lM - li;
//			xM = ((ax*dl + bx)*dl + cx)*dl + x[i];
//			rM = ((ar*dl + br)*dl + cr)*dl + r[i];
//
//			// replace two points with midpoint
//			for (j = i; j < ngeom+1; j++){
//				j1 = j+1;
//				x[j] = x[j1];
//				r[j] = r[j1];
//			}
//			x[i-1] = xM;
//			r[i-1] = rM;
//
//			// resize containers
//			x.resize(ngeom);
//			r.resize(ngeom);
//
//			// update geometry
//			Surface.updateStorage(nelem, nlocl-1);
//			Surface.setGeomParams(nelem, x.data(), r.data());
//		}
	}
	
	if (counter > 0)
		cout << counter << " node(s) added. Total of "
			<< nelem << " boundary elements." << endl;
	
	if (counter < 0)
		cout << -counter << " node(s) removed. Total of "
			<< nelem << " boundary elements." << endl;

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
