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
void checkAngle(surface, double, int &);
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
//	vector<double> nx (ngeom), nr (ngeom);
//	vector<double> Dfx(ngeom), Dfr(ngeom);
//	vector<double> vx (ngeom), vr (ngeom), vn(ngeom);
	vector<double> v  (2*nglob);


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
    Surface.setVel(i, 0.0, 0.0);
	//	vx[i    ] = 0.;
	//	vr[i    ] = 0.;
		v [2*i  ] = 0.;
		v [2*i+1] = 0.;
	}
	//Surface.setVel(nelem, nlocl-1, vx.data(), vr.data());
	
	// get initial geometry
	Surface.getNode(x .data(), r .data());
//	Surface.getNrml(nx.data(), nr.data());
	area = Surface.getArea();
	vlme = Surface.getVlme();

	// extrema for node redistribution
	thetmax = M_PI/8;
	// lmax = 
	// lmin =
	
	
  /*-----------------------------------------------------*/
  /*----------------- TIME INTEGRATION ------------------*/
  /*-----------------------------------------------------*/
	
	for (istep = 0; istep < nstep; istep++){
		cout << "ts = " << istep << ", V = " << vlme << endl;
		
		// write to file before evolving system
		writeNode(istep, ngeom, x.data(), r.data(), opath);
		
		/*-- Step 1: Solve the boundary integral equation. --*/
		singleLayer(IGF, nquad, Surface, v.data());

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

      // get velocity at the ith node
			// START BACK FROM HERE
			// NEED TO REMOVE THE VX, VR, VN VECTORS LOCALLY
			// MAKE USE OF THE STOKES.H CLASS, HAVE THE FIELDS
			// STORED THERE AND DON'T BOTHER WITH HAVING THE LOCAL FIELDS


			// calculate surface velocity
			vx  = v[2*n  ];
			vr  = v[2*n+1];
			vn  = vx*nx + vr*nr;
			
			// advect geometric nodes using forward Euler scheme
			x[i] += nx*vn*dt;
			r[i] += nr*vn*dt;
			
			// NOTE: NEED TO PROTECT AGAINST NEGATIVE RADIUS!!!!
			// MAYBE REDISTRIBUTE POINTS??
			
			if (r[i] < 0){
				for (j = 0; j < nelem; j++){
					printf("r[%d] = %.4f\n", j, r[j]);
				}
			}


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
		checkAngle(Surface, thetmax, istop);

		if (istop == 1)
			exit; // DON'T KNOW WHAT TO DO HERE
	
		// THEN CHECK ISTOP
		// IF ISTOP == 1, BREAK LOOP //

		/*-- Step 5: Get parameters for next timestep  ------*/
    // NOTE: NEED TO FIX THIS STUFF!!!
		x.resize(ngeom);
    r.resize(ngeom);
		Surface.getArea(area);
		Surface.getVlme(vlme);



		
	}
}

/*- NODE REDISTRIBUTION ---*/

/* Redistribute nodes on the contour if the angle 
 * subtended by an arc is too large.
 *
 * If i = 0, add another point. If i > 0, remove
 * the middle point and add two evenly spaced
 * points. */
void checkAngle  (surface Surface, double thetmax, int &istop){
	// declare variables
  int     i,  j,  j1;
	int     i0, i1;
  int     nelem, nnode;
	int     counter;

	double *s , *l ;
	double *ks, *kp;
	double *ax, *bx, *cx;
	double *ar, *br, *cr;

	double  axM, bxM, cxM;
	double  arM, brM, crM;
	double  ax0, bx0, cx0;
	double  ar0, br0, cr0;
	double  ax1, bx1, cx1;
	double  ar1, br1, cr1;

	double  xM, rM, lM;
	double  x0, r0, l0;
	double  x1, r1, l1;

	double  dlM, dl0, dl1;	
	double  Dthet;

	// constants
	const double PIH = 0.5*M_PI;

	// initialize
	istop = 0;
	counter = 0;

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  nnode = nelem + 1;
	
	/* initialize containers for node coordinates
	 * (these may be resized) */
	vector<double> x(nnode), r(nnode);

	// allocate memory
  s  = (double*) malloc( nnode    * sizeof(double));
  l  = (double*) malloc( nnode    * sizeof(double));
  ks = (double*) malloc( nnode    * sizeof(double));
  kp = (double*) malloc( nnode    * sizeof(double));
  ax = (double*) malloc((nnode-1) * sizeof(double));
  bx = (double*) malloc( nnode    * sizeof(double));
  cx = (double*) malloc((nnode-1) * sizeof(double));
  ar = (double*) malloc((nnode-1) * sizeof(double));
  br = (double*) malloc( nnode    * sizeof(double));
  cr = (double*) malloc((nnode-1) * sizeof(double));

	// get geometric parameters
  Surface.getNode(x.data(), r.data());
  Surface.getArcl(s);
  Surface.getPoly(l);
  Surface.getCurv(ks, kp);
	Surface.getSpln(ax, bx, cx,
	                ar, br, cr);

  for (i = 0; i < nelem; i++){
		i0 = i - 1;
		i1 = i + 1;

  	if (i == 0){
			Dthet = 2*s[1]*ks[0];
			Dthet = fabs(Dthet);
  	}
  	else {
			Dthet = (s[i+1] - s[i-1])*ks[i];
			Dthet = fabs(Dthet);
  	}

		if (Dthet > PIH){
			cout << "checkAngle: an arc is excessive at " <<
				(i+1) << "th index" << endl;
			
			istop = 1;
			return;
		}

		if (Dthet > thetmax){
			counter++;

			// resize containers
			nelem++;
			nnode++;
			x.resize(nnode);
			r.resize(nnode);

			if (i == 0){	// add a new point

				// get spline coefficients
				axM = ax[i]; bxM = bx[i]; cxM = cx[i];
				arM = ar[i]; brM = br[i]; crM = cr[i];

				// get midpoint
				 lM = (l[i] + l[i+1])/2.0;
				dlM = lM - l[i];
				 xM = ((axM*dlM + bxM)*dlM + cxM)*dlM + x[i];
				 rM = ((arM*dlM + brM)*dlM + crM)*dlM + r[i];

				// add point
				for (j = nnode; j > 1; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[1] = xM;
				r[1] = rM;

			}
			else {				/* remove the middle point and add
			               * two evenly spaced points */
				// get spline coefficients
				ax0 = ax[i0]; bx0 = bx[i0]; cx0 = cx[i0];
				ar0 = ar[i0]; br0 = br[i0]; cr0 = cr[i0];
				ax1 = ax[i ]; bx1 = bx[i ]; cx1 = cx[i ];
				ar1 = ar[i ]; br1 = br[i ]; cr1 = cr[i ];

				// compute new points
				 l0 = (l[i0] + 2.0*l[i])/3.0; 
				dl0 = l0 - l[i];
				 x0 = ((ax0*dl0 + bx0)*dl0 + cx0)*dl0 + x[i0];
				 r0 = ((ar0*dl0 + br0)*dl0 + cr0)*dl0 + r[i0];

				 l1 = (2.0*l[i] + l[i1])/3.0; 
				dl1 = l1 - l[i];
				 x1 = ((ax1*dl1 + bx1)*dl1 + cx1)*dl1 + x[i ];
				 r1 = ((ar1*dl1 + br1)*dl1 + cr1)*dl1 + r[i ];

				for(j = nnode; j > i; j--){
					j1    = j+1;
					x[j1] = x[j];
					r[j1] = r[j];
				}
				x[i]   = x1;
				r[i]   = r1;

				x[i-1] = x0;
				r[i-1] = r0;
			}

			// update surface geometry
			Surface.setGeomParams(nelem, x.data(), r.data());
		}
	}

	if (counter > 0)
		cout << counter << " nodes added. Total of "
			<< nelem << " boundary elements." << endl;

	// free memory
	free(s ); 
	free(l );
	free(ks);
	free(kp);
	free(ax);
	free(bx);
	free(cx);
	free(ar);
	free(br);
	free(cr);
}

/* Redistribute nodes on the contour if two nodes 
 * are too close or too far apart. */
void checkSpacing(surface Surface, double lmin, double lmax){
	// declare variables
	
	

}

/*- WRITE FUNCTIONS -------*/

/* Write geometric coordinates to file */
void writeNode(int istep, int nnode, double *x, double *r, string path){
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
	for (i = 0; i < nnode; i++){
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
