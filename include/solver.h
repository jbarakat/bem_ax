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
void writeNode(int, int, double*, double*, string);

/* IMPLEMENTATIONS */

/* Time integration */
void timeInt(int nstep, int nquad, double dt, surface Surface, string opath){

  /*-----------------------------------------------------*/
  /*----------------------- SETUP -----------------------*/
  /*-----------------------------------------------------*/

	// declare variables
	int    istep;
	int    i, j, k, m, n;
	int    nelem, ngeom;
	int    nlocl, nglob;

	double area, vlme;
	
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
	vector<double> nx (ngeom), nr (ngeom);
	vector<double> Dfx(ngeom), Dfr(ngeom);
	vector<double> vx (ngeom), vr (ngeom), vn(ngeom);
	vector<double> v  (2*nglob);


  /*-----------------------------------------------------*/
  /*-------------------- INITIALIZE ---------------------*/
  /*-----------------------------------------------------*/

	IGF = 1;   // tube Green's function
	IGF = 0;   // free-space Green's function

	ISURF = 2; // vesicle
	ISURF = 0; // drop

	// set surface velocity
	for (i = 0; i < nglob; i++){
		vx[i    ] = 0.;
		vr[i    ] = 0.;
		v [2*i  ] = 0.;
		v [2*i+1] = 0.;
	}
	Surface.setVel(nelem, nlocl-1, vx.data(), vr.data());
	
	// get initial geometry
	Surface.getNode(x .data(), r .data());
	Surface.getNrml(nx.data(), nr.data());
	Surface.getArea(area);
	Surface.getVlme(vlme);
	
	
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
	
			// calculate surface velocity
			vx[i]  = v[2*n  ];
			vr[i]  = v[2*n+1];
			vn[i]  = vx[i]*nx[i] + vr[i]*nr[i];
			
			// advect geometric nodes using forward Euler scheme
			x [i] += nx[i]*vn[i]*dt;
			r [i] += nr[i]*vn[i]*dt;
			
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
		
		Surface.getNrml(nx.data(), nr.data());
		Surface.getArea(area);
		Surface.getVlme(vlme);




		// NOTE: FOR NOW, DON'T RECALCULATE TENSION
		// AND MOMENTS, BECAUSE WE'RE USING A DROP.
		
		// NOTE #2: ALSO, WOULD HAVE TO UPDATE SURFACTANT
		// CONCENTRATION FIELD HERE.
		
	}
}

/*- NODE REDISTRIBUTION ---*/

/* Redistribute nodes on the contour if the angle 
 * subtended by an arc is too large. */
void checkAngle  (surface Surface, double thetmax){
	// declare variables
  int     i;
  int     nelem, nnode;
	double *s, *ks, *kp;

  // get number of geometric elements and nodes
  nelem = Surface.getNElem();
  nnode = nelem + 1;

	// allocate memory
  s  = (double*) malloc(nnode * sizeof(double));
  ks = (double*) malloc(nnode * sizeof(double));
  kp = (double*) malloc(nnode * sizeof(double));

	// get meridional arc length and principal curvatures
  Surface.getArcl(s);
  Surface.getCurv(ks, kp);

  for (i = 0; i < nnode; i++){
    if (i == 0){
       
    }
    else {

    }
  }
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
