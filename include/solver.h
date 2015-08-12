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
	
//	double * v ;
//	double * vx, * vr, * vn;
//	double *Dfx, *Dfr;
//
//	double * x, * r;
//	double *nx, *nr;

	int    IGF;
	int    ISURF;
	
  /* get number of boundary elements,
   * geometric nodes, and basis nodes */
  nelem = Surface.getNElem();
  ngeom = nelem + 1;
  nlocl = Surface.getNLocl();
  nglob = Surface.getNGlob();

	// allocate memory (for performance)
//	v .reserve(2*nglob);
//	vx.reserve(  ngeom);
//	vr.reserve(  ngeom);
//	vn.reserve(  ngeom);
//	x .reserve(  ngeom);
//	r .reserve(  ngeom);
//	nx.reserve(  ngeom);
//	nr.reserve(  ngeom);
	
	// initialize vectors
	vector<double> v  (2*nglob);
	vector<double> vx (  ngeom), vr (ngeom), vn(ngeom);
	vector<double> Dfx(  ngeom), Dfr(ngeom);
	vector<double> x  (  ngeom), r  (ngeom);
	vector<double> nx (  ngeom), nr (ngeom);


//	// allocate memory
//	v   = (double*) malloc( 2*nglob * sizeof(double));
//	vx  = (double*) calloc(   ngeom , sizeof(double));
//	vr  = (double*) calloc(   ngeom , sizeof(double));
//	vn  = (double*) calloc(   ngeom , sizeof(double));
//	x   = (double*) malloc(   ngeom * sizeof(double));
//	r   = (double*) malloc(   ngeom * sizeof(double));
//	nx  = (double*) malloc(   ngeom * sizeof(double));
//	nr  = (double*) malloc(   ngeom * sizeof(double));


  /*-----------------------------------------------------*/
  /*-------------------- INITIALIZE ---------------------*/
  /*-----------------------------------------------------*/
	
	IGF = 1;   // tube Green's function
	IGF = 0;   // free-space Green's function

	ISURF = 2; // vesicle
	ISURF = 0; // drop

	// set surface velocity
	for (i = 0; i < nglob; i++){
		Surface.setVel(i, 0., 0.);
	}
	
	// get initial geometry
	Surface.getNode(x .data(), r .data());
	Surface.getNrml(nx.data(), nr.data());
	Surface.getArea(area);
	Surface.getVlme(vlme);
	
	
  /*-----------------------------------------------------*/
  /*----------------- TIME INTEGRATION ------------------*/
  /*-----------------------------------------------------*/
	
	for (istep = 0; istep < nstep; istep++){
		cout << "ts = " << istep << endl;
		
		/*-- Step 1: Solve the boundary integral equation. --*/
		singleLayer(IGF, nquad, Surface, v.data());

		// NOTE: FOR SINGLE LAYER POTENTIAL, THERE IS NO NEED
		// TO CALCULATE A MATRIX INVERSE


		// write to file before evolving system
		// ADD WRITE FUNCTION

	
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

			// NEED TO PROTECT AGAINST NEGATIVE RADIUS!!!!
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

		// CHECK VOLUME
		writeNode(istep, ngeom, x.data(), r.data(), opath);



		// NOTE: FOR NOW, DON'T RECALCULATE TENSION
		// AND MOMENTS, BECAUSE WE'RE USING A DROP.
		
		// NOTE #2: ALSO, WOULD HAVE TO UPDATE SURFACTANT
		// CONCENTRATION FIELD HERE.
		
		
//	// free memory	
//	free(v );
//	free(vx);
//	free(vr);
//	free(vn);
//	free(x );
//	free(r );
//	free(nx);
//	free(nr);
	}
}

//void writeNode(int n, double *x, double*r){
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
