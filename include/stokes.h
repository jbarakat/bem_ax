/* STOKES FLOW
 *  Boundary class that inherits the members of GEOM.H and, in addition,
 *  prescribes the fluid properties (e.g., viscosity) as well as the
 *  kinematic (e.g., velocity) and dynamic (e.g., stress) quantities 
 *  associated with the boundary in a Stokes flow.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  x,r    [input]			nodal coordinates
 *  N      [input]			number of (global) elements
 *  M      [input]			number of (local) subelements
 *  lamb   [input]			viscosity ratio
 *  u      [output]			displacement
 *  v      [output]			velocity
 *  f      [output]			traction
 *  c		   [output]			concentration (or any scalar field)
 */

#ifndef STOKES_H
#define STOKES_H

/* HEADER FILES */
#include "geom.h"

const int RIGID = 0;
const int FLUID = 1;

class stokes: public geom {
friend class surface;
private:
	// indicator for type of boundary
	int            type ;
	/*              = RIGID  (rigid boundary)
	 *              = FLUID  (fluid-fluid interface) */

	/* local and global number of basis nodes 
	 * for density function interpolation */
	int            nlocl, nglob;
	
	// displacement
	vector<double> dispx, dispr;

	// velocity
	vector<double> velx , velr ;

	// traction
	vector<double> trctx, trctr;

	// concentration
	vector<double> conc ;

	// viscosity ratio
	double         visc ;
	/* NOTE: By convention, the denominator
	 *       corresponds to the phase into
	 *       which the normal vector points. */

	/* Reserve enough space for containers */
	void resizeContainers(int n, int m){
    int k = n*m + 1;

		// update geometric parameters
		geom::resizeContainers(n);
		
		// reserve memory
		if (k > dispx.size()){
			dispx.resize(k);
			dispr.resize(k);
			velx .resize(k);
			velr .resize(k);
			trctx.resize(k);
			trctr.resize(k);
			conc .resize(k);
		}
		
		/* set number of local and global
		 * basis nodes */
		nlocl = m + 1;
		nglob = k;
		
	}

public:

	/* PROTOTYPES */
	
	
	/* IMPLEMENTATIONS */



	/*- CONSTRUCTORS-----*/
	stokes() : geom() {
	}

	stokes(int id, int N, int M, 
	       double *x, double *r) : geom(N, x, r) {
		// set viscosity ratio = 1 if not given
		stokes(id, N, M, 1., x, r);
	}

	stokes(int id, int N, int M, 
	       double lamb, double *x, double *r) : geom(N, x, r) {
		// set number of local and global basis nodes
		nlocl = M+1;
		nglob = N*M+1;

		// error flags
		if (id != 0 && id != 1){
			printf("Error: id can only take values of 0 or 1.\n");
			return;
		}
		
		if (M < 1){
			printf("Error: cannot have fewer than 1 subelements.\n");
			return;
		}

		int i, n;

		// set boundary type
		type = id;

		// set viscosity ratio
		visc = lamb;
		
		// resize containers
		dispx.resize(nglob);
		dispr.resize(nglob);
		velx .resize(nglob);
		velr .resize(nglob);
		trctx.resize(nglob);
		trctr.resize(nglob);
		conc .resize(nglob);
	}

	/*- DESTRUCTOR ------*/
	~stokes(){
	}

	/*- PUSH/POP BACK ---*/

	void stksPushBack(){
		dispx.push_back(0);
		dispr.push_back(0);
		velx .push_back(0);
		velr .push_back(0);
		trctx.push_back(0);
		trctr.push_back(0);
		conc .push_back(0);
	}

	void stksPushBack(double n){
		dispx.push_back(n);
		dispr.push_back(n);
		velx .push_back(n);
		velr .push_back(n);
		trctx.push_back(n);
		trctr.push_back(n);
		conc .push_back(n);
	}

	void stksPopBack(){
		dispx.pop_back();
		dispr.pop_back();
		velx .pop_back();
		velr .pop_back();
		trctx.pop_back();
		trctr.pop_back();
		conc .pop_back();
	}


	/*- SET FUNCTIONS ---*/
	
	// set number of local basis nodes
	void setNLocl(int n){
		nlocl = n;
	}
	
	// set number of global basis nodes
	void setNGlob(int n){
		nglob = n;
	}

  // set all transport fields (density functions)
  void setStksParams(int n,                   // element index
                     int m,                   // subelement index
                     double *ux, double *ur,  // displacement
                     double *vx, double *vr,  // velocity
                     double *fx, double *fr,  // traction
                     double *c ){             // concentration
   // // reserve memory
   // resizeContainers(n, m);

    int i;

    for (i = 0; i < nglob; i++){
      dispx[i] = ux[i];
      dispr[i] = ur[i];
      velx [i] = vx[i];
      velr [i] = vr[i];
      trctx[i] = fx[i];
      trctr[i] = fr[i];
      conc [i] = c [i];
    }
  }
	
	/* set velocity at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void setVel(int ielem, int ilocl, double vx, double vr){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in vx, vr.\n");
			return;
		}
	  
		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;
		
		velx[iglob] = vx;
		velr[iglob] = vr;
	}

	// set velocity at the (iglob)th global basis node
	void setVel(int iglob, double vx, double vr){
		if (iglob >= nglob){
			printf("Error: index out of bounds in vx, vr.\n");
			return;
		}
		
		velx[iglob] = vx;
		velr[iglob] = vr;
	}

	// set velocity at all global basis nodes
	void setVel(int n, int m, double *vx, double *vr){
		resizeContainers(n, m);
		
		int iglob;

		if (vx == NULL || vr == NULL){
			printf("Error: no memory allocated for vx, vr.\n");
			return;
		}

		for (iglob = 0; iglob < nglob; iglob++){
			 velx[iglob] = vx[iglob];
			 velr[iglob] = vr[iglob];
		}
	}
	
	/* set traction at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void setTrct(int ielem, int ilocl, double fx, double fr){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in fx, fr.\n");
			return;
		}
	  
		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;
		
		trctx[iglob] = fx;
		trctr[iglob] = fr;
	}

	// set traction at the (iglob)th global basis node
	void setTrct(int iglob, double fx, double fr){
		if (iglob >= nglob){
			printf("Error: index out of bounds in fx, fr.\n");
			return;
		}
		
		trctx[iglob] = fx;
		trctr[iglob] = fr;
	}

	// set traction at all global basis nodes
	void setTrct(int n, int m, double *fx, double *fr){
		resizeContainers(n, m);
		
		int iglob;

		if (fx == NULL || fr == NULL){
			printf("Error: no memory allocated for fx, fr.\n");
			return;
		}

		for (iglob = 0; iglob < nglob; iglob++){
			 trctx[iglob] = fx[iglob];
			 trctr[iglob] = fr[iglob];
		}
	}

	/*- GET FUNCTIONS ---*/
  
	// get local number of basis nodes
	int getNLocl(){
		int n = nlocl;
		return(n);
	}
	
	void getNLocl(int n){
		n = nlocl;
	}

	// get global number of basis nodes
	int getNGlob(){
		int n = nglob;
		return(n);
	}
	
	void getNGlob(int n){
		n = nglob;
	}
	
	/* get displacement at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void getDisp(int ielem, int ilocl, double &ux, double &ur){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in ux, ur.\n");
			return;
		}

		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;

		ux = dispx[iglob];
		ur = dispr[iglob];
	}

	// get displacement at the (iglob)th global basis node
	void getDisp(int iglob, double &ux, double &ur){
		if (iglob >= nglob){
			printf("Error: index out of bounds in ux, ur.\n");
			return;
		}
		
		ux = dispx[iglob];
		ur = dispr[iglob];
	}

	// get displacement at all global basis nodes
	void getDisp(double *ux, double *ur){
		int iglob;

		if (ux == NULL || ur == NULL){
			printf("Error: no memory allocated for ux, ur.\n");
		}

		for (iglob = 0; iglob < nglob; iglob++){
			ux[iglob] = dispx[iglob];
			ur[iglob] = dispr[iglob];
		}
	}
	
	/* get velocity at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void getVel(int ielem, int ilocl, double &vx, double &vr){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in vx, vr.\n");
			return;
		}

		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;

		vx = velx[iglob];
		vr = velr[iglob];
	}

	// get velocity at the (iglob)th global basis node
	void getVel(int iglob, double &vx, double &vr){
		if (iglob >= nglob){
			printf("Error: index out of bounds in vx, vr.\n");
			return;
		}
		
		vx = velx[iglob];
		vr = velr[iglob];
	}

	// get velocity at all global basis nodes
	void getVel(double *vx, double *vr){
		int iglob;

		if (vx == NULL || vr == NULL){
			printf("Error: no memory allocated for vx, vr.\n");
		}

		for (iglob = 0; iglob < nglob; iglob++){
			vx[iglob] = velx[iglob];
			vr[iglob] = velr[iglob];
		}
	}

	/* get traction at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void getTrct(int ielem, int ilocl, double &fx, double &fr){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in fx, fr.\n");
			return;
		}

		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;

		fx = trctx[iglob];
		fr = trctr[iglob];
	}

	// get traction at the (iglob)th global basis node
	void getTrct(int iglob, double &fx, double &fr){
		if (iglob >= nglob){
			printf("Error: index out of bounds in fx, fr.\n");
			return;
		}
		
		fx = trctx[iglob];
		fr = trctr[iglob];
	}

	// get traction at all global basis nodes
	void getTrct(double *fx, double *fr){
		int iglob;

		if (fx == NULL || fr == NULL){
			printf("Error: no memory allocated for fx, fr.\n");
		}

		for (iglob = 0; iglob < nglob; iglob++){
			fx[iglob] = trctx[iglob];
			fr[iglob] = trctr[iglob];
		}
	}

	/* get concentration at the (ielem)th local basis node
	 * of the (ilocl)th boundary element */
	void getConc(int ielem, int ilocl, double &c){
		if ((ielem + 1)*ilocl >= nglob){
			printf("Error: index out of bounds in c.\n");
			return;
		}

		int iglob;
		
		// get index of the global basis node
		iglob = ielem*(nlocl - 1) + ilocl;

		c = conc[iglob];
	}

	// get concentration at the (iglob)th global basis node
	void getTrct(int iglob, double &c){
		if (iglob >= nglob){
			printf("Error: index out of bounds in c.\n");
			return;
		}
		
		c = conc[iglob];
	}

	// get concentration at all global basis nodes
	void getConc(double *c){
		int iglob;

		if (c == NULL){
			printf("Error: no memory allocated for fx, fr.\n");
		}

		for (iglob = 0; iglob < nglob; iglob++){
			c[iglob] = conc[iglob];
		}
	}

	// get viscosity ratio
	double getVisc(){
		double lamb;
		lamb = visc;
		return(lamb);
	}
	
	void getVisc(double &lamb){
		lamb = visc;
	}
};

#endif
