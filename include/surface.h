/* SURFACE MECHANICS
 *  Boundary class that inherits the members of STOKES.H with type = 1
 *  (fluid-fluid interface) and, in addition, prescribes a surface
 *  constitutive model for the in-plane tensions (shear and dilatation),
 *  transverse shear tensions, and in-plane moments (bending and twisting).
 *
 * REFERENCES
 *  Mollman, John Wiley & Sons (1982)
 *  
 * PARAMETERS
 *  gamm  [input]			mean tension
 *  ES    [input]			shear modulus
 *  ED    [input]			dilatational modulus
 *  EB    [input]			bending modulus
 *  ET    [input]			twisting modulus
 *  tau   [output]		in-plane tensions
 *  q     [output]		transverse shear tension
 *  m			[output]		moments
 */

#ifndef SURFACE_H
#define SURFACE_H

/* HEADER FILES */
#include "stokes.h"

class surface: public stokes {
private:
  // indicator for constitive model
  int model;
  /*  = 0  (drop w/constant surfactant)
   *  = 1  (drop w/varying surface tension)
   *  = 2  (vesicle w/fixed surface area - Helfrich model)
   *  = 3  (vesicle w/variable surface area)
   *  = 4  (red blood cell - Skalak model) */
  
  // in-plane tensions
  vector<double> tenss, tensp;

  // transverse shear tensions
  vector<double> tensn;

  // moments
  vector<double> mmnts, mmntp;

	// elastic constants
	double shear, dilat;
	double bendg, twist;

	// mean tension
	double tensM;
 
 	/* Reserve memory for containers */
  void resizeContainers(int n, int m){
    int k = n*m + 1;

    stokes::resizeContainers(n,m);
        
    // reserve memroy
		if (k > tenss.size()){
			tenss.resize(k);
			tensp.resize(k);
			tensn.resize(k);
			mmnts.resize(k);
			mmntp.resize(k);
		}
  }
 
public:
  /* PROTOTYPES */
  
  
  /* IMPLEMENTATIONS */


	/*- CONSTRUCTORS ----*/
  surface() : stokes() {
  }

  surface(int id, int N, int M,
          double lamb, double gamm,
          double ES, double ED, double EB, double ET,
          double *x, double *r) : stokes(1, N, M, lamb, x, r) {
    // set constitutive model
    model = id;
		
		// resize containers
    tenss.resize(nglob);
    tensp.resize(nglob);
    tensn.resize(nglob);
    mmnts.resize(nglob);
    mmntp.resize(nglob);

		// set tensions and moments
		setSurfParams(N, M, gamm, ES, ED, EB, ET);
  }
	
	/*- DESTRUCTOR ------*/
	~surface(){
	}

	/*- UPDATE STORAGE --*/
	void updateStorage(int n, int m){
		resizeContainers(n, m);
	}
  
  /*- SET FUNCTIONS ---*/
  
	// set all internal tensions and moments
	void setSurfParams(int n,                 // element index
                     int m,                 // subelement index
                     double gamm,           // mean tension
                     double ES, double ED,  // moduli for tensions
                     double EB, double ET){ // moduli for moments
    // declare variables
    int i;

		// reserve memory
		resizeContainers(n, m);
		
		// set elastic constants
		shear = ES;
		dilat = ED;
		bendg = EB;
		twist = ET;

		// set mean tension
		tensM = gamm;

    if (model == 0){ // drop w/constant surface tension
      for (i = 0; i < nglob; i++){
        tenss[i] = tensM;
        tensp[i] = tensM;
        tensn[i] = 0;
        mmnts[i] = 0;
        mmnts[i] = 0;
      }
    }
    
    if (model == 1){ // drop w/varying surface tension
      for (i = 0; i < nglob; i++){
         
        // NEED TO FINISH THIS
    
      }
    }
    
    if (model == 2){ // vesicle w/constant surface area
      for (i = 0; i < nglob; i++){
        
        // NEED TO FINISH THIS

      }
    }
    
    // ADD IF STATMENTS FOR id == 3 (vesicle w/non-constant area)
    // and id == 4 (red blood cell)
	}

  // set parameters for a drop
  void setSurfParams(int n, int m){
		if (model == 0)
	    setSurfParams(n, m, tensM, 0.0, 0.0, 0.0, 0.0);
  }

  void setSurfParams(int n, int m, double gamm){
    if (model != 0 && model != 1)
      printf("Error: not enough parameters given for the surface constitutive model\n");

    setSurfParams(n, m, gamm, 0.0, 0.0, 0.0, 0.0);
  }

  /*- GET FUNCTIONS ---*/
  
	// get mean tension
	double getMeanTens(){
		return(tensM);
	}

  void getMeanTens(double &gamm){
		gamm = tensM;
  }
  
  // get tension components at the (iglob)th global basis node
  void getTens(int iglob, double &taus, double &taup, double &q){
    if (iglob >= nglob){
      printf("Error: index out of bounds in taus, taup, q.\n");
      return;
    }

    taus = tenss[iglob];
    taup = tensp[iglob];
    q    = tensn[iglob];
  }

  // get tension components at all global basis nodes
  void getTens(double *taus, double *taup, double *q){
    int iglob;

    if (taus == NULL || taup == NULL || q == NULL){
      printf("Error: no memory allocated for taus, taup, q.\n");
    }

    for (iglob = 0; iglob < nglob; iglob++){
      taus[iglob] = tenss[iglob];
      taup[iglob] = tensp[iglob];
      q   [iglob] = tensn[iglob];
    }
  }

  /* Function to calculate the force resultants normal and tangential
   * (in the meridional direction) to a free surface with constitutive
   * tensions and moments. The force resultants are evaluated at the
	 * basis nodes (collocation points) of the surface. 
	 *
	 * The function also returns the surface geometric parameters
	 * evaluated at the collocation points: the principal curvatures,
	 * tangent vector, and normal vector. */
  void calcForce(double *fx, double *fr,
	               double *ks, double *kp,
	               double *tx, double *tr,
	               double *nx, double *nr){
		// declare variables
		int     i ,  j ,  n ;
		double *x , *r ;
		double *fs, *fn;
		double  ax,  bx,  cx;
		double  ar,  br,  cr;
		double  x0,  x1;
		double  r0,  r1;
		double  l0,  l1;
		double  lM,  lD;

		double  dxdl  , drdl  , dsdl;
		double  d2xdl2, d2rdl2;

		double *zlocl, *L, *dLdx;
		double  z , l,  dl;
		double  cf;

		// allocate memory
		x     = (double*) malloc(nglob * sizeof(double));
		r     = (double*) malloc(nglob * sizeof(double));
		fs    = (double*) malloc(nglob * sizeof(double));
		fn    = (double*) malloc(nglob * sizeof(double));
		zlocl = (double*) malloc(nlocl * sizeof(double));

		/* calculate Gauss-Lobatto points on the
		 * interval [-1,1] */
		cf = M_PI/(nlocl-1);
		for (i = 0; i < nlocl; i++){
			zlocl[nlocl-i-1] = gsl_sf_cos(cf*i);
		}
		// NOTE: THIS IS A REDUNDANT CALCULATION - MAYBE INCORPORATE INTO STOKES CLASS
		// AS A MEMBER?? ONCE NLOCL IS KNOWN, IT IS STRAIGHTFORWARD TO
		// CALCULATE ZLOCL ON THE INTERVAL [-1,1]... THE SAME PROBLEM EXISTS IN QUAD.H
		// WHERE WE NEED TO MAKE USE OF LISTS



		/*----------------------------------------------*/
		/*---------------     setup     ----------------*/
		/*----------------------------------------------*/
		
		/* interpolate to local element nodes and
		 * calculate curvature, tangent, and normal */
		for (i = 0; i < nelem; i++){		// loop over boundary elements
			// get nodal coordinates
			x0 = nodex[i  ];
			x1 = nodex[i+1];
			r0 = noder[i  ];
			r1 = noder[i+1];

			// get spline coefficients
			ax = splnax[i];
			bx = splnbx[i];
			cx = splncx[i];
			ar = splnar[i];
			br = splnbr[i];
			cr = splncr[i];

			// get polygonal arc length
			l0 = poly[i  ];
			l1 = poly[i+1];
			lM = 0.5*(l1 + l0);
			lD = 0.5*(l1 - l0);
			
			for (j = 0; j < nlocl; j++){	// loop over local basis nodes
				// get global index for basis node
				n = i*(nlocl - 1) + j;
				
				if (j == 0){							  /* basis node coincides with
				                             * lower geometric node */
					x [n] = x0      ; 
					r [n] = r0      ;
					ks[n] = curvs[i];
					kp[n] = curvp[i];
					tx[n] = tangx[i];
					tr[n] = tangr[i];
					nx[n] = nrmlx[i];
					nr[n] = nrmlr[i];
				}
				else if (j == nlocl - 1){		/* basis node coincides with
				                             * upper geometric node */
					x [n] = x1        ;
					r [n] = r1        ;
					ks[n] = curvs[i+1];
					kp[n] = curvp[i+1];
					tx[n] = tangx[i+1];
					tr[n] = tangr[i+1];
					nx[n] = nrmlx[i+1];
					nr[n] = nrmlr[i+1];

					/* NOTE: This assignment is redundant for
					 * intermediate nodes. */
				}
				else {
					// map Gauss-Lobatto point onto polygonal arc length
					l  = lM + lD*zlocl[j];
					dl = l  - l0;

					// interpolate to basis node
					x [n]  = ((  ax*dl +    bx)*dl + cx)*dl + x0;
					r [n]  = ((  ar*dl +    br)*dl + cr)*dl + r0;

					// calculate derivatives
					dxdl   = (3.*ax*dl + 2.*bx)*dl + cx ;
					drdl   = (3.*ar*dl + 2.*br)*dl + cr ;
					d2xdl2 =  6.*ax*dl + 2.*bx          ;
					d2rdl2 =  6.*ar*dl + 2.*br          ;
					dsdl   = sqrt(dxdl*dxdl + drdl*drdl);
					
					// calculate curvatures
					ks[n]  = -(dxdl*d2rdl2 - d2xdl2*drdl)/pow(dsdl, 3);
					kp[n]  =   dxdl/(r[n]*dsdl);

					// calculate tangent components
					tx[n] = dxdl/dsdl;
					tr[n] = drdl/dsdl;

					// calculate normal components
					nx[n] =  tr[n];
					nr[n] = -tx[n];
				}
			}
		}
		

		// WILL NEED TO CALCULATE THE LAGRANGE POLYNOMIALS
		// AAAAND THEIR DERIVATIVES FOR THE SPATIAL DERIVATIVE TERMS!
		// (NOT FOR THE DROP, BUT FOR A BOUSSINESQ SURFACE FOR
		// SURE).
		//lagrange(nlocl-1, nquad, zlocl, zquad);

		/* NOTE: Calculating the Lagrange polynomials here
		 * is redundant; that is, they are calculated else-
		 * where for the same interval. Could improve per-
		 * formance by eliminating redundancy, but since
		 * this is a small calculation, we'll just keep it
		 * here for now */


		
		for (i = 0; i < nelem; i++){		// loop over boundary elements
			for (j = 0; j < nlocl; j++){	// loop over local basis nodes
				// global element node
				n = i*(nlocl - 1) + j;
				
				// FOR NOW, JUST FOCUS ON THE DROP...
				if (model == 0){
					// surface tension force
					fs[n] = 0;
					fn[n] = -(ks[n]*tenss[n] + kp[n]*tensp[n]);

					// add buoyance force (assume drho = 1, acceleration due to gravity = 1)
					double rhod = 1.;
					double rhoa = 0.;
					double gac  = 1.;
					fn[n] -= (rhod - rhoa)*gac*x[n];
					
					// calculate x and r components
					fx[n] = fs[n]*tx[n] + fn[n]*nx[n];
					fr[n] = fs[n]*tr[n] + fn[n]*nr[n];

				}
				
				// ...AND THEN LATER, THE VESICLE
			}
		}

		// release memory
		free(x    );
		free(r    );
		free(fs   );
		free(fn   );
		free(zlocl);
  }

	/* Same as the previous function, but only returns the force components. */
  void calcForce(double *fx, double *fr){
		// declare variables
		double *ks, *kp;
		double *tx, *tr;
		double *nx, *nr;
		
		// allocate memory
		ks = (double*) malloc(nglob * sizeof(double));
		kp = (double*) malloc(nglob * sizeof(double));
		tx = (double*) malloc(nglob * sizeof(double));
		tr = (double*) malloc(nglob * sizeof(double));
		nx = (double*) malloc(nglob * sizeof(double));
		nr = (double*) malloc(nglob * sizeof(double));

		// calculate forces
		calcForce(fx, fr, ks, kp, tx, tr, nx, nr);

		// release memory
		free(ks);
		free(kp);
		free(tx);
		free(tr);
		free(nx);
		free(nr);
	}
};

#endif
