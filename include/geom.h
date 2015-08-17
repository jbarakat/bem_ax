/* GEOMETRY
 *  Boundary class that prescribes the spline elements and geometrical
 *  parameters of an axisymmetric contour with N+1 nodes.
 *
 * REFERENCES
 *  N/A
 *  
 * PARAMETERS
 *  x,r    [input]		nodal coordinates
 *  N	     [input]		number of boundary elements
 *  a,b,c  [output]		cubic spline coefficients (for interpolation)
 *  l      [output]		polygonal arc length (for interpolation)
 *  s      [output]		meridional arc length
 *  A      [output]		total area
 *  V      [output]		total volume
 *  ks,kp  [output]		principal curvatures
 *  t      [output]		meridional tangent vector
 *  n      [output] 	normal vector
 */

#ifndef GEOM_H
#define GEOM_H

/* HEADER FILES */
#include <vector>
#include <math.h>
#include "interp.h"

using namespace std;

class geom {
friend class stokes;
friend class surface;
private:
	// number of nodes and elements
	int            nnode , nelem  ;

	// nodal coordinates
	vector<double> nodex , noder ;

	// cubic spline coefficients
	vector<double> splnax, splnbx, splncx;
	vector<double> splnar, splnbr, splncr;

	// meridional and polygonal arc lengths
	vector<double> arcl  , poly  ;
	
	// area and volume
	double         area  , vlme  ;

	// curvatures, tangent and normal vectors
	vector<double> curvs , curvp ;
	vector<double> tangx , tangr ;
	vector<double> nrmlx , nrmlr ;
	
	/* Reserve enough memory for containers */
	void resizeContainers(int n){
    int k = n + 1;
		
		// reserve memory
		if (k > nodex.size()){
			nodex .resize(k  );
			noder .resize(k  );
			splnax.resize(k-1);
			splnbx.resize(k  );
			splncx.resize(k-1);
			splnar.resize(k-1);
			splnbr.resize(k  );
			splncr.resize(k-1);
			poly  .resize(k  );
			arcl  .resize(k  );
			nodex .resize(k  );
			noder .resize(k  );
			curvs .resize(k  );
			curvp .resize(k  );
			tangx .resize(k  );
			tangr .resize(k  );
			nrmlx .resize(k  );
			nrmlr .resize(k  );
		}

		// set number of nodes and elements
		nelem = n;
		nnode = k;
	}


public:
	/* PROTOTYPES */
	
	/* IMPLEMENTATIONS */

	/*- CONSTRUCTORS ----*/
	geom(){
	}

	geom(int N, double *x, double *r){
		// set number of geometric elements and nodes
		nelem = N  ;
		nnode = N+1;

		// error flags
		if (N < 1){
			printf("Error: cannot have fewer than 1 elements.\n");
			return;
		}
		
		// resize containers
		nodex .resize(nnode  );
		noder .resize(nnode  );
		splnax.resize(nnode-1);
		splnbx.resize(nnode  );
		splncx.resize(nnode-1);
		splnar.resize(nnode-1);
		splnbr.resize(nnode  );
		splncr.resize(nnode-1);
		poly  .resize(nnode  );
		arcl  .resize(nnode  );
		nodex .resize(nnode  );
		noder .resize(nnode  );
		curvs .resize(nnode  );
		curvp .resize(nnode  );
		tangx .resize(nnode  );
		tangr .resize(nnode  );
		nrmlx .resize(nnode  );
		nrmlr .resize(nnode  );
		
		// set geometric parameters
		setGeomParams(N, x, r);
	}

	/*- DESTRUCTOR ------*/
	~geom(){
	}
	
	/*- PUSH/POP BACK ---*/

	void geomPushBack(){
		nodex .push_back(0);
		noder .push_back(0);
		splnax.push_back(0);
		splnbx.push_back(0);
		splncx.push_back(0);
		splnar.push_back(0);
		splnbr.push_back(0);
		splncr.push_back(0);
		poly  .push_back(0);
		arcl  .push_back(0);
		nodex .push_back(0);
		noder .push_back(0);
		curvs .push_back(0);
		curvp .push_back(0);
		tangx .push_back(0);
		tangr .push_back(0);
		nrmlx .push_back(0);
		nrmlr .push_back(0);
	}
	
	void geomPushBack(double n){
		nodex .push_back(n);
		noder .push_back(n);
		splnax.push_back(n);
		splnbx.push_back(n);
		splncx.push_back(n);
		splnar.push_back(n);
		splnbr.push_back(n);
		splncr.push_back(n);
		poly  .push_back(n);
		arcl  .push_back(n);
		nodex .push_back(n);
		noder .push_back(n);
		curvs .push_back(n);
		curvp .push_back(n);
		tangx .push_back(n);
		tangr .push_back(n);
		nrmlx .push_back(n);
		nrmlr .push_back(n);
	}

	void geomPopBack(){
		nodex .pop_back();
		noder .pop_back();
		splnax.pop_back();
		splnbx.pop_back();
		splncx.pop_back();
		splnar.pop_back();
		splnbr.pop_back();
		splncr.pop_back();
		poly  .pop_back();
		arcl  .pop_back();
		nodex .pop_back();
		noder .pop_back();
		curvs .pop_back();
		curvp .pop_back();
		tangx .pop_back();
		tangr .pop_back();
		nrmlx .pop_back();
		nrmlr .pop_back();
	}
	
	/*- SET FUNCTIONS ---*/

	// calculate and set all quantities
	void setGeomParams(int n, double *x, double *r){
		int i;
		
//		// ensure enough space is allocated
//		resizeContainers(n);
		
		// set all parameters
		for (i = 0; i < n+1; i++){
			nodex[i] = x[i];
			noder[i] = r[i];
		}
	
		calcGeomParams(n            , x            , r            ,
		               splnax.data(), splnbx.data(), splncx.data(),
		               splnar.data(), splnbr.data(), splncr.data(),
		               poly  .data(), arcl  .data(),
				    			 area         , vlme         , 
				    			 curvs .data(), curvp .data(), 
		               tangx .data(), tangr .data(),
				    			 nrmlx .data(), nrmlr .data());
	}

	// set number of geometric nodes
	void setNNode(int n){
		nnode = n;
	}
	
	// set number of boundary elements
	void setNElem(int n){
		nelem = n;
	}

	/*- GET FUNCTIONS ---*/

	// get all quantities at all geometric nodes
	void getGeomParams(int    &n , double *x,  double *r,
	                   double *ax, double *bx, double *cx,
						         double *ar, double *br, double *cr,
						         double *l , double *s ,
						         double &A , double &V ,
						         double *ks, double *kp,
						         double *tx, double *tr,
						         double *nx, double *nr){
		int i;

		if (x  == NULL || r  == NULL ||
		    ax == NULL || bx == NULL || cx == NULL ||
		    ar == NULL || br == NULL || cr == NULL ||
				l  == NULL || s  == NULL || 
				ks == NULL || kp == NULL ||
				tx == NULL || tr == NULL ||
				nx == NULL || nr == NULL){
			printf("Error: no memory allocated for pointers.\n");
			return;
		}
		
		n = nelem;
		A = area;
		V = vlme;

		for (i = 0; i < nnode-1; i++){
			x [i] = nodex [i];
			r [i] = noder [i];
			ax[i] = splnax[i];
			bx[i] = splnbx[i];
			cx[i] = splncx[i];
			ar[i] = splnar[i];
			br[i] = splnbr[i];
			cr[i] = splncr[i];
			l [i] = poly  [i];
			s [i] = arcl  [i];
			ks[i] = curvs [i];
			kp[i] = curvp [i];
			tx[i] = tangx [i];
			tr[i] = tangr [i];
			nx[i] = nrmlx [i];
			nr[i] = nrmlr [i];
		}

		i = nnode-1;
			x [i] = nodex [i];
			r [i] = noder [i];
			bx[i] = splnbx[i];
			br[i] = splnbr[i];
			l [i] = poly  [i];
			s [i] = arcl  [i];
			ks[i] = curvs [i];
			kp[i] = curvp [i];
			tx[i] = tangx [i];
			tr[i] = tangr [i];
			nx[i] = nrmlx [i];
			nr[i] = nrmlr [i];

	}
	
	// get all quantities at the ith geometric node
	void getGeomParams(int i, int &n, double &x , double &r,
	                   double &ax,    double &bx, double &cx,
		       				   double &ar,    double &br, double &cr,
		       				   double &l ,    double &s ,
		       				   double &A,     double &V ,
		       				   double &ks,    double &kp,
		       				   double &tx,    double &tr,
		       				   double &nx,    double &nr){
		if (i >= nnode){
			printf("Error: index out of bounds.\n");
			return;
		}
		
		n = nelem;
		A = area;
		V = vlme;

		x  = nodex [i];
		r  = noder [i];
		ax = splnax[i];
		bx = splnbx[i];
		cx = splncx[i];
		ar = splnar[i];
		br = splnbr[i];
		cr = splncr[i];
		l  = poly  [i];
		s  = arcl  [i];
		ks = curvs [i];
		kp = curvp [i];
		tx = tangx [i];
		tr = tangr [i];
		nx = nrmlx [i];
		nr = nrmlr [i];
	}
	
	// get number of geometric nodes
	int getNNode(){
		int n;
		n = nnode;
		return(n);
	}
	
	void getNNode(int &n){
		n = nnode;
	}
	
	// get number of boundary elements
	int getNElem(){
		int n;
		n = nelem;
		return(n);
	}
	
	void getNElem(int &n){
		n = nelem;
	}

	// get position coordinates at all geometric nodes
	void getNode(double *x, double *r){
		int i;
		
		if (x == NULL || r == NULL){
			printf("Error: no memory allocated for x, r.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			x[i] = nodex[i];
			r[i] = noder[i];
		}
	}
	
	// get position coordinates at the ith geometric node
	void getNode(int i, double &x, double &r){
		if (i >= nnode){
			printf("Error: index out of bounds in x, r.\n");
			return;
		}
		
		x = nodex[i];
		r = noder[i];
	}

	// get spline coefficients at all geometric nodes
	void getSpln(double *ax, double *bx, double *cx,
	             double *ar, double *br, double *cr){
		int i;
		
		if (ax == NULL || bx == NULL || cx == NULL ||
		    ar == NULL || br == NULL || cr == NULL){
			printf("Error: no memory allocated for a, b, c.\n");
			return;
		}

		for (i = 0; i < nnode-1; i++){
			ax[i] = splnax[i];
			bx[i] = splnbx[i];
			cx[i] = splncx[i];
			ar[i] = splnar[i];
			br[i] = splnbr[i];
			cr[i] = splncr[i];
		}
		
		i = nnode-1;
			bx[i] = splnbx[i];
			br[i] = splnbr[i];
	}
	
	// get spline coefficients at the ith geometric node
	void getSpln(int i, 
	             double &ax, double &bx, double &cx,
					 		 double &ar, double &br, double &cr){
		if (i >= nnode){
			printf("Error: index out of bounds in a, b, c.\n");
			return;
		}
		
		ax = splnax[i];
		bx = splnbx[i];
		cx = splncx[i];
		ar = splnar[i];
		br = splnbr[i];
		cr = splncr[i];
	}

	// get polygonal arc length at all geometric nodes
	void getPoly(double *l){
		int i;
		
		if (l == NULL){
			printf("Error: no memory allocated for l.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			l[i] = poly[i];
		}
	}
	
	// get polygonal arc length at the ith geometric node
	void getPoly(int i, double &l){
		if (i >= nnode){
			printf("Error: index out of bounds in l.\n");
			return;
		}
		
		l = poly[i];
	}

	// get meridional arc length at all geometric nodes
	void getArcl(double *s){
		int i;
		
		if (s == NULL){
			printf("Error: no memory allocated for s.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			s[i] = arcl[i];
		}
	}
	
	// get meridional arc length at the ith geometric node
	void getArcl(int i, double &s){
		if (i >= nnode){
			printf("Error: index out of bounds in s.\n");
			return;
		}
		
		s = arcl[i];
	}
	
	// get total area
	double getArea(){
		double A;
		A = area;
		return(A);
	}
	
	void getArea(double &A){
		A = area;
	}
	
	// get total volume
	double getVlme(){
		double V;
		V = vlme;
		return(V);
	}

	void getVlme(double &V){
		V = vlme;
	}
	
	// get principal curvatures at all geometric nodes
	void getCurv(double *ks, double *kp){
		int i;
		
		if (ks == NULL || kp == NULL){
			printf("Error: no memory allocated for ks, kp.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			ks[i] = curvs[i];
			kp[i] = curvp[i];
		}
	}
	
	// get principal curvatures at the ith geometric node
	void getCurv(int i, double &ks, double &kp){
		if (i >= nnode){
			printf("Error: index out of bounds in ks, kp.\n");
			return;
		}
		
		ks = curvs[i];
		kp = curvp[i];
	}

	// get meridional tangent components at all geometric nodes
	void getTang(double *tx, double *tr){
		int i;
		
		if (tx == NULL || tr == NULL){
			printf("Error: no memory allocated for tx, tr.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			tx[i] = tangx[i];
			tr[i] = tangr[i];
		}
	}
	
	// get meridional tangent components at the ith geometric node
	void getTang(int i, double &tx, double &tr){
		if (i >= nnode){
			printf("Error: index out of bounds in tx, tr.\n");
			return;
		}
		
		tx = tangx[i];
		tr = tangr[i];
	}

	// get normal components at all geometric nodes
	void getNrml(double *nx, double *nr){
		int i;
		
		if (nx == NULL || nr == NULL){
			printf("Error: no memory allocated for nx, nr.\n");
			return;
		}

		for (i = 0; i < nnode; i++){
			nx[i] = nrmlx[i];
			nr[i] = nrmlr[i];
		}
	}
	
	// get normal components at the ith geometric node
	void getNrml(int i, double &nx, double &nr){
		if (i >= nnode){
			printf("Error: index out of bounds in nx, nr.\n");
			return;
		}
		
		nx = nrmlx[i];
		nr = nrmlr[i];
	}

	/* Function to calculate geometric parameters for a given set of N+1
	 * nodal coordinates (x,r) */
	void calcGeomParams(int      N, double *x,  double * r,	/* <--- inputs  */
	     								double *ax, double *bx, double *cx,	/* <--- outputs */
	     								double *ar, double *br, double *cr,	/* <------|     */
	     								double * l, double * s,             /* <------|     */
	     								double & A, double & V,             /* <------|     */
	     								double *ks, double *kp,             /* <------|     */
	     								double *tx, double *tr,							/* <------|     */
	     								double *nx, double *nr ){						/* <------|     */
		// declare variables
		int i,j, n;
		int na, nt;
		double axi, bxi, cxi;
		double ari, bri, cri;
		double xi, ri, xj, rj, lj;
		double dx, dr, ds, dA, dV, dl;
		double dxdl, drdl, dsdl, dAdl, dVdl;
		double d2xdl2, d2rdl2;
		double ssum, Asum, Vsum;
		const int MAXIT = 20;
		double stol, Atol, Vtol;
		double fc;

		// calculate polygonal arc length
		l[0] = 0;
		for (i = 1; i < N+1; i++){
			dx = x[i] - x[i-1];
			dr = r[i] - r[i-1];
			dl = sqrt(dx*dx + dr*dr);
			l[i] = l[i-1] + dl;
		}

		// set integration tolerances
		dl = l[N] - l[0];
		stol = 0.0000001*dl;
		Atol = 0.0000001*dl*dl;
		Vtol = 0.0000001*dl*dl*dl;
		
		// calculate cubic spline coefficients for interpolation
		spline(N, l, x, 0.,  0., ax, bx, cx);
		spline(N, l, r, 1., -1., ar, br, cr);

		/* calculate meridional arc length, area, and volume 
		 * using the extended trapezoidal rule */
		s[0] = 0.;
		A    = 0.;
		V    = 0.;
		for (i = 0; i < N; i++){
			// update current step
			 xi =  x[i];
			 ri =  r[i];
			axi = ax[i];
			bxi = bx[i];
			cxi = cx[i];
			ari = ar[i];
			bri = br[i];
			cri = cr[i];
			
			// initialize
			ssum = 0.;
			Asum = 0.;
			Vsum = 0.;
			dl   = l[i+1] - l[i];
			
			// evaluate integrand at l[i], l[i+1]
			for (j = i; j < i+2; j++){
				rj    =  r[j];
				dxdl  =  cxi;
				drdl  =  cri;
				dsdl  =  sqrt(dxdl*dxdl + drdl*drdl);
				dAdl  =  rj*dsdl;
				dVdl  = -rj*rj*dxdl;
				
				ssum += 0.5*dsdl;
				Asum += 0.5*dAdl;
				Vsum += 0.5*dVdl;
			}

			na = 2; // number of points added
			nt = 2; // total number of points

			for (n = 1; n < MAXIT; n++){
				// add 2^(n-1) additional points
				if (n == 1)
					na  = 1;
				else
					na *= 2;

				nt += na;

				// store integrands before adding new points
				ds = ssum*dl;
				dA = Asum*dl;
				dV = Vsum*dl;

				// refine grid spacing
				dl /= 2.;

				// interpolate integrand at additional points
				for (j = 0; j < nt; j++){
					if (j % 2 != 0){ /* avoid double counting
					                  * previous grid points */
						lj    = j*dl;

						rj    =  ((ari*lj + bri)*lj + cri)*lj + ri;
						dxdl  =  (3.*axi*lj + 2.*bxi)*lj + cxi;
						drdl  =  (3.*ari*lj + 2.*bri)*lj + cri;
						dsdl  =  sqrt(dxdl*dxdl + drdl*drdl);
						dAdl  =  rj*dsdl;
						dVdl  = -rj*rj*dxdl;
						
						ssum += dsdl;
						Asum += dAdl;
						Vsum += dVdl;
					}
				}

				// calculate change in integrand
				ds -= ssum*dl; ds = fabs(ds);
				dA -= Asum*dl; dA = fabs(dA);
				dV -= Vsum*dl; dV = fabs(dV);

				// break loop when integrals converge
				if (ds < stol && dA < Atol && dV < Vtol)
					break;
			}
		  
			// diagnose level of refinement [uncomment when required]
		 // printf("%d stages of refinement and %d total points\n", n, nt);
			
			// increment the integrals
			ds = ssum*dl;
			dA = Asum*dl;
			dV = Vsum*dl;

			s[i+1] = s[i] + ds;
			A     += dA;
			V     += dV;
		}

		A *= 2*M_PI;
		V *= M_PI;
		
		/* calculate principal curvatures, meridional tangent vector,
		 * and outward normal vector at the nodal points */
		ks[0]     =  0.;
		kp[0]     =  0.;
		
		tx[0]     =  0.;
		tr[0]     =  1.;
		
		nx[0]     =  1.;
		nr[0]     =  0.;

		for (i = 1; i < N; i++){
			ri      =  r[i];
			bxi     =  bx[i];
			bri     =  br[i];
			cxi     =  cx[i];
			cri     =  cr[i];
			dxdl    =  cxi;
			d2xdl2  =  2*bxi;
			drdl    =  cri;
			d2rdl2  =  2*bri;
			dsdl    =  sqrt(cxi*cxi + cri*cri);

			ks[i]   = -(dxdl*d2rdl2 - d2xdl2*drdl)/pow(dsdl, 3);
			kp[i]   =  dxdl/(ri*dsdl);

			tx[i]   =  dxdl/dsdl;
			tr[i]   =  drdl/dsdl;
			
			nx[i]   =  tr[i];
			nr[i]   = -tx[i];
		}
	
		fc        = (s[0]   - s[2])/(s[1] - s[2]);
		ks[0]     = ks[1]   + fc*(ks[1]   - ks[2]);
		kp[0]     = kp[1]   + fc*(kp[1]   - kp[2]);

		fc        = (s[N] - s[N-2])/(s[N-1] - s[N-2]);
		ks[N]     = ks[N-1] + fc*(ks[N-1] - ks[N-2]);
		kp[N]     = kp[N-1] + fc*(kp[N-1] - kp[N-2]);
		
		tx[N]     =  0.;
		tr[N]     = -1.;
		
		nx[N]     = -1.;
		nr[N]     =  0.;
	}
};

#endif
