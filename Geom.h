#ifndef _inc_Geom
#define _inc_Geom
/*
 * 	Geom.h - prisimatic geometry description base class for FEM/FVM
 * 	2D structure should be defined at derived class
 *     Usage:
 *	int extrude(int nz, double dz[], int ispowz[]);
 *	 - extrude 2D structure to 3D
 *	  returns 3D number of points for FEM
 *	  nz : number of axial nodes
 *	  dz : array of z-interval from top of core (meter)
 *	  ispowz : array of flag for power producing node
 *	int extrude(int nz, int zseg[], double dz[], int ispowz[]);
 *	 - expands input and call 3 argumnent function extrude
 *	  returns 3D number of points for FEM
 *	  nz : number of axial zones
 *	  zseg : number of finer axial zone division
 *	  dz : array of z- zone interval from top of core (meter)
 *	  ispowz : array of flag for power producing zone
 *	void printMat3D(); - print material assignments
 *	double getNodes(int t,int z, double x2d[], double y2d[], int node[]);
 *	 - get node coordinates and point numbers to be used for quadratic FEM
 *	  returns dz of the element
 *	  t : triangle numer (2D)
 *	  z : axial number
 *	  output:
 *	  x2d : array of x-coord. of triangle vertices (in ccw)
 *	  y2d : array of y-coord.
 *	  node : point numbers of the triangular prism (FEM15)
 *	int getNtri()	- returns number of triangle elements (in 2D)
 *	int getNz()	- returns number of axial division
 *	double *Dz()	- returns pointer to axial division
 *	double *Zpos()	- returns z-poisition from top
 *	int getNpts3d() - returns number of 3D points(FEM15)
 *	int getNchan()	- returns number of triangle with coolant holes
 *	int getNZpow()	- returns number of power producing nodes in z-dir
 *	double *Zpow()	- returns pointer to mid point position of element
 *	int getNpows()	- returns total number of power producing elements
 *	int getNpts2d() - returns number of triangles vertices in 2D (FEM15)
 *	double *Px()	- returns pointer to x-coord. triangle vert. (FEM15)
 *	double *Py()	- returns pointer to y-coord. triangle vert. (FEM15)
 *	int *Tri()	- returns pointer to triangle point desription (FEM15) 
 *	int getNbdry()	- returns number of boundary curves (FEM15)
 *	int *Btri()	- returns pointer to boundary triangles (FEM15)
 *	int *Bedge()	- returns pointer to boundary edge description (FEM15)
 *	int *Top()	- returns pointer to point no. at top of prism (FEM15)
 *	int *Ispp()	- retunns pointer to pow. prod. channel flag
 *	int *Tchan() 	- returns pointer to tri no. of channel
 *	int *Chtri() 	- returns pointer to channel no. of triangle
 *	int *Ispowz()	- returns pointer to axial pow.prod. flag
 *	== material and control rod ==
 * 	int *Mat2D()	- returns pointer to 2D description of material index
 *	int *Mat3D()	- returns pointer to 3D expanded material index
 *	 - material index is set by derived class
 *	int *Crod2D()	- returns pointer to control rod number array
 *	 - control rod number is set by derived class
 *	== for XS translation  ==
 * 	int getNzinp()	- returns nz given at 4 arg. extrude call
 *	int *Zseg()	- returns pointer to zseg given at 4 arg. extrude call
 *	== for FVM  ==
 *	double getFVM(int t, int z, int tb[], double area[], double dist[])
 *	 - get surface data needed for FVM
 *	  returns volume of the element (t,z)
 *  	  t : triangle number
 *	  z : axial number
 *	  output:
 *	  tb: neighbor element number (tz). size 5 array (top/bottom/3 sides)
 *	  area : surface areas. size 5 array
 *	  dist : distance to the surface from the center of element. size 5 array
 *	double getVol(int t, int z);
 *	 - returns volume of element
 *	  t : triangle number
 *	  z : axial number
 * 	double getDistL(int etb,  int tz);
 * 	 - returns distance from triangle tz toward face etb
 * 	  etb : surface number (0-4)
 * 	  tz : element number (tz=t*nz+z)
 * 	double getChanArea();
 * 	 - returns total triangle area occupied by channel
 *	double getPowVol();
 *	 - return total volume of power producing elements
 */

class Geom
{
public:
  Geom();
  ~Geom();
  void setPrlev(int prl);
  int extrude(int nz, double dz[], int ispowz[]);
  int extrude(int nz, int zseg[], double dz[], int ispowz[]);
  void printMat3D();

  double getNodes(int t,int z, double x2d[], double y2d[], int node[]);
  int getNtri() {  return ntri; };
  int getNz() { return nz; };
  double *Dz() { return dz; };
  double *Zpos() { return zpos; };
  int getNpts3d() { return npts3d; };
  int getNchan() { return nchan; };
  int getNZpow() {return nzpow; }
  double *Zpow() { return zpow; }
  int getNpows() { return npows; };
  int getNpts2d() { return npts2d; };
  double *Px() { return px; };
  double *Py() { return py; };
  int *Tri() { return tri; };
  int getNbdry() { return nbdry; };
  int *Btri() { return btri; };
  int *Bedge() { return bedge; };

  int *Top() { return top; };
  int *Ispp() { return ispp; };
  int *Tchan() { return tchan; };	// tri no. of channel
  int *Chtri() { return chtri; };	// channel no. of tri
  int *Ispowz() { return ispowz; };

  int *Mat2D() { return mat2d; };
  int *Mat3D() { return mat3d; };

  int *Crod2D() { return crod2d; };

  /*  for XS translation  */
  int getNzinp() { return nzinp; };
  int *Zseg() { return zseg; };

  /*  for FVM  */
  double getFVM(int t, int z, int tb[], double area[], double dist[]);
  double getVol(int t, int z);
  double getDistL(int etb,  int tz);

  double getChanArea();
  double getPowVol();

protected:
  void findNeighbor();
  void makeFVM();

  int npts2d,ntri;
  int nbdry;	// no. of boundary triangles
  int nchan;	// no. of channels
  double *px,*py;
  int *tri;
  int *btri,*bedge;
  int *tchan,*ispp;
  int *mat2d,*mat3d;
  int nz;
  int nzinp;	// used for input spec
  int *zseg;	// no. of axial segments
  int *crod2d;	// control rod
  /*  FVM  */
  int *nbtri;
  double *FVM_area,*FVM_edge,*FVM_dist;

private:
  void prepareTH();
  double Area(int t);

  int prlev;
  int npts3d;
  int *coner,*top;
  double *dz,*zpos;
  /*  flow/power production  */
  int nzpow,npows;
  int *ispowz,*chtri;
  double *zpow;
};
#endif
