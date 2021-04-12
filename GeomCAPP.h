#ifndef _inc_GeomCAPP
#define _inc_GeomCAPP
/*
 * 	GeomCAPP.h
 * 	desribe CAPP type geometry. inherited from class Geom
 *     Usage:
 *	GeomCAPP(); 	- constructor
 *	~GeomCAPP();	- destructor
 * 	void setPlane(double side, int nrow, int mat[]); - basic layout
 * 	  side : side length of a triable (we assume eqi-lat. triangle)
 * 	  nrow : no. of rows (size of mat array is 1+nrow*(nrow-1)/2)
 * 	  mat : array of material index to describe the hexagonal geometry
 * 	    -1 for empty place.
 * 	void setChannel(int ischan[]);	- coolant channel flag
 * 	  ischan : flag  for channel (1=coolant, 0=block)
 * 	void setCrod(int crod[]); - set control rod number
 * 	  crod : array of control rods numers (start from 1)
 *	example)
 *  	  int mat[]  = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
 *	  int chan[] = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
 *	  int crod[] = {1, 0, 2,0, 0,0,0, 0,0,0,0, -1,0,0,0,0};
 *	  geom->setPlane(side, 6, mat);
 *	  geom->setChannel(chan);
 *	  geom->setCrod(crod);
 *
 *	void setHexMat(int h, int matz[]); - set axial material index
 *	 to override 2D material index given by setPlane.
 *	 should be called after extrude
 *	  h : the hexagon number
 *	  matz : axial matrial index
 *  	== for XS mapping  ==
 *	int getNrow()	- returns rows given at setPlane
 *	int *Map() - returns pointer to hexagon map
 *      
 */
#include "Geom.h"

class GeomCAPP : public Geom {
public:
  GeomCAPP();
  ~GeomCAPP();
  void setPrlev(int prl);
  void setPlane(double side, int nrow, int mat[]);
  void setChannel(int ischan[]);
  void setCrod(int crod[]);
  void setHexMat(int h, int matz[]);
  /*  setHexMat: set material index in z-direction of hexagon h  */
  /*  h : hexagon number  */
  /*  matz[]: material numbers  */
  /*  for XS mapping  */
  int getNrow() { return nrow; };
  int *Map() { return map; };
private:
  void genMap();
  void genHbdry();
  void makePlane(double side);
  void makeType1(int ptvar[], int tup,
	double cx, double cy, double s, double h);
  void makeType2(int ptvar[], int tleft, int tdn,
        double cx, double cy, double s, double h);
  void makeType3(int ptvar[], int tup, int tleft, int tdown,
        double cx, double cy, double s, double h);
  void makeType4(int ptvar[], int tup, int tleft,
        double cx, double cy, double s, double h);
  void setBdry();
  void setMat2D(int mat[]);
#ifdef XX
  void findNeighbor();
  void makeFVM();
#endif
  int prlev;
  int nrow,nentry,nhex,nhtri;
  double rt3;
  int *hex;	// hex number
  int *map,*hbdrytype;	// rectangular hex map
  int *htri;	// triangle number of hexagon
  int *hexseg;
  double *ctx,*cty;

};
#endif
