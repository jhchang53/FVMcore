/*
 * 	Cond3D.h
 * 	3D triangular prismic geometry thermal conduction/convection
 * 	using porous media model
 *
 * 	  Cond3D();	- constructor
 * 	  ~Cond3D();	- destructor
 * 	  void setGeom(Geom *geom); - specify problem geometry
 * 	    geom : class descrbing 3D prismatic geometry
 * 	  void setFluence(double *flu); - set fluence for Thermal property
 * 	    flu : array of fluence (10^25 n/m^2)
 * 	  void setTHpro(THpro *thpro);	- specify thermal property
 *	    thpro : class for thermal conductivity and heat capacity
 *	      * if thpro is NULL, builtin property values will be used
 * 	  void setMdot(double Mdot); - specifies flows rate
 * 	    Mdot : flow rate (kg/s) of problem domain (nchan triangle plane)
 * 	  void setTin(double Tin);
 * 	    Tin : inlet temperature (used for side BC too) (K)
 * 	  void setup();	- prepare storages and setup matrices
 * 	  void steady(Vec powden);
 * 	    powden : power density (W/m^3) for each power producing volume
 * 	    	* not power density
 * 	  void start();
 * 	    prepare for transient calculation
 * 	  void setLTE(double epsR, double epsA);
 * 	    tolerance TRBDF-2 time stepping
 * 	    epsR : relative tolerance (default 1.0e-5)
 * 	    epsA : absolute tolerance (default 1.0e-5)
 * 	  double step(double dt, Vec powden);
 * 	    returns LTE estimator (time step is when <1: narrow, > 1: wide)
 * 	  void march();  - prepare for next time step
 * 	   * variables are preserved until march is called.
 * 	  double getTpeak();
 * 	    returns maximum temperature
 * 	  double getTPavg();
 * 	    returns avg temperature of power producing volumes
 * 	  double getTexit();
 * 	    returns average of channel exit temperaure
 * 	  void printAx(int t);
 * 	    print axial Temp. of trianle t.
 * 	  double updateTemp()
 * 	    update TK and TB
 * 	  Vec getPhcoef();
 * 	    returns hcoef vector for power generation elements
 * 	  Vec getPTbulk();
 * 	    returns Tbulk vector for power generation elements
 * 	  double* getTK();
 * 	    returns TK in array to be used at Neutronics
 * 	  double* expandPTvec(Vec xvec);
 * 	    returns expanded Pvec into core wide array for Neutronics
 * 	    Pvec : vector of size npows
 * 	  Vec contractPdouble(double x[]);
 * 	    returns vectors in power region using core wide array for SegCond
 * 	    x : array of size ntri*nz
 */
#include <petscksp.h>
#include "Geom.h"
#include "THpro.h"
#include "Solver.h"

class Cond3D {
public:
  Cond3D();
  ~Cond3D();
  void setPrlev(int prlev);
  void setGeom(Geom *geom);
  void setFluence(double *flu);
  void setTHpro(THpro *thpro);
  void setMdot(double Mdot);
  void setTin(double Tin);
  void setup();
  void steady(Vec powden);
  void start();
  void setLTE(double epsR, double epsA);
  double step(double dt, Vec powden);
  void march();
  double getTpeak();
  double getTPavg();
  double getTexit();
  void printAx(int t);
  double updateTemp();
  double *getTK();	// get TK for core
  Vec getPhcoef();	// heat transfer coef. for SegCond
  Vec getPTbulk();	// coolant Bulk Temp. for SegCond
  double* expandPvec(Vec xvec);
  Vec contractPdouble(double *x);
private:
  void prepare();
  void makeChannelMap();
  void makeP2Smat();
  void makeTBmat();
  void makeToutmat();
  void makePowvol();
  void setupAmat();
  void setupMmat();
  void step_gam(double dt);
  void step_new(double dt);
  double getLTE(double dt);

  double updateTK();
  double updateTB();

  int prlev,mpid;
  int do_setup,do_updTK,do_updTB;
  int ntri,nz,nvols,nchan,npows;
  Geom *geom;
  THpro *thpro;
  double *TK,*flu;	// for conductivity calculation
  double *TB;		// for heat transfer coefficient calculation
  double *PTwork;	// core wide array for expandPTvec
  int *kndtri,*kndchn;	// kndex of mixed matrix
  int ndim;
  int jrange,rAbeg,rAend;
  double mdotcp,Tin;
  PetscErrorCode ierr;
  Mat Amat;
  Vec BCvec;
  Mat P2Smat;  // convert pow to src
  Vec TPvec;	// temp. of power producing elements (size npows)
  Vec powvec;
  Vec powvol;	// volume of power producing element
  Vec hcoefvec;	// hcoef for power producing element
  Vec Tbulkvec;	// Tbulk for power producing element
  Vec pxvec;	// work vector for contractPdouble
  Vec srcvec,workvec,Tvec;
  Solver *solver;
  Mat Toutmat;	// collect channel exit temperatures
  Vec Toutvec;	// temp. of channel exit (size nchan)
  Mat TKmat;	// collect volume temperatures
  Vec TKvec;	// size ntri*nz
  Mat TBmat;	// collect coolant temperatures
  Vec TBvec;	// size nchan*nz;
  /*  for transients */
  double epsR,epsA;
  double TRgam,TRC2;	// TRBDF2 parameters
  Mat Mmat,TRmat;
  Vec Tgam,Tnew,srcnew;
  Vec rhsvec;
};
