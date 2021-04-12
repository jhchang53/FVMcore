#ifndef _inc_SegCond
#define _inc_SegCond
/*
 * 	SegCond.h
 * 	1D heat transfer to calculate
 * 	convective coolant / graphite / gap / compact
 * 	Finite Volume Method parameters is given
 *
 *	Usage:
 *	 SegCond();   - constructor
 *	 ~SegCond();  - destructor
 *	 void setGeom(); - set geometry (built in)
 *	   modify this routine to change dimensions
 *	 void setGeom1D(int nzone, int jheat, double size);
 *	   sample 1D geometry
 *	   nzone: number of uniform mesh division
 *	   jheat: starting node number of heat generation zone
 *	   size : length of problem (meter)
 *	 void setup(int npows); - setup required storages
 *	   npows :  number of fuel segments producing power
 *	 void setHcoef(double hcoef);
 *	   hcoef: heat transfer coefficient (W/m^2/K)
 *	 void setHcoef(Vec hcoefvec);
 *	   vector version
 *	   hcoefvec: size npows
 *	 void setCond(int nmats, double cond[], double rhocp[]);
 *	   nmats: no. of materials
 *	   cond :  conductivity array (graphite/block/compact) (W/m/K)
 *	   rhocp: heat cpacity array (J/m^3/K)
 *	 void setCond(int nmats, Vec condvec[], Vec rhocpvec[]);
 *	   vector version
 *	   condvec, rhocpvec is nmats vector array of size npows
 *
 *	 void setPfactor(double pfact)
 *	   pfact: multiplier to input power density
 *	     ( used to convert core power density to compact power density)
 *	 void steady(double Tbulk, double qden); - solve steady state problem
 *	   Tbulk: coolant bulk temperature (K)
 *	   qden : compact power density (W/m^3)
 *	 void steady(Vec Tbulkvec, Vec qdenvec);
 *	   vector version
 *	   qdenvec : size npows
 *
 *	 void start(); - prepare storage for transient run
 *	 void setLTE(double epsR, double epsA);
 *	   - set converg. criteria of transient run using TRBDF-2
 *	   epsR : relative error (default 1.0e-5)
 *	   epsA : absolute error (default 1.0e-5)
 *
 *	 double step(double dt, double Tbulk, double qden); - time stepping
 *	    returns LTE estimator (time step is when <1: narrow, > 1: wide)
 *	    dt   : time step size (sec)
 *	    Tbulk: coolant bulk temperature (K) at the end of the time step
 *	    qden : compact power density (W/m^3) at the end of the time step
 *	 double step(double dt, Vec Tbulkvec, Vec qdenvec);
 *	   vector version
 *	   Tbulk,qden : size npows
 *	 void march(); - prepare for next time step
 *	   * variables are preserved until march is called.
 *
 *	 double getTpeak(); - returns highest temperature of the problem
 *	 double getTpin();  - returns average of compact region temperarture
 *	 Vec getTpinVec();  - compact region avg. temperature in vector
 *	   * add this with the DTker to obtain the fuel temperature
 *	 Vec getTgraVec();  - graphite region avg. temp.
 *
 */
#include <petscksp.h>
#include "Solver.h"

class SegCond {
public:
  SegCond();
  ~SegCond();
  void setPrlev(int prlev);
  void setGeom();
  void setGeom1D(int nzone, int jheat, double size);
  void setup(int npows);
  void setHcoef(double hcoef);	// same value for all case
  void setHcoef(Vec hcoefvec);	// size npows vector
  void setCond(int nmats, double cond[], double rhocp[]);
  void setCond(int nmats, Vec condvec[], Vec rhocpvec[]);

  void setPfactor(double pfact);
  void steady(double Tbulk, double qden);
  void steady(Vec Tbulkvec, Vec qdenvec);

  void start();
  void setLTE(double epsR, double epsA);
  double step(double dt, double Tbulk, double qden);
  double step(double dt, Vec Tbulkvec, Vec qdenvec);
  void march();

  void dumpT();
  double getTpeakNew();	// peak T at the end of a time step
  double getTpeak();
  double getTpin();
  Vec getTpinVec();
  Vec getTgraVec();
private:
  void printGeom();
  void prepare();
  void setupKPmat();
  void setupWGmat();
  void setupAmat();
  void setupQvec(int jheat);
  void setupMmat();
  void step_gam(double dt);
  void step_new(double dt);
  double getLTE(double dt);

  int prlev,mpid;
  int do_setup;
  int npows;	// for multiple vector
  int nzone,jheat,ndim;
  double PI,rt3;
  double pfact;	// mulitplier to input power density
  int *mat;
  double *vol,*sl,*dL,*dR;
  int nmats;
  double *cond,*rhocp;	// size npows x nmats
  double *hcoef;	// size npows
  double Qvol;

  Solver *solver;
  PetscErrorCode ierr;
  int jrange,rAbeg,rAend;
  Mat Amat;
  Vec BCvec,Qvec;
  Vec srcvec,Tvec;
  Vec workvec,work2vec;
  Mat KPmat;	// for multiple vector  for fuel pin
  Mat WGmat;	// for graphite region weights
  Vec Tpinvec,Tgravec;
  /*  for transients */
  double epsR,epsA;
  double TRgam,TRC2;    // TRBDF2 parameters
  Mat Mmat,TRmat;
  Vec Tgam,Tnew,srcnew;
  Vec rhsvec;
};
#endif
