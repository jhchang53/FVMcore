#ifndef _inc_TriCond
#define _inc_TriCond
/*
 * 	TriCond.h
 * 	spherical FVM problem
 * 	calculate temperature fluctuation due to kernel power.
 * 	by net zero power generattion and net zero average temperature
 *
 * 	Usage:
 * 	 TriCond();	- constructor
 * 	 ~TriCond();	- destructor
 * 	 void setGeom(); - set geometry (built in)
 * 	   modify this routine to change dimensions of kernel and coatings
 *       void setup(int npows);  - setup required storages
 *         npows : number of fuel segments producing power
 * 	 void setCond(int nmats, double cond[], double rhocp[]);
 * 	   nmats: no. of materials (should be 5)
 * 	   cond : conductivity array (kernel/buffer/PyC/SiC/Matrix) (W/m/K)
 * 	   rhocp : heat capacity array (J/m^3/K)
 * 	 void setCond(int nmats, Vec condvec[], Vec rhocpvec[]);
 *	   vector version
 *	   condvec, rhocpvec is nmats vector array of size npows
 *
 *	  void setPfactor(double pfact)
 *	    pfact: multiplier to input power density
 *	     ( used to convert core power density to compact power density)
 * 	 void steady(double pker); - solve steady state problem
 * 	   pker : kernel power (Watt)
 * 	 void steady(Vec pkervec);
 * 	   vector version
 * 	   pkervec: size npows
 *
 * 	 void start();  - prepare storage for transient run
 * 	 void setLTE(double epsR, double epsA);
 * 	   - set converg. criteria of transient run using TRBDF-2
 * 	   epsR : relative error (default 1.0e-5)
 * 	   epsA : absolute error (default 1.0e-5)
 * 	 double step(double dt, double pker); - time stepping
 * 	   returns LTE estimator (time step is when <1: narrow, > 1: wide)
 * 	   dt : time step size (sec)
 * 	   pker : kernel power at the end of the step (Watt)
 * 	 double step(double dt, Vec pkervec);
 * 	   vector version
 * 	   pkervec : size npows
 * 	 void march(); - prepare for next time step
 * 	   * variables are preserved until march is called.
 * 	 void dumpT(); - dump result (Tvec)
 * 	 double getTker();
 * 	   - returns average of kernel temperature rise (all case)
 * 	 Vec getDTker(); - get kernel temperature rise in vector (npows)
 * 	   * add this to the compact temperature from SegCond 
 * 	     to obtain the fuel temperature
 *
 */
#include <petscksp.h>
#include "Solver.h"

class TriCond {
public:
  TriCond();
  ~TriCond();
  void setPrlev(int prlev);
  void setGeom();
  void setup(int npows);
  void setCond(int nmats, double cond[], double rhocp[]);
  void setCond(int nmats, Vec condvec[], Vec rhocpvec[]);

  void setPfactor(double pfact);
  void steady(double pker);
  void steady(Vec pkervec);

  void start();
  void setLTE(double epsR, double epsA);
  double step(double dt, double pker);
  double step(double dt, Vec pkervec);
  void march();

  void dumpT();
  double getTker();
  double getTavg();
  Vec getDTker();
private:
  void prepare();
  void setupKPmat();
  void setupQvec();
  void setupAmat();
  void zeroTavgMat(Mat Amat);
  void zeroTavgRhs(Vec rhs);

  void setupMmat();
  void step_gam(double dt);
  void step_new(double dt);
  double getLTE(double dt);

  int prlev,mpid;
  double PI;
  double pfact;
  int do_setup;
  int npows;
  int nzone,kernel;
  int nmats;
  double vKernel,vParticle;
  int *mat;
  double *vol,*sl,*dL,*dR;
  double *volbig;
  double *cond,*rhocp;
  Solver *solver;
  PetscErrorCode ierr;
  int ndim;
  int jrange,rAbeg,rAend;
  int *cols;
  Mat Amat;
  Vec Qkvec,Qvec,Qfvec;	// power density multiplier to have net zero production
  Vec Tvec,srcvec;
  Mat KPmat;	// for multiple case
  Vec workvec,work2vec;
  Vec kpwork,dTvec;	// size npows vector
  /*  for transients */
  double epsR,epsA;
  double TRgam,TRC2;    // TRBDF2 parameters
  Mat Mmat,TRmat;
  Vec Tgam,Tnew,srcnew;
  Vec rhsvec;
};
#endif
