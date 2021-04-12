/*
 * 	NeutFVM.h
 * 	Neutron solver using FVM
 *    Usage:
 *	NeutFVM(); - constructor
 * 	~NeutFVM(); - destructor
 *	void setGeom(Geom *geom) - set geometry
 *	  geom : class describing 3D primatic geometry
 *  	void setXS(XS *xs, double Tmref, double Tfref);
 *  	 - set element wise cross section using database
 *  	  xs : cross section data base class
 *  	  Tmref : reference moderator temp.
 *  	  Tfref : reference fuel temp.
 *	void setTemp(double *Tm, double *Tf); - set element-wise temp
 *	  Tm : array of moderator temperature (size ntri*nz)
 *	  Tf : array of fuel temperaure (size ntri*nz)
 *	void setRodXS(int nrods, double sigrod[]); - set rod cross section
 *	  nrods : number of rod groups
 *	  sigrod: delta absorption cs. of control rods
 *	void setRodPos(int nrods, double rodpos[]); - set rod position
 *	  nrods : number of rod groups
 *	  rodpos: position of rod from top of problem (meters)
 *	void prepare(int ngrps);
 *	 - set number of neutron groups and prepare storage
 *	  ngrps: number of neutron groups ( = 1)
 *	void setupAB(); - setup system matrix
 *	void setEigShift(double shift); - set shift of eigenvalue solver
 *	  shift: eigenvalue shift (default 0.9)
 *	double solveEigen(); - solver eigenvalue problem and
 *	 returns eigenvalue (keff)
 *	Vec makePower(double totpow); - normalize flux to given power
 *	 returns Vector holding element-wise power
 *	  totpow : total power for problem domain (Watt)
 *	void printPow(); - print channel and axil averaged power
 *	 * it could be very large.
 *
 *	void start(double totpow, double keff, double beta, double gvel[]);
 *	 - prepare storage etc. to start transient calculation
 *	  totpow: initial total power of problem domain (Watt)
 *	  keff  : cricality normalization factor
 *	  beta  : delayed neutrin fraction
 *	  gvel  : neutron speed (m/s)
 *	void setLTE(double epsR, double epsA); - set required tolerances
 *	  epsR : relative tolerance (default: 1.0e-5)
 *	  epsA : absolute tolerance (default: 1.0e-5)
 *	double step(double dt); - perform one time step calculation
 *	  dt : time step (sec)
 *	void march();	- advance time step for next time
 *	double getTotPow(); - returns total power
 *	double getPeakPowDen(); - returns peak power density
 *	Vec getPow(); - returns element-wise power prod.
 *  	Vec getFiss(); - returns element-wise fission 
 *	void setDsrc(Vec dsrcvec); - set delayed neutron source
 *	  dsrcvec : delayed neutrons source (#/s)
 *	Vec getPowden(); - returns power density vector
 */
#include <petscksp.h>
#include "Geom.h"
#include "XS.h"
#include "EigenSIpower.h"
#include "Solver.h"

class NeutFVM {
public:
  NeutFVM();
  ~NeutFVM();
  void setPrlev(int prlev);
  void setGeom(Geom *geom);
  void setXS(XS *xs, double Tmref, double Tfref);

  void setTemp(double *Tm, double *Tf);
  void setRodXS(int nrods, double sigrod[]);
  void setRodPos(int nrods, double rodpos[]);

  void prepare(int ngrps);
  void setupAB();
  void setEigShift(double shift);
  double solveEigen();

  Vec makePower(double totpow);
  void printPow();

  void start(double totpow, double keff, double beta, double gvel[]);
  void setLTE(double epsR, double epsA);
  double step(double dt);
  void march();
  double getTotPow();
  double getPeakPowDen();
  Vec getPow();
  Vec getFiss();
  void setDsrc(Vec dsrcvec);
  Vec getPowden();
private:
  void setupABmat();
  void setupPmat(double fact);
  void setupSmat();
  void makeVolvec();
  void setupMmat();	/* for transient  */
  void setupTransAmat(double fact); /* for transient */
  void makeCinit();

  void step_gam(double dt);
  void step_new(double dt);
  double getLTE(double dt);
  int prlev,mpid;
  int do_setup,do_power;
  int ntri,nz,nvols;
  int npows;
  int ngrps;
  int nrods;
  double Tmref,Tfref;
  Geom *geom;
  XS *xs;
  double *Tm,*Tf;
  double *sigrod,*rodpos;

  PetscErrorCode ierr;
  int ndim;
  int jrange,rAbeg,rAend;
  Mat Amat,Bmat,Pmat;
  Vec phivec,powvec;
  Vec volvec,powden;
  double eig_shift;
  Eigen *eigen;
  /*	transient  */
  Solver *solver;
  double effk;
  double beta;
  double *gvel;
  double epsR,epsA;
  double TRgam,TRC2;    // TRBDF2 parameters

  int ksp_itmax;
  double ksp_rtol,ksp_atol;
  Mat Mmat,ATmat,TRmat;
  Vec srcvec,fissvec;
  Vec phigam,phinew,rhsvec,workvec;
  Mat Smat;
};
