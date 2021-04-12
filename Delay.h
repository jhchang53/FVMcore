/*
 * 	Delay.h
 * 	delayed neutron kinetics
 *     Usage:
 *	Delay(); - constructor
 *      ~Delay(); - destructor
 *      double setParam(int nfiss, int nprec, double dbeta[], double dlamda[]);
 *       - set delayed neutron paramers
 *        nfiss : number of fissioning element ( = npows)
 *        nprec : number of precursor groups
 *        dbeta : array of precursor prod. fraction per fission 
 *          * Sum_nprec dbeta[] = beta (used at Neut class)
 *        dlamda : decay constant of precursor group (sec^-1)
 *      void init(Vec fissvec);	- initialize delayed neutron
 *        fissvec : initial fission source vector (size nfiss)
 *      Vec getNeutron(); - returns neutrons to be used at Neut class
 *      double step(double dt, Vec fissvec);
 *       - compute delayed neutron after time dt
 *        returns relative error estimator for time step control
 *        dt : time step size (sec)
 *        fissvec : fission source vector at time after dt
 *	void march(); - move variable to begining
 *	 to start next time step calculation
 *
 */
#include <petscksp.h>

class Delay {
public:
  Delay();
  ~Delay();
  double setParam(int nfiss, int nprec, double dbeta[], double dlamda[]);
  void init(Vec fissvec);
  Vec getNeutron();
  double step(double dt, Vec fissvec);
  void march();
private:
  void prepare();

  int prlev,mpid;
  int nfiss,nprec;
  double *dbeta,*dlamda;
  PetscErrorCode ierr;
  int ndim;
  int jrange,rCbeg,rCend;
  int krange,rPbeg,rPend;
  int *idxP;
  double *valP;
  Vec fissvec,fissnew;
  Vec cdenvec,cdennew;
  Vec dneutvec;
};
