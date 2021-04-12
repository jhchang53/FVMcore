#ifndef _inc_Solver
#define _inc_Solver
#include <petscksp.h>

class Solver {
public:
  Solver();
  ~Solver();
  void setTol(int itmax, double rtol, double atol);
  int solve(Mat Amat, Vec rhsvec, Vec xvec);
  void ksp_setup(Mat Amat);
  int ksp_solve(Vec rhsvec, Vec xvec);
  void ksp_destroy();
private:
  int mpid,prlev;
  PetscErrorCode ierr;
  int itmax_ksp;
  double rtol_ksp,atol_ksp;
  KSP ksp;
};
#endif
  
  
  
