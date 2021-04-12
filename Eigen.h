#ifndef _inc_Eigen
#define _inc_Eigen
/*
 *	Eigen.h
 *	generic class
*/
#include <petscksp.h>
class Eigen {
public:
  Eigen();
  ~Eigen();
  void setPrlev(int prlev);
  void setKSPtol(int itmax, double rtol);
  void setEigtol(int Eig_itmax, double Eig_rtol);
  virtual void setProblem(Mat Amat, Mat Bmat, double sigma);
  virtual void guess(Vec phi);
  virtual double solve();
  Vec getPhi() { return phi; };
protected:
  void solveAxp(Vec src, Vec rhs);
  int prlev,mpid;
  int ndim,nnz;
  int *irn,*jcn;
  double *a;
  PetscErrorCode ierr;
  int *indunit;
  Mat Amat,Bmat;	// pointer
  int *BCzero;
  KSP ksp;
  int KSP_itmax;
  double KSP_rtol,KSP_atol;
  int Eig_itmax;
  double Eig_rtol;
  Vec rhs;
  Vec phi;	// result
};
#endif
