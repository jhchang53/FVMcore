//
//	Eigen.cpp
//	generic class
//
#include <stdio.h>
#include <stdlib.h>
#include "Eigen.h"
#include "petscmat.h"

Eigen::Eigen()
{
  prlev = 0;
  indunit =NULL;
  a = NULL;
  ksp = NULL;
  KSP_itmax = 500;
  KSP_rtol = 1.0e-8;
  KSP_atol = 1.0e-8;
  Eig_itmax = 500;
  Eig_rtol = 1.0e-8;
  MPI_Comm_rank(MPI_COMM_WORLD,&mpid);
  rhs = NULL;
  phi = NULL;
};	// Eigen::Eigen

Eigen::~Eigen()
{
  delete [] indunit;
  delete [] a;
  KSPDestroy(&ksp);
  VecDestroy(&rhs);
};	// Eigen::~Eigen

void Eigen::setPrlev(int prl)
{
  prlev = prl;
};	// Eigen::setPrlev

void Eigen::setKSPtol(int itmax, double rtol)
{
  KSP_itmax = itmax;
  KSP_rtol = rtol;
};	// Eigen::setKSPtol

void Eigen::setEigtol(int itmax, double rtol)
{
  Eig_itmax = itmax;
  Eig_rtol = rtol;
};      // Eigen::setEigtol

void Eigen::setProblem(Mat Amat_i, Mat Bmat_i, double sigma)
{
  printf("*** %s:%d setProblem not implemented yet.\n",
	__FILE__,__LINE__);
  exit(0);
};	// Eigen::setProblem

void Eigen::guess(Vec phi)
{
  printf("*** %s:%d guess() not implemented.\n",__FILE__,__LINE__);
  exit(0);
};

double Eigen::solve()
{
  printf("= %s:%d\n",__FILE__,__LINE__);
  exit(0);
};	// Eigen::solve

void Eigen::solveAxp(Vec srcvec, Vec rhsvec)
{
  /*  solution  */
  ierr = KSPSolve(ksp,srcvec,rhsvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 2) {
    double snorm,rnorm;
    ierr = VecNorm(srcvec,NORM_2,&snorm);  CHKERRABORT(PETSC_COMM_SELF,ierr);
    ierr = VecNorm(rhsvec,NORM_2,&rnorm);  CHKERRABORT(PETSC_COMM_SELF,ierr);
    if(mpid == 0) {
      printf("solveAxp src=%.2le rhs=%.2le\n",snorm,rnorm);
      fflush(stdin);
    }
  }

  KSPConvergedReason kspreason;
  ierr = KSPGetConvergedReason(ksp,&kspreason);
  if((kspreason < 0) && (mpid == 0)) {
    if(kspreason == KSP_DIVERGED_PCSETUP_FAILED) {
      printf("*** %s:%d preconditioner problem.\n",__FILE__,__LINE__);
      printf("   diaginal may be zero.\n");
      exit(0);
    }
    printf("*** %s:%d notconverged. Failed reason=%d\n",
        __FILE__,__LINE__,kspreason);
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    PetscReal rtol,abstol,dtol;
    PetscInt maxits;
    KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&maxits);
    printf("  rtol=%.1le abstol=%.1le dtol=%.1le maxits=%d\n",
        rtol,abstol,dtol,maxits);
    exit(0);
  }
  if(prlev && (mpid==0)) {
    if(prlev > 1) KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
    int its;
    KSPGetIterationNumber(ksp,&its);
    if(mpid==0) printf("its=%d\n",its);
  }
};      // Eigen::solveAxp

