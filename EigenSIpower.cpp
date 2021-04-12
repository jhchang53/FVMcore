/*
 * 	EigenSIpower.cpp
 * 	shift and invert solver
 */
#include <stdio.h>
#include "EigenSIpower.h"

EigenSIpower::EigenSIpower()
{
  ASmat = NULL;
  Bmat = NULL;
};

EigenSIpower::~EigenSIpower()
{
  MatDestroy(&ASmat);
  MatDestroy(&Bmat);
};	// EigenSIpower::~EigenSIpower

void EigenSIpower::setProblem(Mat Amat_i, Mat Bmat_i, double sigma_i)
{
  Bmat = Bmat_i;
  sigma = sigma_i;
  ierr = MatDuplicate(Amat_i,MAT_COPY_VALUES,&ASmat);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(ASmat,-sigma,Bmat,SUBSET_NONZERO_PATTERN);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  /*  create KSP  */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = KSPSetOperators(ksp,ASmat,ASmat); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = KSPSetTolerances(ksp,KSP_rtol,KSP_atol,PETSC_DEFAULT,KSP_itmax);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  PC pc;
  KSPGetPC(ksp,&pc);
  ierr = KSPSetInitialGuessNonzero(ksp,PETSC_TRUE);
  ierr = KSPSetUp(ksp); CHKERRABORT(PETSC_COMM_SELF,ierr);

  /*    prepare index for vecgetrows */
  PetscInt nrows,ncols;
  MatGetSize(ASmat,&nrows,&ncols);
  ndim = nrows;
  indunit = new int[ndim];
  for(int i=0; i < ndim; i++) indunit[i] = i;
  /*    prepare storage  */
  ierr = VecCreate(PETSC_COMM_WORLD,&phi);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(phi,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(phi); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&rhs);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(rhs,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(rhs); CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// EigenSIpower::setProblem

void EigenSIpower::guess(Vec phi_i)
{
  ierr = VecCopy(phi_i,phi);    CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // EigenSIpower::guess


double EigenSIpower::solve()
{
  double phinorm,phinormold;
  double rhsnorm;
  VecNorm(phi,NORM_2,&phinorm);
  if(prlev && (mpid==0)) printf("phinorm=%.3le\n",phinorm);
  if(phinorm < 1.0e-20) {
    printf("*** %s:%d phi is not proper. norm=%.2le\n",__FILE__,__LINE__,
        phinorm);
    exit(0);
  }
  phinormold = 99.99;
  int converged = 0;
  for(int iter=0;!converged && (iter < Eig_itmax); iter++) {
    VecScale(phi,1.0/phinorm);
    ierr = MatMult(Bmat,phi,rhs);  CHKERRABORT(PETSC_COMM_SELF,ierr);
    solveAxp(rhs,phi);
    VecNorm(phi,NORM_2,&phinorm);
    if(prlev && (mpid==0)) printf("= iter=%d phinorm=%.5le err=%.2le\n",
	iter,phinorm,
	phinorm-phinormold);
    if(fabs(phinorm-phinormold) < Eig_rtol) converged = 1;
    else phinormold = phinorm;
  }
  if(!converged) {
    printf("*** %s:%d eigenvalue Not converged \n",__FILE__,__LINE__);
    exit(0);
  }
  double eig = 1.0/(1.0/phinorm + sigma);
  if(prlev && (mpid==0))printf(" eig=%.5le\n",eig);
  return eig;
};	//  EigenSIpower::solve
