/*
 *   solver
 */
#include "Solver.h"

Solver::Solver()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  ksp = NULL;
  itmax_ksp = 1000;
  rtol_ksp = 1.0e-9;
  atol_ksp = 1.0e-9;
};

Solver::~Solver()
{
  KSPDestroy(&ksp);
  ksp = NULL;
};	// Solver::~Solver

void Solver::setTol(int itmax, double rtol, double atol)
{
  itmax_ksp = itmax;
  rtol_ksp = rtol;
  atol_ksp = atol;
};	// Solver::setTol

int Solver::solve(Mat Amat, Vec rhsvec, Vec xvec)
{
  ksp_setup(Amat);
  int iter = ksp_solve(rhsvec,xvec);
  ksp_destroy();
  return iter;
};	// Solver::solve

void Solver::ksp_setup(Mat Amat)
{
  /*  solve  X * phi = rhs  */
  /*   setup KSP solver */
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp); CHKERRABORT(PETSC_COMM_SELF,ierr);
   ierr = KSPSetType(ksp,KSPBCGS);
  ierr = KSPSetOperators(ksp,Amat,Amat); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = KSPSetTolerances(ksp,rtol_ksp,atol_ksp,PETSC_DEFAULT,itmax_ksp);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  PC pc;
  KSPGetPC(ksp,&pc);
  PCSetType(pc,PCSOR);
  ierr = KSPSetUp(ksp); CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // Solver::ksp_setup

int Solver::ksp_solve(Vec rhsvec, Vec xvec)
{
  /*  solve   A * x = rhs   */
  ierr = KSPSolve(ksp,rhsvec,xvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  KSPConvergedReason kspreason;
  ierr = KSPGetConvergedReason(ksp,&kspreason);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(kspreason < 0) {
    if(mpid==0) {
      printf("\n");
      printf("*** %s:%d notconverged. Failed reason=%d\n",
        __FILE__,__LINE__,kspreason);
      if(kspreason==KSP_CONVERGED_RTOL) printf("*** KSP_CONVERGED_RTOL\n");
      else if(kspreason == KSP_DIVERGED_ITS) printf("*** KSP_DIVERGED_ITS\n");
      else if(kspreason == KSP_DIVERGED_NANORINF)
        printf("*** KSP_DIVERGED_NANORINF\n");
    }
    PetscReal rtol,abstol,dtol;
    PetscInt itmax;
    KSPGetTolerances(ksp,&rtol,&abstol,&dtol,&itmax);
    if(mpid==0) {
      printf("  itmax=%d rtol=%.1le abstol=%.1le dtol=%.1le\n",
        itmax,rtol,abstol,dtol);
    }

    if(prlev > 1) {
      printf("== %s:%d rhsvec:\n",__FILE__,__LINE__);
      VecView(rhsvec,PETSC_VIEWER_STDOUT_WORLD);
    }
    return -1;
  }
  if(prlev > 1) {
    printf("== %s:%d phivec:\n",__FILE__,__LINE__);
    VecView(xvec,PETSC_VIEWER_STDOUT_WORLD);
    KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
  }
  int its;
  KSPGetIterationNumber(ksp,&its);
  if(prlev > 1) printf("== ksp_solve: its=%d\n",its);
  if(prlev > 5) exit(0);
  return its;
};      // Sol_kspver::ksp_solve

void Solver::ksp_destroy()
{
  /*  destroy  */
  ierr = KSPDestroy(&ksp);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ksp = NULL;
};      // Neut::ksp_destroy
