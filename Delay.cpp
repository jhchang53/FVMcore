/*
 * 	Delay.cpp
 * 	delayed neutron kinetics
 */
#include <stdio.h>
#include <assert.h>
#include "Delay.h"

Delay::Delay()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  dbeta = NULL;
  dlamda = NULL;
  cdenvec = NULL;
  idxP = NULL;
  valP = NULL;
  fissvec = NULL;
  fissnew = NULL;
  cdennew = NULL;
};

Delay::~Delay()
{
  delete [] dbeta;
  delete [] dlamda;
  VecDestroy(&cdenvec);
  delete [] idxP;
  delete [] valP;
  VecDestroy(&fissvec);
  VecDestroy(&fissnew);
  VecDestroy(&cdennew);
};

double Delay::setParam(int nfiss_i, int nprec_i, double dbeta_i[], double dlamda_i[])
{
  nfiss = nfiss_i;
  nprec = nprec_i;
  dbeta = new double[nprec];
  dlamda = new double[nprec];
  double beta = 0;
  for(int d=0; d < nprec; d++) {
    dbeta[d] = dbeta_i[d];
    beta += dbeta[d];
    dlamda[d] = dlamda_i[d];
  }
  prepare();
  if(mpid==0 && prlev) printf(" nfiss=%d nprec=%d\n",nfiss,nprec);
  return beta;
};	// Delay::setParam

void Delay::prepare()
{
  /*  prepare working areas */
  ndim = nfiss*nprec;
  ierr = VecCreate(PETSC_COMM_WORLD,&cdenvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(cdenvec,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(cdenvec);     CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecGetOwnershipRange(cdenvec,&rCbeg,&rCend);
  jrange = rCend-rCbeg;

  ierr = VecDuplicate(cdenvec,&cdennew);

  /*  prepare storage for delayed neutrons  */
  ierr = VecCreate(PETSC_COMM_WORLD,&fissvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(fissvec,PETSC_DECIDE,nfiss);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(fissvec);     CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecGetOwnershipRange(fissvec,&rPbeg,&rPend);
  krange = rPend-rPbeg;
  /*  index array  */
  idxP = new int[nfiss];
  valP = new double[nfiss];
  /*  save fissvec etc  */
  ierr = VecDuplicate(fissvec,&dneutvec);
  ierr = VecDuplicate(fissvec,&fissnew);
};

void Delay::init(Vec fissvec_i)
{
  /*  check size  */
  int szPow;
  ierr = VecGetSize(fissvec_i,&szPow);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  assert(szPow == nfiss);
  ierr = VecCopy(fissvec_i,fissvec);
  double *Fiss;
  ierr = VecSet(cdenvec,0.0);
  ierr = VecGetArray(fissvec,&Fiss);
  for(int k=0; k < krange; k++) {
    int n = rPbeg+k;
    int indr[8];	// max. delayed group = 8
    double val[8];
    for(int d=0; d < nprec; d++) {
      indr[d] = n*nprec+d;
 
      val[d] = dbeta[d]*Fiss[k]/dlamda[d];
    }
    ierr = VecSetValues(cdenvec,nprec,indr,val,INSERT_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = VecAssemblyBegin(cdenvec);

  ierr = VecRestoreArray(fissvec,&Fiss);

  ierr = VecAssemblyEnd(cdenvec);
  if(prlev > 1) {
    if(mpid==0) printf("=== %s:%d cdenvec\n",__FILE__,__LINE__);
    VecView(cdenvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};	// Delay::init

Vec Delay::getNeutron()
{
  ierr = VecSet(dneutvec,0.0);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int n=0; n < nfiss; n++) {
    idxP[n] = n;
    valP[n] = 0.0;
  }
  double *C;
  ierr = VecGetArray(cdenvec,&C);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int j=0; j < jrange; j++) {
    int nj = rCbeg +j;
    int n = nj/nprec;
    int d = nj%nprec;
    double lamdaC = dlamda[d]*C[j];
    valP[n] += lamdaC;
  }
  ierr = VecRestoreArray(cdenvec,&C);	CHKERRABORT(PETSC_COMM_SELF,ierr);

  int ibeg = -1;
  int iend = 0;
  for(int n=0; n < nfiss; n++) {
    double fval = fabs(valP[n]);
    if(fval > 0.0) iend = n;
    if(ibeg < 0) if(fval > 0) ibeg = n;
  }
  int idim = iend-ibeg+1;
  ierr = VecSetValues(dneutvec,idim,idxP+ibeg,valP+ibeg,ADD_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyBegin(dneutvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(dneutvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    if(mpid==0) printf("=== %s:%d dneutvec\n",__FILE__,__LINE__);
    ierr = VecView(dneutvec,PETSC_VIEWER_STDOUT_WORLD);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  return dneutvec;
};	//  Delay::getNeutron

double Delay::step(double dt, Vec fissvec_i)
{
  VecCopy(fissvec_i,fissnew);
  /*  prepare exponential multipliers  */
  double xmax = 0.0;
  double fc[8],fo[8],fn[8];	// max. 8 groups
  for(int d=0; d < nprec; d++) {
    double x = dlamda[d]*dt;
    if(xmax < x) xmax = x;
    double blam = dbeta[d]/dlamda[d];
    if(x < 0.01) {
      fc[d] = 1+x*(-1+x*(0.5-x/6));
      fo[d] = x*(0.5+x*(-1.0/3.0+x/6))*blam;
      fn[d] = x*(0.5+x*(-1.0/6.0+x/24))*blam;
    }
    else {
      double ex = exp(-x);
      fc[d] = ex;
      fo[d] = (1-ex-(x-1+ex)/x)*blam;
      fn[d] = (x-1+ex)/x*blam;
    }
  }
  if((mpid==0) && prlev) {
    for(int d=0; d < nprec; d++) {
      printf(" %d: fc=%.4le fo=%.4le fn=%.4le\n",d,fc[d],fo[d],fn[d]);
    }
  }
  /*  c_new = fc[d]*c_old + fo[d]*fiss_old + fn[d]*fiss_new  */
  /* process c_new  = fc[d]*c_o  */
  ierr = VecSet(cdennew,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  double *C;
  ierr = VecGetArray(cdenvec,&C);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int j=0; j < jrange; j++) {
    int nj = rCbeg+j;
    /* int n = nj/nprec;  */
    int d = nj%nprec;
    double cden = C[j]*fc[d];
    ierr = VecSetValues(cdennew,1,&nj,&cden,ADD_VALUES);
  }
  ierr = VecRestoreArray(cdenvec,&C);

  /*  process  */
  double *Po,*Pn;
  ierr = VecGetArray(fissvec,&Po);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(fissnew,&Pn);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int k=0; k < krange; k++) {
    int n = rPbeg+k;
    int indr[8];
    double val[8];
    for(int d=0; d < nprec; d++) {
      indr[d] = n*nprec+d;
      val[d] = fo[d]*Po[k]+fn[d]*Pn[k];
    }
    ierr = VecSetValues(cdennew,nprec,indr,val,ADD_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = VecAssemblyBegin(cdennew);

  ierr = VecRestoreArray(fissvec,&Po);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(fissnew,&Pn);   CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecAssemblyEnd(cdennew);
  if(prlev) {
    printf("=== %s:%d cdennew\n",__FILE__,__LINE__);
    VecView(cdennew,PETSC_VIEWER_STDOUT_WORLD);
  }
  /*  we return  nonlinear estimate  */
  return xmax*xmax/2;
};	// Delay::step

void Delay::march()
{
  VecCopy(fissnew,fissvec);
  VecCopy(cdennew,cdenvec);
};	// Delay::march
