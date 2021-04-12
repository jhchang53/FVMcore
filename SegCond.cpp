/*
 * 	SegCond.cpp
 * 	one dimension conduction problem
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "SegCond.h"


SegCond::SegCond()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  PI = 2*acos(0.0);
  rt3 = sqrt(3.0);
  pfact = 1.0;
  jheat=0;
  mat = NULL;
  vol = NULL;
  sl = NULL;
  dL = NULL;
  dR = NULL;
  cond = NULL;
  rhocp = NULL;
  solver = new Solver();
  do_setup = 1;
  hcoef = NULL;	// it is array
  Amat = NULL;
  BCvec = NULL;
  Qvec = NULL;	// heat volume vector
  Tvec = NULL;
  srcvec = NULL;
  workvec = NULL;
  KPmat = NULL;	// for multiple case
  Tpinvec = NULL;
  WGmat = NULL;
  Tgravec = NULL;
  work2vec = NULL;
  /* transient  */
  TRgam = 2.0-sqrt(2.0);
  TRC2 =  2*fabs(((-3*TRgam+4)*TRgam-2)/(12*(2-TRgam)));
  epsR = 1.0e-5;
  epsA = 1.0e-5;
  Mmat = NULL;
  TRmat = NULL;
  rhsvec = NULL;
  Tgam = NULL;
  Tnew = NULL;
  srcnew = NULL;
};	// SegCond::SegCond

SegCond::~SegCond()
{
  delete [] mat;
  delete [] vol;
  delete [] sl;
  delete [] dL;
  delete [] dR;
  delete [] cond;
  delete [] rhocp;
  delete solver;
  delete [] hcoef;
  MatDestroy(&Amat);
  VecDestroy(&BCvec);
  VecDestroy(&Qvec);
  VecDestroy(&Tvec);
  VecDestroy(&srcvec);
  MatDestroy(&KPmat);
  VecDestroy(&Tpinvec);
  MatDestroy(&WGmat);
  VecDestroy(&Tgravec);

  VecDestroy(&workvec);
  VecDestroy(&work2vec);

  MatDestroy(&Mmat);
  MatDestroy(&TRmat);
  VecDestroy(&rhsvec);
  VecDestroy(&Tgam);
  VecDestroy(&Tnew);
  VecDestroy(&srcnew);

};	// SegCond::~SegCond

void SegCond::setPrlev(int prl)
{
  prlev = prl;
};	// SegCond::setPrlev

void SegCond::setGeom()
{
  /*  arrange from coolant hole to compact center  */
  double pitch = 0.035556;
  double rhole = 0.007;
  int nring = 4;
  double rpin = 0.0115;
  double gap = 0.0005;
  /*  this is generated using Gmsh2d/Cond2D  */
  /*  relative to graphite volume  */
  double gravol[] = {0.063,0.122,0.237,0.579};
  /*  distance from center to left/right side relative to pitch */
  double gradL[] = {0.041,0.055,0.078,0.135};
  double gradR[] = {0.042,0.054,0.073,0.132};
  /*  interface length  */
  double grasl[] = {0.079,0.167,0.279,0.097};
  int ngra = sizeof(gravol)/sizeof(double);

  nzone = ngra+1+nring;
  mat = new int[nzone];
  vol = new double[nzone];
  dL = new double[nzone];
  dR = new double[nzone];
  sl = new double[nzone];	// left side interface length
  /*  allocate graphite region */
  double volgra = rt3/4*pitch*pitch - PI*rhole*rhole/12
	- PI*(rpin+gap)*(rpin+gap)/6;
  sl[0] = PI*rhole/6;
  for(int j=0; j < ngra; j++) {
    mat[j] = 0;	// material is graphite
    vol[j] = volgra*gravol[j];
    dL[j] = gradL[j]*pitch;
    dR[j] = gradR[j]*pitch;
    sl[j+1] = grasl[j]*pitch;
  }
  /*  helium gap  */
  mat[ngra] = 1;	// helium
  vol[ngra] = PI*((rpin+gap)*(rpin+gap)-rpin*rpin)/6;
  dL[ngra] = gap/2;
  dR[ngra] = gap/2;
  sl[ngra] = PI*(rpin+gap)/3;
  /*  fuel region - divide equal volume */
  double ringvol = rpin*rpin/nring;
  double r0 = 0;
  for(int k=0; k < nring; k++) {
    int j = nzone-k-1;
    mat[j] = 2;	// compact
    vol[j] = PI/6*ringvol;
    double r = sqrt((k+1)*ringvol);
    sl[j] = PI*r/6;
    double rc = sqrt((k+0.5)*ringvol);
    dL[j] = r-rc;
    dR[j] = rc-r0;
    r0 = r;
  }
  /*  set conductivity  and specific heat */
  nmats = 3;
  cond = new double[nmats];
  rhocp = new double[nmats];
  cond[0] = 5.0;  rhocp[0] = 1.0e+6;      // graphite
  cond[1] = 0.25; rhocp[1] = 100.0;       // helium gap
  cond[2] = 5.0;  rhocp[2] = 1.0+6;      // compact
  /*  set heat generation region */
  jheat = ngra+1;
};	// SegCond::setGeom

void SegCond::setGeom1D(int nzone_i, int jheat_i, double size)
{
  /* simple 1D geometry  */
  nzone = nzone_i;
  jheat = jheat_i;
  double dx = size/nzone;
  double dy = 1.0;
  mat = new int[nzone];
  vol = new double[nzone];
  sl  = new double[nzone];
  dL  = new double[nzone];
  dR  = new double[nzone];
  for(int j=0; j < nzone; j++) {
    mat[j] = 0;
    vol[j] = dx*dy;
    sl[j] = dy;
    dL[j] = 0.5*dx;
    dR[j] = 0.5*dx;
  }
  /*  set conductivity  */
  cond = new double[1];
  rhocp = new double[1];
  cond[0] = 5.0;	rhocp[0] = 3000.0;
};

void SegCond::printGeom()
{
  printf("mat:");
  for(int j=0; j < nzone; j++) printf(" %8d",mat[j]);
  printf("\n");
  printf("vol:");
  for(int j=0; j < nzone; j++) printf(" %.1le",vol[j]);
  printf("\n");
  printf("sl :");
  for(int j=0; j < nzone; j++) printf(" %.1le",sl[j]);
  printf("\n");
  printf("dL :");
  for(int j=0; j < nzone; j++) printf(" %.1le",dL[j]);
  printf("\n");
  printf("dR :");
  for(int j=0; j < nzone; j++) printf(" %.1le",dR[j]);
  printf("\n");
};	// SegCond::printGeom

void SegCond::prepare()
{
  /*  prepare for npows case with same geometry  */
  ndim = nzone*npows;
  int maxdiag = ndim;
  int maxoff = ndim;
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &Amat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Amat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Amat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = MatGetOwnershipRange(Amat,&rAbeg,&rAend);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
  jrange = rAend-rAbeg;

  ierr = VecCreate(PETSC_COMM_WORLD,&Tvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(Tvec,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(Tvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(Tvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(Tvec,&BCvec);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecDuplicate(Tvec,&Qvec);
  VecDuplicate(Tvec,&srcvec);
  VecDuplicate(Tvec,&workvec);
  VecDuplicate(Tvec,&work2vec);

  /*  for multiple case  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &KPmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(KPmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(KPmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &WGmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(WGmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(WGmat); CHKERRABORT(PETSC_COMM_SELF,ierr);


  ierr = VecCreate(PETSC_COMM_WORLD,&Tpinvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(Tpinvec,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(Tpinvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(Tpinvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(Tpinvec,&Tgravec);

  setupKPmat();
  setupWGmat();

  /*  hcoef array  */
  hcoef = new double[npows];
  for(int kp=0; kp < npows; kp++) hcoef[kp] = 1000.0;

  if(npows > 1) {
    /* expand cond and rhocp array  */
    double *cond_old  = cond;
    double *rhocp_old = rhocp;
    cond = new double[npows*nmats];
    rhocp = new double [npows*nmats];
    for(int kp=0; kp < npows; kp++) {
      for(int m=0; m < nmats; m++) {
        cond [kp*nmats+m] = cond_old[m];
        rhocp[kp*nmats+m] = rhocp_old[m];
      }
    }	// for kp
    delete [] cond_old;
    delete [] rhocp_old;
  }
};	// SegCond::prepare

void SegCond::setupKPmat()
{
  /*  matrix to translate multiple case input */
  ierr = MatZeroEntries(KPmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  int *indr = new int[nzone];
  double *val = new double[nzone];
  for(int n=0; n < nzone; n++) val[n] = 1.0;
  int indc[1];
  for(int kp=0; kp < npows; kp++) {
    /*  when all rows of a case is out side we skip  */
    int indfrom = kp*nzone;       int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    for(int n=0; n < nzone; n++) {
      indr[n] = kp*nzone+n;
      /*  avoid duplication */
      if((indr[n] < rAbeg) || (indr[n] >= rAend)) indr[n] = -1;
    }
    indc[0] = kp;
    ierr = MatSetValues(KPmat,nzone,indr,1,indc, val,INSERT_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = MatAssemblyBegin(KPmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(KPmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  delete [] indr;
  delete [] val;
};	// SegCond::setupKPmat

void SegCond::setupWGmat()
{
  /*  setup matrix to get avg. temp. of graphite zone (up to jheat-1) */
  int jgra = jheat-1;	// except helium gap
  double volgra = 0;
  for(int j=0; j < jgra; j++) volgra += vol[j];
  double *wgt = new double[jheat-1];
  for(int j=0; j < jgra; j++) {
    wgt[j] = vol[j]/volgra;
  }
  ierr = MatZeroEntries(WGmat);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  int *indr = new int[jgra];
  int indc[1];
  for(int kp=0; kp < npows; kp++) {
    /*  when all rows of a case is out side we skip  */
    int indfrom = kp*nzone;       int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    for(int j=0; j < jgra; j++) {
      indr[j] = kp*nzone+j;
      /*  avoid duplication */
      if((indr[j] < rAbeg) || (indr[j] >= rAend)) indr[j] = -1;
    }
    indc[0] = kp;
    ierr = MatSetValues(WGmat,jgra,indr,1,indc, wgt,INSERT_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = MatAssemblyBegin(WGmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(WGmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  delete [] indr;
  delete [] wgt;
};	// SegCond::setupWGmat

void SegCond::setHcoef(double hcoef_i)
{
  assert(hcoef != NULL);
  for(int kp=0; kp < npows; kp++) hcoef[kp] = hcoef_i;
  /*  setupAmat when hcoef is changed */
  do_setup = 1;
};	// SegCond::setHcoef (scaler)

void SegCond::setHcoef(Vec hcoefvec)
{
  assert(hcoef != NULL);
  Vec Hwork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(hcoefvec,&ctx,&Hwork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,hcoefvec,Hwork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,hcoefvec,Hwork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  double *h_vec;
  ierr = VecGetArray(Hwork,&h_vec);
  for(int kp=0; kp < npows; kp++) hcoef[kp] = h_vec[kp];
  ierr = VecRestoreArray(Hwork,&h_vec);

  do_setup = 1;
};	// SegCond::setHcoef (vector)

void SegCond::setCond(int nmats_i, double cond_i[], double rhocp_i[])
{
  /*  set same value for all case  */
  if(nmats_i != nmats) {
    printf("*** %s:%d mismatch in no. of materials, require %d but %d.\n",
        __FILE__,__LINE__,nmats,nmats_i);
    printf("    3 materils: graphite/gap/compact\n");
  }
  assert(nmats_i == nmats);
  for(int m=0; m < nmats; m++) {
    for(int kp=0; kp < nmats; kp++) {
      cond[kp*nmats+m] = cond_i[m];
      rhocp[kp*nmats+m] = rhocp_i[m];
    }
  }

  do_setup = 1;
};	// SegCond::setCond (scalar)

void SegCond::setCond(int nmats_i, Vec condvec[], Vec rhocpvec[])
{
  assert(nmats_i == nmats);
  for(int m=0; m < nmats; m++) {
    /*  confirm size of vectors */
    int szvec;
    ierr = VecGetSize(condvec[m],&szvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
    assert(szvec == npows);
    ierr = VecGetSize(rhocpvec[m],&szvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
    assert(szvec == npows);

    Vec condwork;
    VecScatter ctx_cond;
    VecScatterCreateToAll(condvec[m],&ctx_cond,&condwork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
    VecScatterBegin(ctx_cond,condvec[m],condwork,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx_cond,condvec[m],condwork,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx_cond);
    double *condkp;
    ierr = VecGetArray(condwork,&condkp);
    for(int kp=0; kp < npows; kp++) cond[kp*nmats+m] = condkp[kp];
    ierr = VecRestoreArray(condwork,&condkp);

    Vec rhocpwork;
    VecScatter ctx_rcp;
    VecScatterCreateToAll(rhocpvec[m],&ctx_rcp,&rhocpwork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
    VecScatterBegin(ctx_rcp,rhocpvec[m],rhocpwork,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(ctx_rcp,rhocpvec[m],rhocpwork,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterDestroy(&ctx_rcp);
    double *rhocpkp;
    ierr = VecGetArray(rhocpwork,&rhocpkp);
    for(int kp=0; kp < npows; kp++) rhocp[kp*nmats+m] = rhocpkp[kp];
    ierr = VecRestoreArray(rhocpwork,&rhocpkp);
  }	// for m

  do_setup = 1;
};      // SegCond::setCond (scalar)

void SegCond::setup(int npows_i)
{
  npows = npows_i;
  prepare();
  setupQvec(jheat);
};	// SegCond::setup

void SegCond::setupAmat()
{
  ierr = MatZeroEntries(Amat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSet(BCvec,0.0);	//

  for(int kp=0; kp < npows; kp++) {
    /*  when all rows of a case is out side we skip  */
    int indfrom = kp*nzone;       int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;

    double tauBulk = 0;
    for(int n=0; n < nzone; n++) {
      /* interface n-1:n  */
      double condK,tauL;
      if(n > 0) {
        condK = cond[kp*nmats+mat[n]];
        double condL = cond[kp*nmats+mat[n-1]];
        tauL = sl[n]/(dL[n]/condK+dR[n-1]/condL);
      }
      else {
        condK = cond[mat[n]];
        tauL = sl[n]*(condK/dL[n]+hcoef[kp]);
        tauBulk = tauL;
      }
      /* interface n:n+1  */
      double condR,tauR;
      if(n < nzone-1) {
        condR = cond[kp*nmats+mat[n+1]];
        tauR = sl[n+1]/(dR[n]/condK+dL[n+1]/condR);
      }
      else tauR = 0;
      int indr[1],indc[3];
      double A[3];
      indr[0] = kp*nzone+n;
      /* avoid multiple addition to a row  */
      if((indr[0] < rAbeg) || (indr[0] >= rAend)) indr[0] = -1;
      indc[0] = kp*nzone+n-1;	A[0] = -tauL;
      if(n < 1) indc[0] = -1;
      indc[1] = kp*nzone+n;	A[1] =  tauL+tauR;
      indc[2] = kp*nzone+n+1;	A[2] = -tauR;
      if(n+1 >= nzone) indc[2] = -1;
      ierr = MatSetValues(Amat,1,indr,3,indc, A,ADD_VALUES);
     
    }	// for n
    int indBC[1];
    indBC[0] = kp*nzone;
    if((indBC[0] >= rAbeg) && (indBC[0] < rAend)) {
      ierr = VecSetValues(BCvec,1,indBC,&tauBulk, ADD_VALUES);
    }
  }	// for kp
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyBegin(BCvec);
  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(BCvec);

  do_setup = 0;
};	// SegCond::setupAmat

void SegCond::setupQvec(int jheat)
{
  /* set volume for heat generation */
  ierr = VecSet(Qvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  int nheat = nzone-jheat;
  int *indr = new int[nheat];
  double *Q = new double[nheat];
  for(int kp=0; kp < npows; kp++) {
    int indfrom = kp*nzone;	int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    Qvol = 0.0;
    for(int j=0; j < nheat; j++) {
      indr[j] = kp*nzone + jheat+j;
      if((indr[j] < rAbeg) || (indr[j] >= rAend)) indr[j] = -1;
      /*  we use identical geometry  */
      Q[j] = vol[jheat+j];
      Qvol += Q[j];
    }
    ierr = VecSetValues(Qvec,nheat,indr,Q, INSERT_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }	// for kp
  ierr = VecAssemblyBegin(BCvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(BCvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  delete [] indr;
  delete [] Q;
};	// SegCond::setupQvec

void SegCond::setPfactor(double pfact_i)
{
  /*  multiplier to input power density  */
  pfact = pfact_i;
};	// SegCond::setPfactor

void SegCond::steady(double Tbulk, double qden)
{
  if(do_setup) setupAmat();
  ierr = VecCopy(BCvec,srcvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(srcvec,Tbulk);
  VecAXPY(srcvec,pfact*qden,Qvec);
  solver->solve(Amat,srcvec,Tvec);
  if(prlev > 1) VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
};	// SegCond::steady (scalar)

void SegCond::steady(Vec Tbulkvec, Vec qdenvec)
{
  if(do_setup) setupAmat();
  /*  apply case wise appplication of Tbulk and qden */
  ierr = MatMult(KPmat,Tbulkvec,workvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecPointwiseMult(srcvec,BCvec,workvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatMult(KPmat,qdenvec,work2vec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecPointwiseMult(workvec,work2vec,Qvec);
  ierr = VecAXPY(srcvec,pfact,workvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("%s:%d srcvec.\n",__FILE__,__LINE__);
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
    printf("%s:%d Amat.\n",__FILE__,__LINE__);
    MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
  }
  solver->solve(Amat,srcvec,Tvec);
  if(prlev > 1) {
    if(mpid==0) printf("%s:%d Tvec.\n",__FILE__,__LINE__);
    VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};      // SegCond::steady (vector)

void SegCond::start()
{
  /*  prepare for transient  */
  int maxdiag = ndim;
  int maxoff = ndim;
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &Mmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Mmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Mmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &TRmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(TRmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(TRmat); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyBegin(TRmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(TRmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);


  ierr = VecDuplicate(Tvec,&Tgam);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(Tvec,&Tnew);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(Tvec,&srcnew);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(Tvec,&rhsvec);   CHKERRABORT(PETSC_COMM_SELF,ierr);

  setupMmat();
};	//  SegCond::start

void SegCond::setupMmat()
{
  /*  setup rhocp  */
  ierr = MatZeroEntries(Mmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int kp=0; kp < npows; kp++) {
    int indfrom = kp*nzone;       int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;

    for(int n=0; n < nzone; n++) {
      int indr[1],indc[1];
      double M[1];
      indr[0] = kp*nzone + n;
      if((indr[0] < rAbeg) || (indr[0] >= rAend)) indr[0] = -1;
      indc[0] = kp*nzone + n;
      M[0] = rhocp[kp*nmats+mat[n]]*vol[n];
      ierr = MatSetValues(Mmat,1,indr,1,indc, M,ADD_VALUES);
    }
  }	// for kp
  ierr = MatAssemblyBegin(Mmat,MAT_FINAL_ASSEMBLY);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// SegCond::setupMmat

void SegCond::setLTE(double epsR_i, double epsA_i)
{
  /*  set LTE error estimator */
  epsR = epsR_i;
  epsA = epsA_i;
};      // SegCond::setLTE


double SegCond::step(double dt, double Tbulk, double qden)
{
  if(do_setup) setupAmat();
  /*  prepare new source */
  ierr = VecCopy(BCvec,srcnew);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(srcnew,Tbulk);
  VecAXPY(srcnew,qden,Qvec);
  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  return LTE;
};	// SegCond::step (scalar)

double SegCond::step(double dt, Vec Tbulkvec, Vec qdenvec)
{
  if(do_setup) setupAmat();
  /*  prepare new source */
  /*  apply case wise appplication of Tbulk and qden */
  ierr = MatMult(KPmat,Tbulkvec,workvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecPointwiseMult(srcnew,BCvec,workvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatMult(KPmat,qdenvec,work2vec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecPointwiseMult(workvec,work2vec,Qvec);
  ierr = VecAXPY(srcnew,pfact,workvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  return LTE;
};      // SegCond::step (vector)

void SegCond::march()
{
  ierr = VecCopy(Tnew,Tvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecCopy(srcnew,srcvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// SegCond::march

void SegCond::step_gam(double dt)
{
  /*  Gamma step  */
  /*  TRgam = 2*rhoM + gam*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,TRgam*dt,Amat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  solver->ksp_setup(TRmat);
  /*  rhs = 2*M*T_o - gam*dt*A*T_o + gam*dt*((1-gam)*p_o+gam*p_n)  */
  MatMult(Mmat,Tvec,rhsvec);
  VecScale(rhsvec,2.0);
  MatMult(Amat,Tvec,workvec);
  VecAXPY(rhsvec,-TRgam*dt,workvec);

  ierr = VecCopy(srcvec,workvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(workvec,(2-TRgam));
  VecAXPY(workvec,TRgam,srcnew);
  VecAXPY(rhsvec,TRgam*dt,workvec);

  solver->ksp_solve(rhsvec,Tgam);
  solver->ksp_destroy();
};      // SegCond::step_gam

void SegCond::step_new(double dt)
{
  /*  new step  */
  /*  TRnew = (2-gam)*M + (1-gam)*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0-TRgam);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,(1.0-TRgam)*dt,Amat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  solver->ksp_setup(TRmat);
  /*  rhs = gam^-1*(M*T_g - (1-gam)^2*M*T_o) - (1-gam)*dt*p_n  */
  MatMult(Mmat,Tgam,rhsvec);
  MatMult(Mmat,Tvec,workvec);
  VecAXPY(rhsvec,-(1-TRgam)*(1-TRgam),workvec);
  VecScale(rhsvec,1.0/TRgam);
  VecAXPY(rhsvec,(1-TRgam)*dt,srcnew);

  solver->ksp_solve(rhsvec,Tnew);
  solver->ksp_destroy();
};      // SegCond::step_new

double SegCond::getLTE(double dt)
{
  VecWAXPY(workvec,-1.0/(1-TRgam),Tgam,Tvec);
  VecScale(workvec,1.0/TRgam);
  VecAXPY(workvec,1.0/(1-TRgam),Tnew);
  double wnorm;
  VecNorm(workvec,NORM_2,&wnorm);
  /*  now prepare divider  */
  double rnorm;
#ifdef MMAT
  MatMult(Mmat,Tnew,rhsvec);
  VecNorm(rhsvec,NORM_2,&rnorm);
#else
  VecNorm(Tnew,NORM_2,&rnorm);
#endif
  double LTE = wnorm/(epsR*rnorm+epsA);
  return LTE;
};      // SegCond::getLTE

void SegCond::dumpT()
{
  if(mpid==0) printf("== Tvec\n");
  VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
};

double SegCond::getTpeakNew()
{
  int ppos;
  double Tpeak;
  VecMax(Tnew,&ppos,&Tpeak);
  return Tpeak;
};	//  SegCond::getTpeak

double SegCond::getTpeak()
{
  int ppos;
  double Tpeak;
  VecMax(Tvec,&ppos,&Tpeak);
  return Tpeak;
};      //  SegCond::getTpeak

double SegCond::getTpin()
{
  /*  get volume average temp of fuel pin */
  /*  we know that Qvec is volume of the power prod. zone */
  VecPointwiseMult(workvec,Tvec,Qvec);
  double wsum;
  VecSum(workvec,&wsum);
  return wsum/Qvol/npows;
};	// SegCond::getTpin

Vec SegCond::getTpinVec()
{
  /* returns vector similar to Tpin  */
  VecPointwiseMult(workvec,Tvec,Qvec);
  ierr = MatMultTranspose(KPmat,workvec,Tpinvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScale(Tpinvec,1.0/Qvol);
  return Tpinvec;
};	// Vec SegCond::getTpinVec

Vec SegCond::getTgraVec()
{
  /*  returns average temperature of graphite region */
  ierr = MatMultTranspose(WGmat,Tvec,Tgravec);
  return Tgravec;
};
