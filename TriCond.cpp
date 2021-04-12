/*
 * 	TriCond.cpp
 * 	spherical FVM problem
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "TriCond.h"

TriCond::TriCond()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  PI = 2*acos(0.0);
  pfact = 1.0;
  do_setup = 1;
  mat = NULL;
  vol = NULL;
  volbig = NULL;
  sl = NULL;
  dL = NULL;
  dR = NULL;
  cond = NULL;
  rhocp = NULL;
  solver = new Solver();
  solver->setTol(5000,1.0e-8,1.0e-8);

  Amat = NULL;
  Qvec = NULL;
  Qkvec = NULL;
  Qfvec = NULL;
  Tvec = NULL;
  srcvec = NULL;
  workvec = NULL;
  work2vec = NULL;
  cols = NULL;
  npows = 1;
  KPmat = NULL;
  dTvec = NULL;
  
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
};	// TriCond::TriCond

TriCond::~TriCond()
{
  delete [] mat;
  delete [] vol;
  delete [] volbig;
  delete [] sl;
  delete [] dL;
  delete [] dR;
  delete [] cond;
  delete [] rhocp;

  delete  solver;
  MatDestroy(&Amat);
  VecDestroy(&Qvec);
  VecDestroy(&Qkvec);
  VecDestroy(&Qfvec);
  delete [] cols;
  VecDestroy(&Tvec);
  VecDestroy(&srcvec);
  VecDestroy(&workvec);
  VecDestroy(&work2vec);
  MatDestroy(&KPmat);
  VecDestroy(&dTvec);

  MatDestroy(&Mmat);
  MatDestroy(&TRmat);
  VecDestroy(&rhsvec);
  VecDestroy(&Tgam);
  VecDestroy(&Tnew);
  VecDestroy(&srcnew);
};	// TriCond::~TriCond

void TriCond::setPrlev(int prl)
{
  prlev = prl;
};	// TriCond::setPrlev

void TriCond::setGeom()
{
  /*  arrange center to out  */
  double um = 1.0e-6;
  double Dker = 500.0*um;  int nker = 4;
  double tBuf = 100.0*um;  int nbuf = 4;
  double tiPyC = 40.0*um;  int nipy = 1;
  double thSiC = 35.0*um;  int nsic = 1;
  double toPyC = 40.0*um;  int nopy = 1;
  double packfr = 0.4;	   int nmat = 3;
  double rtriso = 0.5*Dker+tBuf+tiPyC+thSiC+toPyC;
  double vtriso = rtriso*rtriso*rtriso;
  double vFCM = vtriso/packfr;
  double rFCM = pow(vFCM,0.33333);
  double tMat = rFCM-rtriso;
  nzone = nker+nbuf+nipy+nsic+nopy+nmat;
  if(prlev) printf("rtriso=%.2lf rFCM=%.2lf %d zones\n",
	rtriso/um,rFCM/um,nzone);
  double *r = new double[nzone+1];
  mat = new int[nzone];
  r[0] = 0.0;
  int n=1;
  double dr = Dker/2/nker;
  for(int j=0; j < nker; j++) {
    mat[n-1] = 0;       // kernel
    r[n] = r[n-1] + dr;
    n++;
  }
  dr = tBuf/nbuf;
   for(int j=0; j < nbuf; j++) {
    mat[n-1] = 1;       // buffer
    r[n] = r[n-1] + dr;
    n++;
  }
  dr = tiPyC/nipy;
   for(int j=0; j < nipy; j++) {
    mat[n-1] = 2;       // iPyC
    r[n] = r[n-1] + dr;
    n++;
  }
  dr = thSiC/nsic;
   for(int j=0; j < nsic; j++) {
    mat[n-1] = 3;       // SiC
    r[n] = r[n-1] + dr;
    n++;
  }
  dr = toPyC/nopy;
   for(int j=0; j < nopy; j++) {
    mat[n-1] = 4;       // oPyC
    r[n] = r[n-1] + dr;
    n++;
  }
  dr = tMat/nmat;
   for(int j=0; j < nmat; j++) {
    mat[n-1] = 5;       // Matrix
    r[n] = r[n-1] + dr;
    n++;
  }
  assert(n == nzone+1);
  if(prlev) {
    printf("r:");
    for(int j=0; j < nzone; j++) {
      printf(" %.2f",r[j+1]/um);
    }
    printf("\n");
    printf("mat:");
    for(int j=0; j < nzone; j++) {
      printf("%5d  ",mat[j]);
    }
    printf("\n");
  }
  kernel = nker;
  double rK = r[nker];
  double rP = r[nzone];
  vKernel = 4*PI/3*rK*rK*rK;
  vParticle = 4*PI/3*rP*rP*rP;
  if(prlev) printf("rK=%.2lf rP=%.2lf\n",rK/um,rP/um);
  vol = new double[nzone];
  sl  = new double[nzone];
  dL  = new double[nzone];
  dR  = new double[nzone];
  double v0 = 0.0;
  for(int j=0; j < nzone; j++) {
    double vj = 4*PI/3*r[j+1]*r[j+1]*r[j+1];
    vol[j] = vj-v0;
    v0 = vj;
    sl[j] = 4*PI*r[j+1]*r[j+1];
    double dr = (r[j+1]-r[j]);
    dL[j] = 0.5*dr;
    dR[j] = 0.5*dr;
  }
  /*  introduce volbig for reasonable matrix  */
  volbig = new double[nzone];
  for(int j=0; j < nzone; j++) volbig[j] = vol[j]*1.0e+10;
  /*  enforce kernel region volume to zero  */
  for(int j=0; j < nker; j++) volbig[j] = 0.0;
  /*	assign conductivity */
  nmats = 6;
  cond = new double[nmats];
  rhocp = new double[nmats];
  cond[0] = 3.7;	rhocp[0] = 1.0e+6;	// kernel
  cond[1] = 0.5;	rhocp[1] = 1.0e+6;	// buf
  cond[2] = 4.0;	rhocp[2] = 1.0e+6;	// iPyC
  cond[3] = 16.0;	rhocp[3] = 1.0e+6;	// SiC
  cond[4] = 4.0;        rhocp[4] = 1.0e+6;	// oPyC
  cond[5] = 5.0;	rhocp[5] = 1.0e+6;	// matrix
  delete [] r;
};	// TriCond::setGeom

void TriCond::prepare()
{
  /*  for mutliple case  */
  ndim = npows*nzone;
  assert(ndim > 0);
  assert(Amat == NULL);
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
  VecDuplicate(Tvec,&Qvec);
  VecDuplicate(Tvec,&Qkvec);
  VecDuplicate(Tvec,&Qfvec);
  VecDuplicate(Tvec,&srcvec);
  VecDuplicate(Tvec,&workvec);
  VecDuplicate(Tvec,&work2vec);

  setupQvec();

  /* unit cols */
  cols = new int[ndim];

  ierr = MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &KPmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(KPmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(KPmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&dTvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(dTvec,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(dTvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);

  setupKPmat();

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
    }   // for kp
    delete [] cond_old;
    delete [] rhocp_old;
  }
};	// TriCond::prepare

void TriCond::setupKPmat()
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
};	// TriCond::setupKPmat

void TriCond::setup(int npows_i)
{
  npows = npows_i;
  prepare();
};	// TriCond::setup

void TriCond::setCond(int nmats_i, double cond_i[], double rhocp_i[])
{
  /*  set same value for all case  */
  if(nmats_i != nmats) {
    printf("*** %s:%d mismatch in no. of materials, require %d but %d.\n",
	__FILE__,__LINE__,nmats,nmats_i);
    printf("    5 materils: kernel/buffer/PyC/SiC/Matrix\n");
  }
  assert(nmats_i == nmats);
  for(int m=0; m < nmats; m++) {
    for(int kp=0; kp < nmats; kp++) {
      cond[kp*nmats+m] = cond_i[m];
      rhocp[kp*nmats+m] = rhocp_i[m];
    }
  }

  do_setup = 1;
};      // TriCond::setCond (scalar)

void TriCond::setCond(int nmats_i, Vec condvec[], Vec rhocpvec[])
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
  }     // for m

  do_setup = 1;
};      // TriCond::setCond (scalar)

void TriCond::setupQvec()
{
  /*  make volume vectors  for net zero power production */
  /*  Qfvec :=  1.0/vKer (for kernel region) - 1.0/vPart  */
  VecSet(Qvec,0.0);
  VecSet(Qfvec,0.0);
  VecSet(Qfvec,0.0);
  double *frac = new double[nzone];
  for(int n=0; n < nzone; n++) {
    frac[n] = -vol[n]/vParticle;
    if(n < kernel) frac[n] += vol[n]/vKernel;
  }

  int *indr = new int[nzone];
  for(int kp=0; kp < npows; kp++) {
    int indfrom = kp*nzone;     int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    for(int n=0; n < nzone; n++) {
      indr[n] = kp*nzone + n;
      if((indr[n] < rAbeg) || (indr[n] >= rAend)) indr[n] = -1;
    }
    VecSetValues(Qvec, nzone, indr,vol,INSERT_VALUES);
    VecSetValues(Qkvec,kernel,indr,vol,INSERT_VALUES);
    VecSetValues(Qfvec,nzone, indr,frac,INSERT_VALUES);
  }	// for kp
  VecAssemblyBegin(Qvec);
  VecAssemblyBegin(Qkvec);
  VecAssemblyBegin(Qfvec);
  VecAssemblyEnd(Qvec);
  VecAssemblyEnd(Qkvec);
  VecAssemblyEnd(Qfvec);
  delete [] frac;
  delete [] indr;
};	// TriCond::setupQvec

void TriCond::setupAmat()
{
  /*  we apply Reflective condition at both side */
  /*  we apply average zero temperature condition */
  ierr = MatZeroEntries(Amat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  for(int kp=0; kp < npows; kp++) {
    int indfrom = kp*nzone;     int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    /*  sl is right side surface area  */
    for(int n=0; n < nzone; n++) {
      /* interface n-1:n  */
      double condK = cond[kp*nmats+mat[n]];
      double tauL = 0;
      if(n > 0) {
        double condL = cond[kp*nmats+mat[n-1]];
        tauL = sl[n-1]/(dL[n]/condK+dR[n-1]/condL);
      }
      /* interface n:n+1  */
      double condR,tauR;
      if(n < nzone-1) {
        condR = cond[kp*nmats+mat[n+1]];
        tauR = sl[n]/(dR[n]/condK+dL[n+1]/condR);
      }
      else tauR = 0;
      int indr[1],indc[3];
      double A[3];
      indr[0] = kp*nzone+n;
      if((indr[0] < rAbeg) || (indr[0] >= rAend)) indr[0] = -1;
      indc[0] = kp*nzone+n-1;      A[0] = -tauL;
      if(n < 1) indc[0] = -1;
      indc[1] = kp*nzone+n;        A[1] =  tauL+tauR;
      indc[2] = kp*nzone+n+1;      A[2] = -tauR;
      if(n+1 >= nzone) indc[2] = -1;
      ierr = MatSetValues(Amat,1,indr,3,indc, A,ADD_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
      if(prlev > 1) printf("n%d: %.2le %.2le\n",n,tauL,tauR);
    }	// for n
  }	// for kp
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  do_setup = 0;
};	// TriCond:;setupAmat

void TriCond::setPfactor(double pfact_i)
{
  /*  set multiplier to input pker  */
  pfact = pfact_i;
};	// TriCond::setPfactor(

void TriCond::steady(double pker)
{
  if(do_setup) setupAmat();
  /*  split qker to kernel and others  */
  ierr = VecCopy(Qfvec,srcvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(srcvec,pfact*pker);

  zeroTavgMat(Amat);
  zeroTavgRhs(srcvec);
  if(prlev > 1) VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  int kiter = solver->solve(Amat,srcvec,Tvec);
  if(kiter < 0) {
     printf("Amat:");
     MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
     printf("srcvec:");
     VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
     exit(0);
  }
  if(prlev > 1) {
    printf("Tvec:\n");
    VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};	// TriCond::steady (scalar)

void TriCond::steady(Vec pkervec)
{
  if(do_setup) setupAmat();
  /*  split qker to kernel and others  */
  /*  src = pker*Qfvec  */
  ierr = MatMult(KPmat,pkervec,workvec);
  VecScale(workvec,pfact);
  ierr = VecPointwiseMult(srcvec,Qfvec,workvec);

  zeroTavgMat(Amat);
  zeroTavgRhs(srcvec);
  if(prlev > 1) {
    MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  int kiter = solver->solve(Amat,srcvec,Tvec);

  if(kiter < 0) {
     printf("Amat:");
     MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
     printf("srcvec:");
     VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
     printf("*** nzone=%d npows=%d\n",nzone,npows);
     exit(0);
  }

  if(prlev > 1) {
    printf("Tvec:\n");
    VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};      // TriCond::steady (vector)


void TriCond::zeroTavgMat(Mat Amat)
{
  /*  enforce zero Tavg  */
  int lastrow[1];
  for(int kp=0; kp < npows; kp++) {
    lastrow[0] = kp*nzone+nzone-1;
    for(int n=0; n < nzone; n++) cols[n] = kp*nzone+n;
    ierr = MatSetValues(Amat,1,lastrow,nzone,cols,volbig, INSERT_VALUES);
  	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// TriCond::zeroTavgMat

void TriCond::zeroTavgRhs(Vec rhs)
{
  /*  enforce zero Tavg  */
  for(int kp=0; kp < npows; kp++) {
    int lastrow = kp*nzone+nzone-1;
    ierr = VecSetValue(rhs,lastrow,0.0,INSERT_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = VecAssemblyBegin(rhs);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(rhs);	CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // TriCond::zeroTavgRhs


void TriCond::start()
{
  /*  prepare for transient  */
  assert(Mmat==NULL);
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
};      //  TriCond::start

void TriCond::setupMmat()
{
  /*  setup rhocp  */
  ierr = MatZeroEntries(Mmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  for(int kp=0; kp < npows; kp++) {
    int indfrom = kp*nzone;     int indto = (kp+1)*nzone;
    if((indfrom > rAend) || (indto < rAbeg)) continue;
    for(int n=0; n < nzone; n++) {
      int indr[1],indc[1];
      double M[1];
      indr[0] = kp*nzone + n;
      if((indr[0] < rAbeg) || (indr[0] >= rAend)) indr[0] = -1;
      indc[0] = kp*nzone + n;
      M[0] = rhocp[mat[n]]*vol[n];
      ierr = MatSetValues(Mmat,1,indr,1,indc, M,ADD_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  }
  ierr = MatAssemblyBegin(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // TriCond::setupMmat

void TriCond::setLTE(double epsR_i, double epsA_i)
{
  /*  set LTE error estimator */
  epsR = epsR_i;
  epsA = epsA_i;
};      // TriCond::setLTE


double TriCond::step(double dt, double pker_i)
{
  if(do_setup) setupAmat();
  double pker = pfact*pker_i;
  /*  split qker to kernel and others  */
  double qden = pker/vKernel;
  double qpden = pker/vParticle;
  double qkden = qden-qpden;
  if(prlev > 1) printf(" pker=%.2le qkden=%.2le qpden=%.2le\n",
	pker,qkden,qpden);
  ierr = VecCopy(Qkvec,srcnew);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(srcnew,qkden);
  VecAXPY(srcnew,-qpden,Qvec);
  VecSetValue(srcnew,ndim-1,0.0,INSERT_VALUES);
  VecAssemblyBegin(srcnew);
  VecAssemblyEnd(srcnew);

  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  return LTE;
};      // TriCond::step (scalar)

double TriCond::step(double dt, Vec pkervec)
{
  if(do_setup) setupAmat();
  ierr = MatMult(KPmat,pkervec,workvec);
  VecScale(workvec,pfact);
  ierr = VecPointwiseMult(srcnew,Qfvec,workvec);

  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  return LTE;
};      // TriCond::step (vector)

void TriCond::march()
{
  ierr = VecCopy(Tnew,Tvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecCopy(srcnew,srcvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // TriCond::march

void TriCond::step_gam(double dt)
{
  /*  Gamma step  */
  /*  TRgam = 2*rhoM + gam*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,TRgam*dt,Amat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  zeroTavgMat(TRmat);
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

  zeroTavgRhs(rhsvec);
  solver->ksp_solve(rhsvec,Tgam);
  solver->ksp_destroy();
};      // TriCond::step_gam

void TriCond::step_new(double dt)
{
  /*  new step  */
  /*  TRnew = (2-gam)*M + (1-gam)*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0-TRgam);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,(1.0-TRgam)*dt,Amat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  zeroTavgMat(TRmat);
  solver->ksp_setup(TRmat);
  /*  rhs = gam^-1*(M*T_g - (1-gam)^2*M*T_o) - (1-gam)*dt*p_n  */
  MatMult(Mmat,Tgam,rhsvec);
  MatMult(Mmat,Tvec,workvec);
  VecAXPY(rhsvec,-(1-TRgam)*(1-TRgam),workvec);
  VecScale(rhsvec,1.0/TRgam);
  VecAXPY(rhsvec,(1-TRgam)*dt,srcnew);

  zeroTavgRhs(rhsvec);
  solver->ksp_solve(rhsvec,Tnew);
  solver->ksp_destroy();
};      // TriCond::step_new

double TriCond::getLTE(double dt)
{
  VecWAXPY(workvec,-1.0/(1-TRgam),Tgam,Tvec);
  VecScale(workvec,1.0/TRgam);
  VecAXPY(workvec,1.0/(1-TRgam),Tnew);
  double wnorm;
  VecNorm(workvec,NORM_1,&wnorm);
  /*  now prepare divider  */
  /*  Multiply to Mmat is too strict  */
  double rnorm;
  VecNorm(Tnew,NORM_1,&rnorm);
  double LTE = wnorm/(epsR*rnorm+epsA);
  if(prlev) {
    int ppos;
    double Tpeak;
    VecMax(Tnew,&ppos,&Tpeak);
    printf("LTE=%.3le %.3le %.2le  %.2lf\n",wnorm,rnorm,LTE, Tpeak);
  }
  return LTE;
};      // TriCond::getLTE

void TriCond::dumpT()
{
  if(mpid==0) printf("== Tvec\n");
  VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);
};

double TriCond::getTker()
{
  /*  average of kernel  */
  VecPointwiseMult(workvec,Tvec,Qkvec);
  double Tksum;
  VecSum(workvec,&Tksum);
  return Tksum/(npows*vKernel);
};	// TriCond::getTker

Vec TriCond::getDTker()
{
  /*  get dTkernel  vector  */
  VecPointwiseMult(workvec,Tvec,Qkvec);
  ierr = MatMultTranspose(KPmat,workvec,dTvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecScale(dTvec,1.0/vKernel);
  return dTvec;
};	//

double TriCond::getTavg()
{
  /*  average of problem  */
  ierr = VecPointwiseMult(workvec,Tvec,Qvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  double Tksum;
  VecSum(workvec,&Tksum);
  return Tksum/vParticle;
};      // TriCond::getTavg

