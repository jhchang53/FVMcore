/*
 * 	NeutFVM.cpp
 * 	Neutron solver using FVM
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "NeutFVM.h"

NeutFVM::NeutFVM()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  do_setup = 1;
  do_power = 1;
  Amat = NULL;
  Bmat = NULL;
  Pmat = NULL;
  phivec = NULL;
  powvec = NULL;
  volvec = NULL;
  powden = NULL;
  eig_shift = 0.9;
  Tm = NULL;
  Tf = NULL;
  nrods = 0;
  sigrod = NULL;
  rodpos = NULL;
  ksp_itmax = 1000;
  ksp_rtol = 1.0e-8;
  ksp_atol = 1.0e-9;
  gvel = NULL;
  srcvec = NULL;
  fissvec = NULL;
  solver = new Solver();
  /*  transient */
  TRgam = 2.0-sqrt(2.0);
  TRC2 =  2*fabs(((-3*TRgam+4)*TRgam-2)/(12*(2-TRgam)));
  epsR = 1.0e-5;
  epsA = 1.0e-5;

  Mmat = NULL;
  ATmat = NULL;
  TRmat = NULL;
  Smat = NULL;
  phigam = NULL;
  phinew = NULL;
};	// NeutFVM::NeutFVM

NeutFVM::~NeutFVM()
{
  MatDestroy(&Amat);
  MatDestroy(&Bmat);
  MatDestroy(&Pmat);
  VecDestroy(&phivec);
  VecDestroy(&powvec);
  VecDestroy(&volvec);
  VecDestroy(&powden);
  delete [] Tm;
  delete [] Tf;
  delete [] sigrod;
  delete [] rodpos;
  delete [] gvel;
  VecDestroy(&srcvec);
  VecDestroy(&fissvec);
  delete solver;
  MatDestroy(&Mmat);
  MatDestroy(&ATmat);
  MatDestroy(&TRmat);
  MatDestroy(&Smat);
  VecDestroy(&phigam);
  VecDestroy(&phinew);
};	// NeutFVM::~NeutFVM

void NeutFVM::setPrlev(int prl)
{
  prlev = prl;
};	// NeutFVM::setPrlev

void NeutFVM::setGeom(Geom *geom_i)
{
  geom = geom_i;
  ntri = geom->getNtri();
  nz   = geom->getNz();
  nvols = ntri*nz;
  npows = geom->getNpows();
};	// NeutFVM::setGeom

void NeutFVM::setXS(XS *xs_i, double Tmref_i, double Tfref_i)
{
  xs = xs_i;
  Tmref = Tmref_i;
  Tfref = Tfref_i;
};	// NeutFVM::setXS

void NeutFVM::setTemp(double *Tmod, double *Tfuel)
{
  assert(Tm != NULL);
  assert(Tfuel != NULL);
  for(int tz=0; tz < nvols; tz++) {
    Tm[tz] = Tmod[tz];
    Tf[tz] = Tfuel[tz];
  }
  do_setup = 1;  // setup required
  do_power = 1;
};      // NeutFVM::setTemp

void NeutFVM::setRodXS(int nrods_i, double sigrod_i[])
{
  if(nrods < 1) {
    nrods = nrods_i;
    sigrod = new double[nrods];
    rodpos = new double[nrods];
  }
  else assert(nrods_i == nrods);
  for(int r=0; r < nrods; r++) sigrod[r] = sigrod_i[r];
};      // NeutFVM::setRodXS

void NeutFVM::setRodPos(int nrods_i, double rodpos_i[])
{
  assert(nrods == nrods_i);
  for(int r=0; r < nrods; r++) rodpos[r] = rodpos_i[r];
  do_setup = 1;
  do_power = 1;
};      // NeutFVM::setRodPo


void NeutFVM::prepare(int ngrps_i)
{
  ngrps = ngrps_i;
  ndim = ngrps*nvols;

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

  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &Bmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Bmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Bmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &ATmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(ATmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(ATmat); CHKERRABORT(PETSC_COMM_SELF,ierr);


  ierr = VecCreate(PETSC_COMM_WORLD,&phivec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(phivec,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(phivec);     CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecDuplicate(phivec,&srcvec);

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
  ierr = MatAssemblyEnd(TRmat,MAT_FINAL_ASSEMBLY);

  ierr = VecDuplicate(phivec,&rhsvec);
  ierr = VecDuplicate(phivec,&workvec);
  ierr = VecDuplicate(phivec,&phigam);
  ierr = VecDuplicate(phivec,&phinew);

  /*  temperatures for XS set default to Tref  */
  Tm  = new double[nvols];
  Tf = new double[nvols];
  for(int n=0; n < nvols; n++) {
    Tm[n]  = Tmref;
    Tf[n] = Tfref;
  }

  /*  power matrix  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &Pmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Pmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Pmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&powvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(powvec,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(powvec);     CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecDuplicate(powvec,&fissvec);
  ierr = VecDuplicate(powvec,&volvec);
  ierr = VecDuplicate(powvec,&powden);

  /* inverse fissvec to phivec  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &Smat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Smat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Smat); CHKERRABORT(PETSC_COMM_SELF,ierr);
  setupSmat();
  makeVolvec();
};	// NeutFVM::prepare

void NeutFVM::setupAB()
{
  setupABmat();
  setupPmat(1.0);
};

void NeutFVM::setupABmat()
{
  ierr = MatZeroEntries(Amat);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatZeroEntries(Bmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  double *dz = geom->Dz();
  double *zpos = geom->Zpos();
  int *crod2d = geom->Crod2D();

  for(int t=0; t < ntri; t++) {
    int r = crod2d[t]-1;	// crod2d is offseted by 1
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      /* *** This logic NEED modification for multi group  *** */
      if((tz < rAbeg) || (tz >= rAend)) continue;
      int indr[1];
      indr[0] = tz;
      /* 5 faces for a tri-prism  */
      int tb[5];
      double area[5],dist[5];
      double vol = geom->getFVM(t,z,tb,area,dist);
      if(prlev > 1) {
        printf("t=%d z=%d  A=%.3lf:",t,z,vol);
        for(int e=0; e < 5; e++) printf(" %d(%.1le,%.1le)",
	  tb[e],area[e],dist[e]);
        printf("\n");
      }
      double sig[3];
      xs->getSig(tz,Tm[tz],Tf[tz], sig);
      double sigtrK = sig[0]*100;	// into SI unit
      double sigabK = sig[1]*100;
      double signfK = sig[2]*100;
      double tauKL,sigtrL;
      int indc[6];
      double A[6];
      indc[0] = tz;
      A[0] = sigabK*vol;
      for(int i=1; i < 6; i++) {
        indc[i] = -1;
        A[i] = 0.0;
      }
      for(int e=0; e < 5; e++) {
        int etb = tb[e];
        if(etb >= 0) {
          double distL = geom->getDistL(etb,tz);
          sigtrL = xs->getSigtr(etb,Tm[tz],Tf[tz]);
          tauKL = area[e]/(3*sigtrK*dist[e]+3*sigtrL*distL);
        }
        else if(etb == -2) {	// zero BC
          tauKL = area[e]/(3*sigtrK*dist[e]);
        }
        else tauKL = 0.0;	// reflective BC
        indc[1+e] = etb;
        A[1+e] = -tauKL;
        A[0] += tauKL;
      }	// for e
      /*  check if control rod exist at this triangle  */
      if(r >= 0) {
        if(rodpos[r] > zpos[z]) {
          double sigr =  sigrod[r]*100.0;
          double zins = rodpos[r]-zpos[z];
          double zweit = zins/dz[z];
          if(zweit > 1.0) zweit = 1.0;
          A[0] += zweit*sigr*vol;
        }
      }

      if(prlev > 1) {
        printf("A %d:",tz);
        for(int i=0; i < 6; i++) printf(" %.2le(%d)",A[i],indc[i]);
        printf("\n");
      }
      ierr = MatSetValues(Amat,1,indr,6,indc, A,ADD_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
      int indb[1];
      double B[1];
      indb[0] = tz;
      B[0] = signfK*vol;
      ierr = MatSetValues(Bmat,1,indr,1,indb, B,ADD_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
    }	// for z
  }	// for t
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyBegin(Bmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Bmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  if(prlev > 1) {
    printf("Amat:\n");
    MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
  }
  do_setup = 0;
};	//  NeutFVM::setupAmat

void NeutFVM::setupPmat(double fact)
{
  ierr = MatZeroEntries(Pmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *chtri = geom->Chtri();
  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();

  int kf = 0;
  for(int t=0; t < ntri; t++) {
    int fisflag = 0;
    int c = chtri[t];
    if((c >= 0) && ispp[c]) fisflag = 1;
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      if(fisflag && ispowz[z]) {
       /*  for counter kf  */
       if((tz >= rAbeg) && (tz < rAend)) {
         int indr[1];
         indr[0] = tz;
         double signf = xs->getSignf(tz,Tm[tz],Tf[tz]) * 100.0*fact;
         double vol = geom->getVol(t,z);
         int indc[1];
         double P[1];
         indc[0] = kf;
         P[0] = signf*vol;
         ierr = MatSetValues(Pmat,1,indr,1,indc,P, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
       }
       kf++;
      }
    }	// for z
  }	// for t
  assert(kf == npows);
  ierr = MatAssemblyBegin(Pmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Pmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// NeutFVM::setupPmat

void NeutFVM::setupSmat()
{
  /*  prepare dsrc to flux */
  ierr = MatZeroEntries(Smat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *chtri = geom->Chtri();
  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();

  int kf = 0;
  for(int t=0; t < ntri; t++) {
    int fisflag = 0;
    int c = chtri[t];
    if((c >= 0) && ispp[c]) fisflag = 1;
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      if(fisflag && ispowz[z]) {
       if((tz >= rAbeg) && (tz < rAend)) {
         int indr[1];
         indr[0] = tz;
         int indc[1];
         double S[1];
         indc[0] = kf;
         double vol = geom->getVol(t,z);
         S[0] = 1.0;
         ierr = MatSetValues(Smat,1,indr,1,indc,S, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
       }
       kf++;
      }
    }   // for z
  }     // for t
  assert(kf == npows);
  ierr = MatAssemblyBegin(Smat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Smat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("Smat:\n");
    MatView(Smat,PETSC_VIEWER_STDOUT_WORLD);
  }
};      // NeutFVM::setupSmat

void NeutFVM::makeVolvec()
{
  /*  prepare volume vector for power density calculation */
  ierr = VecSet(volvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *chtri = geom->Chtri();
  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();

  int kf = 0;
  for(int t=0; t < ntri; t++) {
    int fisflag = 0;
    int c = chtri[t];
    if((c >= 0) && ispp[c]) fisflag = 1;
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      if(fisflag && ispowz[z]) {
        int indr[1];
        indr[0] = kf;
        double V[1];
        double vol = geom->getVol(t,z);
        V[0] = vol;
        ierr = VecSetValues(volvec,1,indr, V,INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
        kf++;
      }
    }   // for z
  }     // for t
  assert(kf == npows);
  ierr = VecAssemblyBegin(volvec); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(volvec); CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("volvec:\n");
    VecView(volvec,PETSC_VIEWER_STDOUT_WORLD);
    if(prlev > 5) exit(0);
  }
};	//  NeutFVM::makeVolvec

void NeutFVM::setEigShift(double shift)
{
  eig_shift = shift;
};

double NeutFVM::solveEigen()
{
  if(do_setup) {
    setupAB();
  }
  eigen = new EigenSIpower();
  eigen->setKSPtol(2000,1.0e-10);
  eigen->setEigtol(2000,1.0e-10);
  eigen->setProblem(Amat,Bmat,eig_shift);
  VecSet(phivec,1.0);
  eigen->guess(phivec);

  double keff = eigen->solve();
  ierr = VecCopy(eigen->getPhi(),phivec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  do_power = 1;
  return keff;
};      // FVMn::solveEigen


/*  transient routines  */
Vec NeutFVM::makePower(double totpow)
{
  /* normalize flux and power to give totpow  */
  ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
  double ptotal;
  VecSum(powvec,&ptotal);
  double factor = totpow/ptotal;
  VecScale(phivec,factor);
  VecScale(powvec,factor);
  do_power = 0;
  return powvec;
};      //  Neut::makePower

void NeutFVM::printPow()
{
  /*  scatter to local storage for print  */
  /*  it is collective action so ...  */
  getPowden();

  Vec Pwork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(powden,&ctx,&Pwork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,powden,Pwork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,powden,Pwork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);
  if(mpid != 0) return;
  /*  we scattered vector to all local processors */
  /*  we know that T is in ntri*nz arrangement  */
  double *paxial = new double[nz];
  for(int z=0; z < nz; z++) paxial[z] = 0.0;
  int npchan = 0;
  PetscScalar *P;
  ierr = VecGetArray(Pwork,&P);
  int nchan = geom->getNchan();
  int *ispp = geom->Ispp();
  int *ispowz = geom->Ispowz();
  int kf = 0;
  for(int c=0; c < nchan; c++) {
    if(ispp[c]) {
      printf("p%d:",c);
      for(int z=0; z < nz; z++) {
        if(ispowz[z]) {
          printf(" %.2le",P[kf]);
          paxial[z] += P[kf];
          kf++;
        }
      }
      printf("\n");
      npchan++;
    }
  }     // for c
  VecRestoreArray(Pwork,&P);
  VecDestroy(&Pwork);
  printf("avg:");
  for(int z=0; z < nz; z++) {
    if(paxial[z] > 0) printf(" %.2le",paxial[z]/npchan);
  }
  printf("\n");
};      // NeutFVM::printPow

/*  transient routines  */


void NeutFVM::start(double totpow, double keff, double beta_i, double gvel_i[])
{
  effk = keff;
  beta = beta_i;
  gvel = new double[ngrps];
  for(int g=0; g < ngrps; g++) gvel[g] = gvel_i[g];
  setupMmat();
  setupTransAmat((1-beta)/effk);
  setupPmat(1.0/effk);
  do_power = 1;
  do_setup =1;
  makePower(totpow);
  makeCinit();
};	//  NeutFVM::start

void NeutFVM::setupMmat()
{
  /*  make diagonal 1/v matrix for transient  */
  ierr = MatZeroEntries(Mmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);


  for(int t=0; t < ntri; t++) {
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      if((tz < rAbeg) || (tz >= rAend)) continue;
      int indr[1],indc[1];
      indr[0] = tz;
      indc[0] = tz;
      double vol = geom->getVol(t,z);
      double Vinv[1];
      Vinv[0] = vol/gvel[0];
      ierr = MatSetValues(Mmat,1,indr,1,indc,Vinv, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
    }   // for z
  }     // for t
  ierr = MatAssemblyBegin(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// NeutFVM::setupMmat


void NeutFVM::setupTransAmat(double fact)
{
  /*  setup  A = -D+sigab-(1-beta)*signf/effk   */
  /*	dialect of setupAmat  */
  ierr = MatZeroEntries(ATmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  double *dz = geom->Dz();
  double *zpos = geom->Zpos();
  int *crod2d = geom->Crod2D();

  for(int t=0; t < ntri; t++) {
    int r = crod2d[t]-1;        // crod2d is offseted by 1
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;
      if((tz < rAbeg) || (tz >= rAend)) continue;
      int indr[1];
      indr[0] = tz;
      /* 5 faces for a tri-prism  */
      int tb[5];
      double area[5],dist[5];
      double vol = geom->getFVM(t,z,tb,area,dist);
      if(prlev > 1) {
        printf("t=%d z=%d V=%.2le :",t,z,vol);
        for(int e=0; e < 5; e++) printf(" %d(%.1le,%.1le)",
          tb[e],area[e],dist[e]);
        printf("\n");
      }
      double sig[3];
      xs->getSig(tz,Tm[tz],Tf[tz], sig);
      double sigtrK = sig[0]*100;       // into SI unit
      double sigabK = sig[1]*100;
      double signfK = sig[2]*100;
      double tauKL,sigtrL;
      int indc[6];
      double A[6];
      indc[0] = tz;
      A[0] = (sigabK - fact*signfK)*vol;
      for(int i=1; i < 6; i++) {
        indc[i] = -1;
        A[i] = 0.0;
      }
      for(int e=0; e < 5; e++) {
        int etb = tb[e];
        if(etb >= 0) {
          double distL = geom->getDistL(etb,tz);
          sigtrL = xs->getSigtr(etb,Tm[tz],Tf[tz]);
          tauKL = area[e]/(3*sigtrK*dist[e]+3*sigtrL*distL);
        }
        else if(etb == -2) {    // zero BC
          tauKL = area[e]/(3*sigtrK*dist[e]);
        }
        else tauKL = 0.0;       // reflective BC
        indc[1+e] = etb;
        A[1+e] = -tauKL;
        A[0] += tauKL;
      } // for e
      /*  check if control rod exist at this triangle  */
      if(r >= 0) {
        if(rodpos[r] > zpos[z]) {
          double sigr =  sigrod[r]*100.0;
          double zins = rodpos[r]-zpos[z];
          double zweit = zins/dz[z];
          if(zweit > 1.0) zweit = 1.0;
          A[0] += zweit*sigr*vol;
        }
      }

      if(prlev > 1) {
        printf("A %d:",tz);
        for(int i=0; i < 6; i++) printf(" %.2le(%d)",A[i],indc[i]);
        printf("\n");
      }
      ierr = MatSetValues(ATmat,1,indr,6,indc, A,ADD_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
    }   // for z
  }     // for t
  ierr = MatAssemblyBegin(ATmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = MatAssemblyEnd(ATmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  if(prlev > 1) {
    printf("ATmat:\n");
    MatView(ATmat,PETSC_VIEWER_STDOUT_WORLD);
  }
  do_setup = 0;
};      //  NeutFVM::setupTransAmat

void NeutFVM::makeCinit()
{
  /*  compute steady C = (-D +siga - (1-beta)*signf/effk) * phi   */
  ierr = MatMult(ATmat,phivec,srcvec);	 CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("cinit src\n");
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};	// NeutFVM::makeCinit

void NeutFVM::setLTE(double epsR_i, double epsA_i)
{
  /*  parameter for LTE estimation  */
  epsR = epsR_i;
  epsA = epsA_i;
};	// NeutFVM::setLTE

double NeutFVM::step(double dt)
{
  if(do_setup) {
    setupTransAmat((1-beta)/effk);
    setupPmat(1.0/effk);
    do_setup = 0;
  }
#ifdef CRANK
  /*  simple Crank-Nicholson  */
  /*  solve   (M+0.5*dt*G) * phinew = M*phio - 0.5*dt*A *phio + dt*src  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,0.5*dt,ATmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  solver->ksp_setup(TRmat);
  /*  rhs = M * phivec - 0.5*dt*A*phivec + dt*src  */
  ierr = MatMult(Mmat,phivec,rhsvec);
  ierr = MatMult(ATmat,phivec,workvec);
  ierr = VecAXPY(rhsvec,-0.5*dt,workvec);
  ierr = VecAXPY(rhsvec,dt,srcvec);
  solver->ksp_solve(rhsvec,phinew);
  solver->ksp_destroy();
  /*  find difference  */
  ierr = VecWAXPY(workvec,-1.0,phivec,phinew);
  ierr = VecPointwiseDivide(rhsvec,workvec,phivec);
  int perr;
  double rerr;
  ierr = VecMax(rhsvec,&perr,&rerr);
  CHKERRABORT(PETSC_COMM_SELF,ierr);
  do_power = 1; // invalidate powvec
  return fabs(rerr);
#else
  /*  TRBDF  */
  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  do_power = 1;
  return LTE;
#endif
};	//  NeutFVM::step

void NeutFVM::march()
{
  ierr = VecCopy(phinew,phivec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// NeutFVM::march

void NeutFVM::step_gam(double dt)
{
  /*  Gamma step : we use ATmat for transient */
  /*  TRgam = 2*rhoM + gam*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,TRgam*dt,ATmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  solver->ksp_setup(TRmat);
  /*  rhs = 2*M*phi_o - gam*dt*A*phi_o + 2*gam*dt*src  */
  MatMult(Mmat,phivec,rhsvec);
  VecScale(rhsvec,2.0);
  MatMult(ATmat,phivec,workvec);
  VecAXPY(rhsvec,-TRgam*dt,workvec);

  VecAXPY(rhsvec,2*TRgam*dt,srcvec);

  solver->ksp_solve(rhsvec,phigam);
  solver->ksp_destroy();
};      // NeutFVM::step_gam

void NeutFVM::step_new(double dt)
{
  /*  new step   we use ATmat for transient */
  /*  TRnew = (2-gam)*M + (1-gam)*dt*A  */
  ierr = MatCopy(Mmat,TRmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatScale(TRmat,2.0-TRgam);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAXPY(TRmat,(1.0-TRgam)*dt,ATmat,DIFFERENT_NONZERO_PATTERN);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  solver->ksp_setup(TRmat);
  /*  rhs = gam^-1*(M*phi_g - (1-gam)^2*M*phi_o) - (1-gam)*dt*src  */
  MatMult(Mmat,phigam,rhsvec);
  MatMult(Mmat,phivec,workvec);
  VecAXPY(rhsvec,-(1-TRgam)*(1-TRgam),workvec);
  VecScale(rhsvec,1.0/TRgam);
  VecAXPY(rhsvec,(1-TRgam)*dt,srcvec);

  solver->ksp_solve(rhsvec,phinew);
  solver->ksp_destroy();
};      // NeutFVM::step_new

double NeutFVM::getLTE(double dt)
{
  VecWAXPY(workvec,-1.0/(1-TRgam),phigam,phivec);
  VecScale(workvec,1.0/TRgam);
  VecAXPY(workvec,1.0/(1-TRgam),phinew);
  double wnorm;
  VecNorm(workvec,NORM_2,&wnorm);
  /*  now prepare divider  */
  double rnorm;
#ifdef PREC
  MatMult(Mmat,phinew,rhsvec);
  VecNorm(rhsvec,NORM_2,&rnorm);
#else
  VecNorm(phinew,NORM_2,&rnorm);
#endif
  double LTE = wnorm/(epsR*rnorm+epsA);
  if(prlev > 1) {
    int ppos;
    double Tpeak;
    VecMax(phinew,&ppos,&Tpeak);
    printf("LTE=%.3le %.3le %.2le  %.2le\n",wnorm,rnorm,LTE, Tpeak);
  }
  return LTE;
};      // NeutFVM::getLTE

double NeutFVM::getTotPow()
{
  if(do_power) {
    ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    do_power = 0;
  }
  double Totpow;
  VecSum(powvec,&Totpow);
  return Totpow;
};

double NeutFVM::getPeakPowDen()
{
  if(do_power) {
    ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    do_power = 0;
  }
  ierr = VecPointwiseDivide(powden,powvec,volvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  PetscInt mpos;
  PetscScalar xmax;
  VecMax(powden,&mpos,&xmax);
  return xmax;
};	// NeutFVM::getPeakPowDen

Vec NeutFVM::getPow()
{
  if(do_power) {
    ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    do_power = 0;
  }
  return powvec;
};	// NeutFVM::getPow

Vec NeutFVM::getPowden()
{
  if(do_power) {
    ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    do_power = 0;
  }
  ierr = VecPointwiseDivide(powden,powvec,volvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("powden:");
    VecView(powden,PETSC_VIEWER_STDOUT_WORLD);
  }
  return powden;
};

Vec NeutFVM::getFiss()
{
  if(do_power) {
    ierr = MatMultTranspose(Pmat,phivec,powvec);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    do_power = 0;
  }
  ierr = VecCopy(powvec,fissvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  return fissvec;
};      // NeutFVM::getPow


void NeutFVM::setDsrc(Vec dsrcvec)
{
  /*  set src (ndim)  using dsrcvec (npows) */
  if(prlev > 1) {
    printf("dsrcvec:\n");
    VecView(dsrcvec,PETSC_VIEWER_STDOUT_WORLD);

    printf("oldsrc:\n");
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  ierr = MatMult(Smat,dsrcvec,srcvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("srcvec:\n");
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  }
};	// NeutFVM::setDsrc
