/*
 * 	Cond3D.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Cond3D.h"


Cond3D::Cond3D()
{
  prlev = 0;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  mdotcp = 100.0;
  Tin = 0.0;
  do_setup = 1;
  do_updTK = 1;
  do_updTB = 1;
  kndtri = NULL;
  kndchn = NULL;
  Amat = NULL;
  P2Smat = NULL;
  TPvec  = NULL;
  powvol = NULL;
  powvec = NULL;
  pxvec = NULL;
  hcoefvec = NULL;
  Tbulkvec = NULL;
  Toutmat = NULL;
  Toutvec = NULL;
  TKmat = NULL;
  TKvec = NULL;
  TBmat = NULL;
  TBvec = NULL;
  solver = new Solver();
  thpro = NULL;
  TK = NULL;
  flu = NULL;
  TB = NULL;
  PTwork = NULL;
  /*  transient */
  TRgam = 2.0-sqrt(2.0);
  TRC2 =  2*fabs(((-3*TRgam+4)*TRgam-2)/(12*(2-TRgam)));
  epsR = 1.0e-5;
  epsA = 1.0e-5;
  Mmat = NULL;
  TRmat = NULL;
};	// Cond3D::Cond3D

Cond3D::~Cond3D()
{
  delete [] kndtri;
  delete [] kndchn;
  delete solver;
  delete [] TK;
  delete [] flu;
  delete [] TB;
  delete [] PTwork;
  MatDestroy(&Amat);
  MatDestroy(&P2Smat);
  VecDestroy(&TPvec);
  VecDestroy(&powvol);
  VecDestroy(&powvec);
  VecDestroy(&pxvec);
  VecDestroy(&hcoefvec);
  VecDestroy(&Tbulkvec);
  MatDestroy(&Toutmat);
  VecDestroy(&Toutvec);
  MatDestroy(&TKmat);
  VecDestroy(&TKvec);
  MatDestroy(&TBmat);
  VecDestroy(&TBvec);

  MatDestroy(&Mmat);
  MatDestroy(&TRmat);
};	// Cond3D::~Cond3D

void Cond3D::setPrlev(int prl)
{
  prlev = prl;
};	// Cond3D::setPrlev

void Cond3D::setGeom(Geom *geom_i)
{
  geom = geom_i;
  ntri = geom->getNtri();
  nz   = geom->getNz();
  nvols = ntri*nz;
  nchan = geom->getNchan();
  npows = geom->getNpows();
};      // Core3D::setGeom

void Cond3D::setFluence(double flu_i[])
{
  assert(flu_i != NULL);
  assert(ntri*nz > 0);
  if(flu == NULL) flu = new double[ntri*nz];
  for(int tz=0; tz < ntri*nz; tz++) flu[tz] = flu_i[tz];
};	// Cond3D::setFluence

void Cond3D::setTHpro(THpro *thpro_i)
{
  thpro = thpro_i;
};	// Cond3D::setTHpro

void Cond3D::setup()
{
  assert(Amat == NULL);
  prepare();
};	// Cond3D::setup

void Cond3D::prepare()
{
  makeChannelMap();
  ndim = (ntri+nchan)*nz;
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

  ierr = VecCreate(PETSC_COMM_WORLD,&srcvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(srcvec,PETSC_DECIDE,ndim);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(srcvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(srcvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(srcvec,&BCvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  VecDuplicate(srcvec,&Tvec);

  /*  prepare src conversion matrix  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,npows, maxdiag,NULL, maxoff,NULL, &P2Smat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(P2Smat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(P2Smat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&TPvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(TPvec,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(TPvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(TPvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecDuplicate(TPvec,&powvol);
  ierr = VecDuplicate(TPvec,&powvec);
  ierr = VecDuplicate(TPvec,&pxvec);
  ierr = VecDuplicate(TPvec,&hcoefvec);
  ierr = VecDuplicate(TPvec,&Tbulkvec);

  /*  prepare Tout extraction matrix */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,nchan, maxdiag,NULL, maxoff,NULL, &Toutmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Toutmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Toutmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&Toutvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(Toutvec,PETSC_DECIDE,nchan);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(Toutvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(Toutvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  /*	prepare TK extraction matrix  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ntri*nz, maxdiag,NULL, maxoff,NULL, &TKmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(TKmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(TKmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&TKvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(TKvec,PETSC_DECIDE,ntri*nz);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(TKvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(TKvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  /*    prepare TB extraction matrix  */
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,nchan*nz, maxdiag,NULL, maxoff,NULL, &TBmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(TBmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(TBmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&TBvec);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(TBvec,PETSC_DECIDE,nchan*nz);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(TBvec);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(TBvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  makeToutmat();
  makeP2Smat();		// we will prepare P2S and TK
  makeTBmat();
  makePowvol();

  if(thpro == NULL) {
    if(mpid==0) printf("*** We use built in TH properties.\n");
  }
  /*  temperatures for conductivity calculation  */
  TK = new double[ntri*nz];
  for(int tz=0; tz < ntri*nz; tz++) {
    TK[tz] = Tin;
  }
  if(flu == NULL) {
    flu = new double[ntri*nz];
    for(int tz=0; tz < ntri*nz; tz++) {
      flu[tz] = 1.0;	// default fluence
    }
  }

  /*  for heat transfer coefficients  */

  TB = new double[nchan*nz];
  for(int cz=0; cz < nchan*nz; cz++) TB[cz] = Tin;
  /*  work area for expandPTvec  */
  PTwork = new double[ntri*nz];
  for(int tz=0; tz < ntri*nz; tz++) PTwork[tz] = 0.0;

};	// Cond3D::prepare

void Cond3D::makeChannelMap()
{
  /*  create index array for channel */
  kndtri = new int[ntri];
  kndchn = new int[nchan];
  int *chtri = geom->Chtri();
  int k=0;
  for(int t=0; t < ntri; t++) {
    kndtri[t] = k;
    int c = chtri[t];
    if(c >= 0) {
      k++;
      kndchn[c] = k;
    }
    k++;
  }
  if(prlev) {
    printf("kndtri=");
    for(int t=0; t < ntri; t++) printf(" %d",kndtri[t]);
    printf("\n");
    printf("kndchn=");
    for(int c=0; c < nchan; c++) printf(" %d",kndchn[c]);
    printf("\n");
    printf(" all=%d\n",ntri+nchan);
  }
};	//  Cond3D::makeChannelMap


void Cond3D::makeP2Smat()
{
  /*  P2S : matrix to convert power into source  */
  /*  TK  : mapping Tvec into volume quantitis */
  ierr = MatZeroEntries(P2Smat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatZeroEntries(TKmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *chtri = geom->Chtri();
  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();
  double one[1];
  one[0] = 1.0;
  int kf = 0;
  for(int t=0; t < ntri; t++) {
    int fisflag = 0;
    int c = chtri[t];
    if((c >= 0) && ispp[c]) fisflag = 1;
    for(int z=0; z < nz; z++) {
      int ktz = kndtri[t]*nz+z;
      if(fisflag && ispowz[z]) {
        /*  we cannot continue due to counter kf  */
        if((ktz >= rAbeg) && (ktz < rAend)) {
          int indr[1],indc[1];
          indr[0] = kndtri[t]*nz+z;
          indc[0] = kf;
          ierr = MatSetValues(P2Smat,1,indr,1,indc,one, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
        }
        kf++;
      }
    }   // for z
    /* for volume quantities */
    for(int z=0; z < nz; z++) {
      int ktz = kndtri[t]*nz+z;
      if((ktz < rAbeg) || (ktz>= rAend)) continue;
      int indr[1],indc[1];
      indr[0] = ktz;
      indc[0] = t*nz+z;
      ierr = MatSetValues(TKmat,1,indr,1,indc,one, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
    }	// for z
  }     // for t
  assert(kf == npows);
  ierr = MatAssemblyBegin(P2Smat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyBegin(TKmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(P2Smat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(TKmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // Cond3D::makeP2Smat

void Cond3D::makeTBmat()
{
  /*  make TBmat together */
  ierr = MatZeroEntries(TBmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  double one[1];
  one[0] = 1.0;
  for(int c=0; c < nchan; c++) {
    for(int z=0; z < nz; z++) {
      int indr[1],indc[1];
      indr[0] = kndchn[c]*nz+z;        // index of channel exit
      if((indr[0] < rAbeg) || (indr[0] >= rAend)) continue;
      indc[0] = c*nz+z;
      ierr = MatSetValues(TBmat,1,indr,1,indc,one, INSERT_VALUES);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
    }	// for z
  }	// for c
  ierr = MatAssemblyBegin(TBmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(TBmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};      // Cond3D::makeTBmat

void Cond3D::makeToutmat()
{
  /*  vector to obtain channel exit temperatures */
  ierr = MatZeroEntries(Toutmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  double one[1];
  one[0] = 1.0;
  for(int c=0; c < nchan; c++) {
    int indr[1],indc[1];
    indr[0] = kndchn[c]*nz+nz-1;	// index of channel exit
    if((indr[0] < rAbeg) || (indr[0] >= rAend)) continue;
    indc[0] = c;
    ierr = MatSetValues(Toutmat,1,indr,1,indc,one, INSERT_VALUES);
	 CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  ierr = MatAssemblyBegin(Toutmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Toutmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// Cond3D::makeToutmat

void Cond3D::makePowvol()
{
  /*  make a vector for volume of power producing element */
  /*	to be used to convert input power density to volume quantity  */
  ierr = VecSet(powvol,0.0);  CHKERRABORT(PETSC_COMM_SELF,ierr);

  int rPbeg,rPend;
  ierr = VecGetOwnershipRange(powvol,&rPbeg,&rPend);
         CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *chtri = geom->Chtri();
  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();
  double vol[1];
  int kf = 0;
  for(int t=0; t < ntri; t++) {
    int fisflag = 0;
    int c = chtri[t];
    if((c >= 0) && ispp[c]) fisflag = 1;
    for(int z=0; z < nz; z++) {
      if(fisflag && ispowz[z]) {
        /*  we cannot continue due to counter kf  */
        if((kf >= rPbeg) && (kf < rPend)) {
          int indr[1];
          indr[0] = kf;
          vol[0] = geom->getVol(t,z);
          ierr = VecSetValues(powvol,1,indr,vol, ADD_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
        }
        kf++;
      }
    }   // for z
  }	// for t
  ierr = VecAssemblyBegin(powvol);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(powvol);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("* %s:%d powvol:\n",__FILE__,__LINE__);
    VecView(powvol,PETSC_VIEWER_STDOUT_WORLD);
  }
};

void Cond3D::setMdot(double Mdot)
{
  /*  Mdot(kg/s)  */
  /*  mdot * cp density (W/m^2/K)  */
  double charea = geom->getChanArea();
  double heCp = 5195;	// (J/kg/K)
  mdotcp = Mdot*heCp/charea;	
  if(mpid==0) printf("== mdotcp=%.2le (W/m^2/K)\n",mdotcp);
  do_setup = 1;	// require to setup system matrix
};	// Cond3D::setMdot

void Cond3D::setTin(double Tin_i)
{
  Tin = Tin_i;
};	// Cond3D::setTin

void Cond3D::setupAmat()
{
  ierr = MatZeroEntries(Amat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSet(BCvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  double *dz = geom->Dz();
  int *chtri = geom->Chtri();
  int *mat3d = geom->Mat3D();
  double *BC = new double[nz];
  int *indbc = new int[nz];
  double hAcoefz = 1.6e+4;	// (W/m^3/K)
  hAcoefz = 100.0;
  if(prlev) printf("ntri=%d nz=%d\n",ntri,nz);
  /* int noside = 1;  // to check reflective side  */
  double xheatz;
  for(int t=0; t < ntri; t++) {
    int isbc = 0;
    for(int z=0; z < nz; z++) {
      BC[z] = 0;
      indbc[z] = -1;
    }
    int c = chtri[t];
    for(int z=0; z < nz; z++) {
      int tz = t*nz+z;	// for material property
      int ktz = kndtri[t]*nz+z;	// for matrix
      if((ktz < rAbeg) || (ktz >= rAend)) continue;
      int indr[1];
      indr[0] = ktz;
      /* 5 faces for a tri-prism  */
      int tb[5];
      double area[5],dist[5];
      double vol = geom->getFVM(t,z,tb,area,dist);
      if(prlev > 3) {
        printf("t=%d z=%d  A=%.3lf:",t,z,vol);
        for(int e=0; e < 5; e++) printf(" %d(%.1le,%.1le)",
          tb[e],area[e],dist[e]);
        printf("\n");
      }
      double condK = 5.0;
      if(thpro) condK = thpro->cond(mat3d[tz],TK[tz],flu[tz]);
      double tauKL,condL;
      int indc[6];
      double A[6];
      indc[0] = ktz;
      if(c >= 0) {
        if(thpro) hAcoefz = thpro->hAcoef(TB[c*nz+z]);
        A[0] = hAcoefz*vol;
      }
      else A[0] = 0;

      for(int i=1; i < 6; i++) {
        indc[i] = -1;
        A[i] = 0.0;
      }
      for(int e=0; e < 5; e++) {
        int tzL = tb[e];
        /* if(noside && (e > 1) && (tzL == -2)) tzL = -1;  */
        /* enforce Neumann for bottom  */
        if((e==1) && (z == nz-1)) tzL = -1;
        if(tzL >= 0) {
          double distL = geom->getDistL(tzL,t*nz+z);
          condL = condK;	// assume same conductivity
          if(thpro) condL = thpro->cond(mat3d[tzL],TK[tzL],flu[tzL]);
          tauKL = area[e]/(dist[e]/condK+distL/condL);
        }
        else if(tzL == -2) {    // Dirichlet BC
          tauKL = area[e]/(dist[e]/condK);
          isbc = 1;
          BC[z] += tauKL;	// check
          indbc[z] = ktz;
          if(prlev > 1) printf("t=%d z=%d BC=%.2le\n",t,z,BC[z]);
        }
        else tauKL = 0.0;       // reflective BC
        /* translate back */
        int et = tzL/nz;
        int ez = tzL%nz;
        indc[1+e] = kndtri[et]*nz+ez;
        A[1+e] = -tauKL;
        A[0] += tauKL;
      } // for e
      if(prlev > 3) {
        printf("A %d:z",ktz);
        for(int i=0; i < 6; i++) printf(" %.2le(%d)",A[i],indc[i]);
        printf("\n");
      }
      ierr = MatSetValues(Amat,1,indr,6,indc, A,ADD_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
      if(c >= 0) {
        int indcoef[2];
        double valch[2];
        indcoef[0] = -1;
        if(z==0) {
          valch[0] = hAcoefz/2*vol;
          ierr = VecSetValues(BCvec,1,indr,valch,ADD_VALUES);
        }
        else {
          indcoef[0] = kndchn[c]*nz+z-1;
          valch[0] = -hAcoefz/2*vol;
        }
        indcoef[1] = kndchn[c]*nz+z;
        valch[1] =  -hAcoefz/2*vol;
        ierr = MatSetValues(Amat,1,indr,2,indcoef, valch,ADD_VALUES);
            CHKERRABORT(PETSC_COMM_SELF,ierr);
      }
    }   // for z
    if(isbc) {
      /* modify for MPI */
      for(int z=0; z < nz; z++) {
        if((indbc[z] < rAbeg) || (indbc[z] >= rAend)) indbc[z] = -1;
      }
      ierr = VecSetValues(BCvec,nz,indbc, BC,ADD_VALUES);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
    if(c >= 0) {
      /*	for channel  */
      for(int z=0; z < nz; z++) {
        int indr[1],indc[3];
        double valch[3];
        indr[0] = kndchn[c]*nz+z;
        if((indr[0] < rAbeg) || (indr[0] >= rAend)) continue;
        if(thpro) hAcoefz = thpro->hAcoef(TB[c*nz+z]);
        xheatz = hAcoefz*dz[z]/mdotcp;
        /*  upside */
        if(z == 0) {
          indc[0] = -1;
          valch[0] = (1-xheatz/2);
          ierr = VecSetValues(BCvec,1,indr,valch,ADD_VALUES);
	    CHKERRABORT(PETSC_COMM_SELF,ierr);
        }
        else {
          indc[0] = kndchn[c]*nz+z-1;
          valch[0] = -1+xheatz/2;
        }
        /*  down side */
        indc[1] = kndchn[c]*nz+z;
        valch[1] = 1+xheatz/2;
        /*  volume element side */
        indc[2] = kndtri[t]*nz+z;
        valch[2] = -xheatz;
        ierr = MatSetValues(Amat,1,indr,3,indc,valch,ADD_VALUES);
	  CHKERRABORT(PETSC_COMM_SELF,ierr);
      }	// for z
    }	// for c
  }     // for t
  ierr = MatAssemblyBegin(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyBegin(BCvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Amat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(BCvec);
  if(prlev > 1) {
    printf("%s:%d Amat:\n",__FILE__,__LINE__);
    MatView(Amat,PETSC_VIEWER_STDOUT_WORLD);
    printf("BCvec:\n");
    VecView(BCvec,PETSC_VIEWER_STDOUT_WORLD);
    if(prlev > 5) exit(0);
  }
  delete [] BC;
  delete [] indbc;

  do_setup = 0;
};	// Cond3D::setupAmat

void Cond3D::steady(Vec powden)
{
  ierr = VecPointwiseMult(powvec,powden,powvol);
  if(do_setup) setupAmat();
  /*  convert powvec to srcvec  */
  ierr = MatMult(P2Smat,powvec,srcvec);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAXPY(srcvec,Tin,BCvec);
  if(prlev > 1) {
    printf("%s:%d srcvec:\n",__FILE__,__LINE__);
    VecView(srcvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  solver->solve(Amat,srcvec,Tvec);
  /* VecView(Tvec,PETSC_VIEWER_STDOUT_WORLD);  */
  do_updTK = 1;
  do_updTB = 1;	// temperature array need update
};	// Cond3D::steady

void Cond3D::start()
{
  /*  prepare for transient  */
  int maxdiag = ndim;
  int maxoff = ndim;
  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &Mmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(Mmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(Mmat); CHKERRABORT(PETSC_COMM_SELF,ierr);

  setupMmat();

  ierr =  MatCreateAIJ(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,
        ndim,ndim, maxdiag,NULL, maxoff,NULL, &TRmat);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatSetOption(TRmat,MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
  ierr = MatSetUp(TRmat); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyBegin(TRmat,MAT_FINAL_ASSEMBLY);
	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(TRmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);

  /*  vectors */
  VecDuplicate(Tvec,&Tgam);
  VecDuplicate(Tvec,&Tnew);
  VecDuplicate(Tvec,&srcnew);
  VecDuplicate(Tvec,&workvec);
  VecDuplicate(Tvec,&rhsvec);
};	// Cond3D::start

void Cond3D::setupMmat()
{
  /*  create matrix for transient */
  ierr = MatZeroEntries(Mmat);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  double rhocp = 3000.0;
  for(int t=0; t < ntri; t++) {
    for(int z=0; z < nz; z++) {
      int ktz = kndtri[t]*nz+z;
      if((ktz < rAbeg) || (ktz >= rAend)) continue;
      int indr[1],indc[1];
      indr[0] = ktz;
      indc[0] = ktz;
      double M[1];
      M[0] = rhocp;
      ierr = MatSetValues(Mmat,1,indr,1,indc, M,ADD_VALUES);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  }
  ierr = MatAssemblyBegin(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = MatAssemblyEnd(Mmat,MAT_FINAL_ASSEMBLY);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
};	// Cond3D::setupMmat

void Cond3D::setLTE(double epsR_i, double epsA_i)
{
  /*  set LTE error estimator */
  epsR = epsR_i;
  epsA = epsA_i;
};	// Cond3D::setLTE

double Cond3D::step(double dt, Vec powden)
{
  ierr = VecPointwiseMult(powvec,powden,powvol);
  if(do_setup) setupAmat();
  /*  prepare new source  */
  ierr = MatMult(P2Smat,powvec,srcnew);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAXPY(srcnew,Tin,BCvec);

  step_gam(dt);
  step_new(dt);
  double LTE = getLTE(dt);
  return LTE;
};	// Cond3D::step

void Cond3D::march()
{
  ierr = VecCopy(Tnew,Tvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecCopy(srcnew,srcvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  do_updTK = 1;
  do_updTB = 1; // temperature array need update
};	//  Cond3D::march

void Cond3D::step_gam(double dt)
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
};	// Cond3D::step_gam

void Cond3D::step_new(double dt)
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
  if(prlev > 2) {
    ierr = VecCopy(Tnew,Tvec);	CHKERRABORT(PETSC_COMM_SELF,ierr);
    printAx(3);
    printAx(4);
  }
};	// Cond3D::step_new

double Cond3D::getLTE(double dt)
{
  VecWAXPY(workvec,-1.0/(1-TRgam),Tgam,Tvec);
  VecScale(workvec,1.0/TRgam);
  VecAXPY(workvec,1.0/(1-TRgam),Tnew);
  double wnorm;
  VecNorm(workvec,NORM_2,&wnorm);
  /*  now prepare divider  */
  double rnorm;
#ifdef PREC
  /*  strict definition of divider  */
  MatMult(Mmat,Tnew,rhsvec);
  VecNorm(rhsvec,NORM_2,&rnorm);
#else
  VecNorm(Tnew,NORM_2,&rnorm);
#endif
  double LTE = wnorm/(epsR*rnorm+epsA);
  return LTE;
};	// Cond3D::getLTE

double Cond3D::getTpeak()
{
  /*  get peak temperature of problem  */
  int posmax;
  double Tmax;
  VecMax(Tvec,&posmax,&Tmax);
  return Tmax;
};	// Cond3D::getTpeak

double Cond3D::getTPavg()
{
  /*  get average of all power producing nodes  */
  ierr = MatMultTranspose(P2Smat,Tvec,TPvec);
  double Tpsum;
  VecSum(TPvec,&Tpsum);
  double TPavg = Tpsum/npows;	// assume all power prod. volume is same szie
  return TPavg;
};	// Cond3D::getTPavg

double Cond3D::getTexit()
{
  ierr = MatMultTranspose(Toutmat,Tvec,Toutvec);
	
  if(prlev > 1) VecView(Toutvec,PETSC_VIEWER_STDOUT_WORLD);
  double Texit;
  VecSum(Toutvec,&Texit);
  return Texit/nchan;
};	// Cond3D::getTexit

void Cond3D::printAx(int t)
{
  Vec Twork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(Tvec,&ctx,&Twork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,Tvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,Tvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);
  if(mpid != 0) return;

  double *T;
  ierr = VecGetArray(Twork,&T);
  printf("T%d:",t);
  for(int z=0; z < nz; z++) printf(" %.2lf",T[t*nz+z]);
  printf("\n");
  ierr = VecRestoreArray(Twork,&T);
};	// Cond3D::printAx

double Cond3D::updateTemp()
{
  /*  update TK and TB array  */
  double errTK = updateTK();
  double errTB = updateTB();
  double err = errTK;
  if(err < errTB) err = errTB;
  do_setup = 1;
  return err;
};	// Cond3D::updateTemp

double Cond3D::updateTK()
{
  /*  update TK array  */
  ierr = MatMultTranspose(TKmat,Tvec,TKvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("TKvec");
    VecView(TKvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  /*  scatter and collect to all processor  */
  Vec Twork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(TKvec,&ctx,&Twork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,TKvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,TKvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);
  /*  we make array for all processor  */
  double *T;
  ierr = VecGetArray(Twork,&T);
  double errmax = 0;
  for(int tz=0; tz < ntri*nz; tz++) {
    double err = fabs(TK[tz]-T[tz]);
    if(errmax < err) errmax = err;
    TK[tz] = T[tz];
  }
  if(prlev > 1) printf("*** %s:%d maxerr=%.2le\n",__FILE__,__LINE__,errmax);
  ierr = VecRestoreArray(Twork,&T);
  do_setup = 1;
  do_updTK = 0;
  return errmax;
};	// Cond3D::updateTK

double* Cond3D::getTK()
{
  if(do_updTK) updateTK();
  return TK;
};	//  Cond3D::getTK

double Cond3D::updateTB()
{
  /*  update TB array  */
  ierr = MatMultTranspose(TBmat,Tvec,TBvec);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  if(prlev > 1) {
    printf("TBvec");
    VecView(TBvec,PETSC_VIEWER_STDOUT_WORLD);
  }

  /*  scatter and collect to all processor  */
  Vec Twork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(TBvec,&ctx,&Twork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,TBvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,TBvec,Twork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);
  /*  we make array for all processor  */
  double *T;
  ierr = VecGetArray(Twork,&T);
  double errmax = 0;
  for(int c=0; c < nchan; c++) {
    int cz = c*nz;
    double Tb = (T[cz]+Tin)/2;
    double err = fabs(TB[cz]-Tb);
    if(errmax < err) errmax = err;
    TB[cz] = Tb;
    for(int z=1; z < nz; z++) {
      cz = c*nz+z;
      Tb = (T[cz]+T[cz-1])/2;
      err = fabs(TB[cz]-Tb);
      if(errmax < err) errmax = err;
      TB[cz] = Tb;
    }
  }	// for c
  if(prlev > 1) printf("*** %s:%d maxerr=%.2le\n",__FILE__,__LINE__,errmax);
  ierr = VecRestoreArray(Twork,&T);
  do_setup = 1;
  do_updTB = 0;
  return errmax;
};      // Cond3D::updateTB

Vec Cond3D::getPhcoef()
{
  if(do_updTB) updateTB();
  /*  get hcoef for power producing elements */
  if(thpro==NULL) printf("*** cannot call getPhcoef ***\n");
  ierr = VecSet(hcoefvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  int rPbeg,rPend;
  ierr = VecGetOwnershipRange(hcoefvec,&rPbeg,&rPend);
         CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();

  int kf = 0;
  int indr[1];
  double val[1];
  for(int c=0; c < nchan; c++) {
    if(ispp[c]) {
      for(int z=0; z < nz; z++) {
        if(ispowz[z]) {
          if((kf >= rPbeg) && (kf < rPend)) {
            double hcoef = 1700.0;
            if(thpro) hcoef = thpro->hcoef(TB[c*nz+z]);
            indr[0] = kf;
            val[0] = hcoef;
            ierr = VecSetValues(hcoefvec,1,indr,val,INSERT_VALUES);
		CHKERRABORT(PETSC_COMM_SELF,ierr);
          }
          kf++;
        }
      }	// for z
    }	// if ispp
  }	// for c
  assert(kf == npows);
  ierr = VecAssemblyBegin(hcoefvec);
  ierr = VecAssemblyEnd(hcoefvec);
  return hcoefvec;
};	// Cond3D::getPhcoef

Vec Cond3D::getPTbulk()
{
  if(do_updTB) updateTB();
  /*  get Tbulk for power producting elements */
  ierr = VecSet(Tbulkvec,0.0);  CHKERRABORT(PETSC_COMM_SELF,ierr);
  int rPbeg,rPend;
  ierr = VecGetOwnershipRange(Tbulkvec,&rPbeg,&rPend);
         CHKERRABORT(PETSC_COMM_SELF,ierr);

  int *ispp  = geom->Ispp();
  int *ispowz = geom->Ispowz();

  int kf = 0;
  int indr[1];
  double val[1];
  for(int c=0; c < nchan; c++) {
    if(ispp[c]) {
      for(int z=0; z < nz; z++) {
        if(ispowz[z]) {
          if((kf >= rPbeg) && (kf < rPend)) {
            indr[0] = kf;
            val[0] = TB[c*nz+z];
            ierr = VecSetValues(Tbulkvec,1,indr,val,INSERT_VALUES);
                CHKERRABORT(PETSC_COMM_SELF,ierr);
          }
          kf++;
        }
      } // for z
    }   // if ispp
  }     // for c
  assert(kf == npows);
  ierr = VecAssemblyBegin(Tbulkvec);
  ierr = VecAssemblyEnd(Tbulkvec);
  return Tbulkvec;
};	// Cond3D::getPTbulk

double* Cond3D::expandPvec(Vec Pvec)
{
  /* expand power zone vector to core wide array  */
  Vec Pwork;
  VecScatter ctx;
  ierr = VecScatterCreateToAll(Pvec,&ctx,&Pwork);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecScatterBegin(ctx,Pvec,Pwork,INSERT_VALUES,SCATTER_FORWARD);
  ierr = VecScatterEnd(ctx,Pvec,Pwork,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterDestroy(&ctx);

  double *pT;
  ierr = VecGetArray(Pwork,&pT);
  /* collect pT which is defined on power producing nodes  */
  int *tchan = geom->Tchan();
  int *ispowz = geom->Ispowz();
  int kf = 0;
  for(int c=0; c < nchan; c++) {
    int t = tchan[c];
    for(int z=0; z < nz; z++) {
      if(ispowz[z]) {
        int tz = t*nz+z;
        PTwork[tz] = pT[kf];
        kf++;
      }
    }   // for z
  }     // for t
  assert(kf == npows);
  ierr = VecRestoreArray(Pwork,&pT);

  return PTwork;
};	// Cond3D::expandPvec

Vec Cond3D::contractPdouble(double x[])
{
  /*  contract core wide array into vector in power region  */
  ierr = VecSet(pxvec,0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  int jrange,rPbeg,rPend;
  ierr = VecGetOwnershipRange(pxvec,&rPbeg,&rPend);
         CHKERRABORT(PETSC_COMM_SELF,ierr);
  jrange = rPend-rPbeg;
  int *indr = new int[jrange];
  double *val = new double[jrange];
  int *tchan = geom->Tchan();
  int *ispowz = geom->Ispowz();
  int kf = 0;
   for(int c=0; c < nchan; c++) {
    int t = tchan[c];
    int isin = 0;
    int ii = 0;
    for(int z=0; z < nz; z++) {
      if(ispowz[z]) {
        int tz = t*nz+z;
        indr[ii] = kf;
        if((kf < rPbeg) || (kf >= rPend)) indr[ii] = -1;
        else isin++;
        val[ii] = x[tz];
        ii++;
        kf++;
      }
    }   // for z
    if(isin >0 ) {
      ierr = VecSetValues(pxvec,ii,indr,val,INSERT_VALUES);
	 CHKERRABORT(PETSC_COMM_SELF,ierr);
    }
  }	// for c
  assert(kf == npows);
  ierr = VecAssemblyBegin(pxvec); CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecAssemblyEnd(pxvec); CHKERRABORT(PETSC_COMM_SELF,ierr);
  return pxvec;
};	// Cond3D::contractPdouble
