/*
 * 	MPwork.cpp
 * 	material properties for SegCond and TriCond
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <lapacke.h>
#include "MPwork.h"

MPwork::MPwork()
{
  prlev = 0;
  jmon = -1;
  Psys = 3.0e+6;
  emiss = 0.8;	// graphite emissivity
  nzone = 0;
  zrad = NULL;
  condL = NULL;
  rhocpL = NULL;
  Mat_keff = NULL;
  b_keff = NULL;
  ipiv_keff = NULL;
  TK_triso = 0.0;	// flag for Triso property calculation
};

MPwork::~MPwork()
{
  delete [] zrad;
  delete [] condL;
  delete [] rhocpL;
  delete [] Mat_keff;
  delete [] b_keff;
  delete [] ipiv_keff;
};	// MPwork::~MPwork

void MPwork::setPrlev(int prl)
{
  prlev = prl;
};	// MPwork::setPrlev

void MPwork::setMonitor(int jmon_i)
{
  /*  set monitor node  */
  jmon = jmon_i;
};	// MPwork::setMonitor

void MPwork::set(double Psys_i)
{
  Psys = Psys_i;
};

void MPwork::setEmiss(double emiss_i)
{
  emiss = emiss_i;
};	// MPwork::setEmiss

void MPwork::setTriso(double dKer, int nlays, double tlay[], double PF_i)
{
  PF = PF_i;
  /*  save zone radius for effCond calculation  */
  nzone = nlays+2;
  zrad = new double[nzone];
  zrad[0] = dKer/2;
  for(int k=0; k < nlays; k++) {
    zrad[k+1] = zrad[k] + tlay[k];
  }
  double rP = zrad[nlays];
  double rM = rP/pow(PF,1.0/3.0);
  zrad[nlays+1] = rM;
  condL = new double[nzone];
  Mat_keff = new double[4*nzone*nzone];
  b_keff = new double[2*nzone];
  ipiv_keff = new int[2*nzone];
  rhocpL = new double[nzone];
};	//  MPwork::setTriso

void MPwork::setCondL(double flu, double fima, double TK)
{
  /*  set conductivities and heat capacities  */
  if((fabs(flu_triso - flu) < 1.0e-4) && (fabs(fima_triso - fima) < 1.0e-4) &&
	(fabs(TK_triso - TK) < 1.0e-4)) return;
  assert(nzone == 6);
  /*  kernel/buffer/iPyC/SiC/oPyC/NITE  */
  double porUCO = 0.0;
  condL[0] = cond_UCO_Ben(TK,fima,porUCO);
  // condL[0] = cond_UCO(TK);
  double rho_buf = rho_Buf();
  condL[1] = cond_Buf(rho_buf);
  double porPyC = 0.0;
  condL[2] = cond_PyC(porPyC);
  double porSiC = 0.0;
  double TI = TK;
  condL[3] = cond_SiC_Ben(TK,TI,flu);
  condL[4] = condL[2];
  condL[5] = cond_NITE_Ben(TK,TI,flu);
  double Umass = 238.05*0.8 + 235.04*0.2;
  double Cfrac = 0.4;
  double Ofrac = 1.5;
  rhocpL[0] = dens_UO2(TK)*cp_UCO(TK,Umass,Cfrac,Ofrac);
  rhocpL[1] = rho_Buf()*cp_Buf();
  rhocpL[2] = rho_PyC()*cp_PyC(TK);
  rhocpL[3] = rho_SiC()*cp_SiC(TK);
  rhocpL[4] = rhocpL[3];
  rhocpL[5] = rho_SiC()*cp_SiC(TK);	// use SiC data
  flu_triso = flu;
  fima_triso = fima;
  TK_triso = TK;
};

double MPwork::effCond()
{
  assert(condL[0] > 1.0e-5);
  if(prlev > 2) {
    printf("zrad:");
    for(int j=0; j < nzone; j++) printf(" %.1lf",zrad[j]);
    printf("\n");
    printf("condL:");
    for(int j=0; j < nzone; j++) printf(" %.1lf",condL[j]);
    printf("\n");
  }
  /*  TRISO effective conductivity  */
  /* ref) Stainsby, Investigation of Local Heat Transfer Phenomena      */
  /*    in a Pebble Bed HTGR Core, NR001/RP/002 R01, 2009. App.C.       */
  /*  setup matrix  */
  int ndimx = 2*nzone;
  for(int ij=0; ij < ndimx*ndimx; ij++) Mat_keff[ij] = 0;
  for(int i=0; i < nzone; i++) {
    Mat_keff[i*ndimx+i] = 1;
    if(i > 0) {
      Mat_keff[(i-1)*ndimx+i] = -1;
    }
  }
  double a53 = zrad[nzone-2]*zrad[nzone-2]*zrad[nzone-2];
  for(int i=0; i < nzone-1; i++) {
    double ai3 = zrad[i]*zrad[i]*zrad[i];
    Mat_keff[i*ndimx+nzone+i] = a53/ai3;
    Mat_keff[i*ndimx+nzone+i+1] = - a53/ai3;
  }

  for(int i=1; i < nzone; i++) {
    Mat_keff[(nzone+i)*ndimx+i] = condL[i];
    Mat_keff[(nzone+i)*ndimx+i-1] = -condL[i-1];
  }
  Mat_keff[nzone*ndimx+nzone] = 1;
  for(int i=1; i < nzone; i++) {
    double ai3 = zrad[i-1]*zrad[i-1]*zrad[i-1];
    Mat_keff[(nzone+i)*ndimx+nzone+i] = -2*condL[i]*a53/ai3;
    Mat_keff[(nzone+i)*ndimx+nzone+i-1] = 2*condL[i-1]*a53/ai3;
  }
  for(int i=0; i < ndimx; i++) b_keff[i] = 0;
  b_keff[nzone-1] = 1;
  if(prlev > 2) {
    for(int j=0; j < ndimx; j++) {
      printf("M%d:",j);
      for(int i=0; i < ndimx; i++) printf(" %.1le",Mat_keff[j*ndimx+i]);
      printf(" = %.1le\n",b_keff[j]);
    }
  }
  lapack_int info;
  info = LAPACKE_dgesv(LAPACK_ROW_MAJOR,ndimx,1,Mat_keff,ndimx,ipiv_keff,
        b_keff,1);
  if(prlev > 1) {
    printf("bk:");
    for(int j=0; j < ndimx; j++) printf(" %.2le",b_keff[j]);
    printf("\n");
    printf("condL:");
    for(int j=0; j < nzone; j++) printf(" %.3lf",condL[j]);
    printf("\n");
  }
  double km = condL[nzone-1];
  if(prlev > 1) printf("  kmatrix=%.3lf\n",km);
  double Bl = b_keff[ndimx-1];
  if(prlev > 1) printf("km=%.3le Bl=%.3le PF=%.4lf\n",km,Bl,PF);
  double kp = km*(1-2*Bl)/(1+Bl);
  double kcompact = km*(1-2*PF*Bl)/(1+PF*Bl);
  if(prlev > 1) printf("k (particle)=%.6lf (compact=%.6lf) (W/K/m)\n",kp,kcompact);
  return kcompact;
};      //  MPwork::effCond

double MPwork::effRhoCp()
{
  /* find volume average of heat capacity */
  double volsum = 0;
  double rhocpsum = 0;
  double vol0 = 0.0;
  for(int j=0; j < nzone; j++) {
    double vol = zrad[j]*zrad[j]*zrad[j];
    double dvol = vol-vol0;
    rhocpsum += rhocpL[j]*dvol;
    volsum += dvol;
    vol0 = vol;
  }
  return rhocpsum/volsum;
};	// MPwork::effRhoCp

void MPwork::makeSegProp(double flu, double fima, int nmats, double TK[],
	double Rad, double cond[], double rhocp[])
{
  /*	prepare conductivities for SegCond  */
  assert(nmats == 3);
  /*  graphite/helium/compact  */
  cond[0]  = cond_IG110(TK[0],flu);
  rhocp[0] = rho_gra(TK[0]);
  cond[1]  = cond_He(TK[1],Psys) + emiss*Rad;
  rhocp[1] = 5195.0*rho_He(TK[1],Psys);
  /*  compute compact conductivity using dimensions */
  setCondL(flu,fima,TK[2]);
  cond[2] = effCond();
  rhocp[2]  = effRhoCp();
};	// MPwork::makeSegProp

void MPwork::makeTriProp(double flu, double fima, double TK,
        int nmats, double cond[], double rhocp[])
{
  /*	prepare conductivities for TriCond  */
  assert(nmats == 6);
  /*  kernel/buffer/ipyc/sic/opyc/nite	*/
  setCondL(flu,fima,TK);
  /*  we have condL (ker/buf/ipyc/sic/opyc/mat)	*/
  for(int j=0; j < nmats; j++) {
    cond[j] = condL[j];
    rhocp[j] = rhocpL[j];
  }
};	// MPwork::makeTriProp

/*  vector version  */
void MPwork::makeSegProp(Vec fluvec, Vec fimavec, int nmats, Vec Tempvec[],
        Vec radvec, Vec condvec[], Vec rhocpvec[])
{
  /*  prepares condvec and rhocpvec to be used for SegCond  */
  int npows;
  VecGetSize(fluvec,&npows);
  int rTbeg,rTend;
  VecGetOwnershipRange(fluvec,&rTbeg,&rTend);
  int jrange = rTend-rTbeg;
  double *flu,*fima;
  PetscErrorCode ierr;
  ierr = VecGetArray(fluvec,&flu);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(fimavec,&fima);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  assert(nmats == 3);
  double *T0,*T1,*T2;
  ierr = VecGetArray(Tempvec[0],&T0);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(Tempvec[1],&T1);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(Tempvec[2],&T2);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  double *rad;
  ierr = VecGetArray(radvec,&rad);       CHKERRABORT(PETSC_COMM_SELF,ierr);
  /*  prepare condvec, rhocpvec  */
  for(int m=0; m < nmats; m++) {
    ierr = VecSet(condvec[m],0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
    ierr = VecSet(rhocpvec[m],0.0);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  /*  get temperatures  */
  int *indr = new int[jrange];
  double *condM = new double[jrange*nmats];
  double *rocpM = new double[jrange*nmats];
  for(int j=0; j < jrange; j++) {
    double TK[3];
    TK[0] = T0[j];  TK[1] = T1[j];  TK[2] = T2[j];
    double cond[3],rhocp[3];
    makeSegProp(flu[j],fima[j], nmats,TK,rad[j],cond,rhocp);
    indr[j] = rTbeg+j;
    for(int m=0; m < nmats; m++) {
      condM[m*jrange+j] = cond[m];
      rocpM[m*jrange+j] = rhocp[m];
    }
    if((prlev && (j==jmon)) || (prlev > 1)) {
      printf("j=%d :",j);
      printf(" %.2lf %.2lf %.2lf",T0[j],T1[j],T2[j]);
      printf(" %.2le %.2le %.2le",cond[0],cond[1],cond[2]);
      printf(" %.2le %.2le %.2le",rhocp[0],rhocp[1],rhocp[2]);
      printf("\n");
    }
  }	// for j
  for(int m=0; m < nmats; m++) {
    ierr = VecSetValues(condvec[m], jrange,indr,condM+m*jrange,INSERT_VALUES);
    ierr = VecSetValues(rhocpvec[m],jrange,indr,rocpM+m*jrange,INSERT_VALUES);
  }	// for m
  for(int m=0; m < nmats; m++) {
    ierr = VecAssemblyBegin(condvec[m]);
    ierr = VecAssemblyBegin(rhocpvec[m]);
  }
  for(int m=0; m < nmats; m++) {
    ierr = VecAssemblyEnd(condvec[m]);
    ierr = VecAssemblyEnd(rhocpvec[m]);
  }
  ierr = VecRestoreArray(fluvec,&flu);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(fimavec,&fima);    CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(Tempvec[0],&T0);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(Tempvec[1],&T1);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(Tempvec[2],&T2);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(radvec,&rad);	CHKERRABORT(PETSC_COMM_SELF,ierr);
  delete [] condM;
  delete [] rocpM;
};	// MPwork::makeSegProp

void MPwork::makeTriProp(Vec fluvec, Vec fimavec, Vec TKvec,
        int nmats, Vec condvec[], Vec rhocpvec[])
{
  /*  prepares condvec and rhocpvec to be used for SegCond  */
  if(prlev > 1) {
    printf("%s:%d TKvec\n",__FILE__,__LINE__);
    VecView(TKvec,PETSC_VIEWER_STDOUT_WORLD);
  }
  int npows;
  VecGetSize(fluvec,&npows);
  int rTbeg,rTend;
  VecGetOwnershipRange(fluvec,&rTbeg,&rTend);
  int jrange = rTend-rTbeg;
  double *flu,*fima,*TK;
  PetscErrorCode ierr;
  ierr = VecGetArray(fluvec,&flu);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(fimavec,&fima);    CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecGetArray(TKvec,&TK);    CHKERRABORT(PETSC_COMM_SELF,ierr);
  assert(nmats == 6);
  /*  prepare condvec, rhocpvec  */
  for(int m=0; m < nmats; m++) {
    ierr = VecSet(condvec[m],0.0);      CHKERRABORT(PETSC_COMM_SELF,ierr);
    ierr = VecSet(rhocpvec[m],0.0);     CHKERRABORT(PETSC_COMM_SELF,ierr);
  }
  /*  get temperatures  */
  int *indr = new int[jrange];
  double *condM = new double[jrange*nmats];
  double *rocpM = new double[jrange*nmats];
  for(int j=0; j < jrange; j++) {
    double cond[6],rhocp[6];
    makeTriProp(flu[j],fima[j],TK[j], nmats,cond,rhocp);
    indr[j] = rTbeg+j;
    for(int m=0; m < nmats; m++) {
      condM[m*jrange+j] = cond[m];
      rocpM[m*jrange+j] = rhocp[m];
    }
    if(prlev > 1) {
      printf("j=%d :",j);
      printf(" %.2lf" ,TK[j]);
      printf(" %.2le %.2le %.2le",cond[0],cond[1],cond[2]);
      printf(" %.2le %.2le %.2le",rhocp[0],rhocp[1],rhocp[2]);
      printf("\n");
    }
  }     // for j
  for(int m=0; m < nmats; m++) {
    ierr = VecSetValues(condvec[m], jrange,indr,condM+m*jrange,INSERT_VALUES);
    ierr = VecSetValues(rhocpvec[m],jrange,indr,rocpM+m*jrange,INSERT_VALUES);
  }     // for m
  for(int m=0; m < nmats; m++) {
    ierr = VecAssemblyBegin(condvec[m]);
    ierr = VecAssemblyBegin(rhocpvec[m]);
  }
  for(int m=0; m < nmats; m++) {
    ierr = VecAssemblyEnd(condvec[m]);
    ierr = VecAssemblyEnd(rhocpvec[m]);
  }
  ierr = VecRestoreArray(fluvec,&flu);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(fimavec,&fima);    CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecRestoreArray(TKvec,&TK);   CHKERRABORT(PETSC_COMM_SELF,ierr);
  delete [] condM;
  delete [] rocpM;
};	// MPwork::makeTriProp
