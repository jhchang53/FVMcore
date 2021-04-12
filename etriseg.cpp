/*
 * 	etriseg.cpp
 * 	test Cond3D and SegCond 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <petscksp.h>
#include "GeomEtri.h"
#include "Cond3D.h"
#include "SegCond.h"


int main(int argc, char **argv)
{
  /* initial protocol to call petsc */
  PetscErrorCode ierr;
  ierr = PetscInitialize(&argc,&argv,(char*)0,NULL);
  if (ierr) {
    printf("*** problem during MPI initialization.\n");
    return ierr;
  }
  auto start = std::chrono::high_resolution_clock::now();

  int mpid;
  MPI_Comm_rank (MPI_COMM_WORLD, &mpid);
  GeomEtri *geom = new GeomEtri();
  if(mpid==0) geom->setPrlev(1);
  geom->setBdry(0);
  double side = 1.0;
  int ntri = geom->generate(side,1);
  int nz = 10;
  double *dz = new double[nz];
  int *ispowz = new int[nz];
  for(int z=0; z < nz; z++) {
    dz[z] = 10.0/nz;
    ispowz[z] = 1;
  }
  geom->extrude(nz,dz,ispowz);
  delete [] dz;
  delete [] ispowz;

  Cond3D *cond = new Cond3D();
  if(mpid==0) cond->setPrlev(1);
  cond->setGeom(geom);
  double Mdot = 1.0e-2;	// kg/s
  cond->setMdot(Mdot);
  double Tin = 600.0;
  cond->setTin(Tin);
  cond->setup();
  double totpow = 1.0e+4;
  int npows = geom->getNpows();
  Vec powden;
  ierr = VecCreate(PETSC_COMM_WORLD,&powden);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(powden,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(powden);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(powden,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  double corevol = geom->getPowVol();
  double powd = totpow/corevol;
  VecSet(powden,powd);
  double Mdotcp = Mdot*5195.0;
  if(mpid==0) printf(" powd=%.3le dT=%.2lf \n",powd,totpow/Mdotcp);
  for(int iter=0; iter < 3; iter++) {
    cond->steady(powden);
    double err = cond->updateTemp();
    if(err < 1.0) break;
  }
  cond->printAx(0);
  cond->printAx(1);

  /*  now call SegCond  */
  Vec hcoefvec = cond->getPhcoef();
  Vec Tbulkvec = cond->getPTbulk();
  
  SegCond *seg = new SegCond();
  seg->setGeom();
  seg->setup(npows);
  seg->setPfactor(1000.0); // 4.5);	// multiplier to convert core power density to compact
  seg->setHcoef(hcoefvec);
  seg->steady(Tbulkvec,powden);
  Vec Tcomp = seg->getTpinVec();
  VecView(Tcomp,PETSC_VIEWER_STDOUT_WORLD);

  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
