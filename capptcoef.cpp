/*
 *	capptcoef.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <petscksp.h>
#include "CoreDesc.h"
#include "GeomCAPP.h"
#include "XScapp.h"
#include "NeutFVM.h"

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

  CoreDesc desc;
  double Totpow = 15.0e+6;      // (Watt)
  double Psys = 3.0e+6;         // (Pa)
  double T0 = 273.15;
  double Tin =  300.0 +T0;          // (degK)
  double Tout = 600.0 +T0;
  desc.setMMR(Totpow,Psys,Tin,Tout);
  // if(mpid == 0) desc.print();
  GeomCAPP *geom = new GeomCAPP();
  int mat[]  = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
  double pitch = 0.322;  // (meter)
  double side = pitch/sqrt(3.0);
  geom->setPlane(side, 6, mat);
  int chan[] = {0, 1, 0,1, 1,1,1, 0,0,0,0, -1,0,0,0,0};
  geom->setChannel(chan);
  /*  set control rod number  */
  int crod[] = {1, 0, 2,0, 0,0,0, 0,0,0,0, -1,0,0,0,0};
  geom->setCrod(crod);

  /* units in meter  */
  double dz[] = {0.350,0.700,0.01736,0.66528,0.03472,0.66528,0.03472,0.66528,
        0.03472,0.66528,0.03472,0.66528,0.03472,0.66528,0.01736,0.70,0.35};
  int zseg[]  = {   2,   4,    1,     4,    1,     4,    1,     4,
            1,     4,    1,     4,    1,     4,    1,   4,   2};
  int ispowz[] = {0,0,0, 1,0,1,0,1, 0,1,0,1,0,1,0,0,0};
  int nz = sizeof(dz)/sizeof(double);
  int nz2 = sizeof(ispowz)/sizeof(int);
  assert(nz==nz2);
  geom->extrude(nz,zseg,dz,ispowz);
  /* set material index for fuel assemblies */
  /* 0: graphite,  1: gra. w. hole,  2: fuel block  */
  int matz[] = {1,1,1, 2,1,2,1,2,1, 2,1,2,1,2, 1,1,1};
  int nz3 = sizeof(matz)/sizeof(int);
  assert(nz3==nz);
  geom->setHexMat(1,matz);
  geom->setHexMat(3,matz);
  geom->setHexMat(4,matz);
  geom->setHexMat(5,matz);
  geom->setHexMat(6,matz);

  int ntri = geom->getNtri();

  XScapp *xs = new XScapp();
  xs->setGeom(geom);
  int efpd =50;
  double Tmref = 800;
  double Tfref = 800;
  xs->open("MMR1g.db",efpd,Tmref,Tfref);

  NeutFVM *neut = new NeutFVM();
  neut->setGeom(geom);
  neut->setXS(xs,Tmref,Tfref);
  neut->prepare(1);     // set no. of neutron groups

  neut->setEigShift(0.7);
  double sigrod[] = {0.0046,0.0046};
  neut->setRodXS(2, sigrod);
  double rodpos[] = {0.0,0.0};
  neut->setRodPos(2,rodpos);

  if(mpid==0) printf("#  Tm and Tf variation\n");
  int nvols = ntri*nz;
  double *Tm = new double[nvols];
  double *Tf = new double[nvols];
  if(mpid==0) printf("#  Tm Tf  keff\n");
  int nk = 10;
  for(int k=0; k <= nk; k++) {
    double Tfvar = 500+k*50;
    double Tmvar = Tmref;
    for(int n=0; n < nvols; n++) {
      Tm[n] = Tmvar;
      Tf[n] = Tfvar;
    }
    neut->setTemp(Tm,Tf);

    neut->setupAB();
    double keff = neut->solveEigen();
    if(mpid==0) printf("%.1lf %.1lf %.5lf\n",Tmvar,Tfvar,keff);
  }
  for(int k=0; k <= nk; k++) {
    double Tfvar = Tfref;
    double Tmvar = 500+k*50;
    for(int n=0; n < nvols; n++) {
      Tm[n] = Tmvar;
      Tf[n] = Tfvar;
    }
    neut->setTemp(Tm,Tf);

    neut->setupAB();
    double keff = neut->solveEigen();
    if(mpid==0) printf("%.1lf %.1lf %.5lf\n",Tmvar,Tfvar,keff);
  }
  for(int k=0; k <= nk; k++) {
    double Tfvar = 500+k*50;
    double Tmvar = 500+k*50;
    for(int n=0; n < nvols; n++) {
      Tm[n] = Tmvar;
      Tf[n] = Tfvar;
    }
    neut->setTemp(Tm,Tf);

    neut->setupAB();
    double keff = neut->solveEigen();
    if(mpid==0) printf("%.1lf %.1lf %.5lf\n",Tmvar,Tfvar,keff);
  }




  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
