/*
 *	cappkin.cpp
 *	test with delayed neutron
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
#include "Delay.h"

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

  XScapp *xs = new XScapp();
  xs->setGeom(geom);
  int efpd = 50;
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
  neut->setupAB();
  double keff = neut->solveEigen();
  if(mpid==0) printf("keff= %.5lf\n",keff);

  double dlam[] = {1.33373E-02, 3.27241E-02, 1.20789E-01,
                   3.03011E-01, 8.50566E-01, 2.85639E+00};
  double dbeta[]= {2.33250E-04, 1.20779E-03, 1.15433E-03,
                   2.59444E-03, 1.07174E-03, 4.48447E-04};
  double gvel[] = {1.042574e+4};
  int nprec = sizeof(dbeta)/sizeof(double);
  double beta = 0;
  for(int d=0; d < nprec; d++) beta += dbeta[d];
  double six = 6.0;
  neut->start(Totpow/six,keff,beta,gvel);
  // neut->makePower(Totpow/six);	// need to set for delayed neutron

  Vec fissvec = neut->getFiss();
  Delay *delay = new Delay();
  int npows = geom->getNpows();
  delay->setParam(npows,nprec,dbeta,dlam);
  delay->init(fissvec);
  Vec dnsrc = delay->getNeutron();
  neut->setDsrc(dnsrc);

  neut->setLTE(1.0e-5,1.0e-10);
  double endtime = 1000.0;
  double dtmin = 1.0e-6;
  double dtmax = 100.0;
  double dt = dtmin;
  double time = 0;
  rodpos[0] = 6.0;
  rodpos[1] = 6.0;
  // neut->setRodPos(2,rodpos);
  printf("# time dt  rLTE pow  peakden\n");
  for(int nt=0; (time < endtime) && (nt < 10000); nt++) {
    double rLTE = neut->step(dt);
    neut->march();
    time += dt;
    printf("%.3le %.1le %.1le",time,dt,rLTE);
    double peakden = neut->getPeakPowDen();
    double totpow = neut->getTotPow();
    printf("  %.6le %.2le", totpow,peakden);
    fissvec = neut->getFiss();
    delay->step(dt,fissvec);
    delay->march();
    dnsrc = delay->getNeutron();
    neut->setDsrc(dnsrc);

    printf("\n");
    /*  control time step  */
    double hstar = dt*pow(fabs(rLTE),-0.333333);
    if(hstar < 0.5*dt) dt = hstar;
    if(hstar > 2.0*dt) dt = hstar;
    if(dt < dtmin) dt = dtmin;
    if(dt > dtmax) dt = dtmax;
    double remtime = endtime-time;
    if(dt > remtime) dt = remtime;
  }


  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
