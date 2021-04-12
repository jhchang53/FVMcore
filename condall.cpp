/*
 * 	condall.cpp
 * 	core-wide conduction usingFVM
 * 	combined with SegCond and TriCond
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <chrono>
#include <petscksp.h>
#include "GeomCAPP.h"

#include "Cond3D.h"
#include "SegCond.h"
#include "TriCond.h"

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

  int prlev = 0;

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
  int nz1 = sizeof(dz)/sizeof(double);
  int nz2 = sizeof(ispowz)/sizeof(int);
  assert(nz1==nz2);
  geom->extrude(nz1,zseg,dz,ispowz);
  /* set material index for fuel assemblies */
  /* 0: graphite,  1: gra. w. hole,  2: fuel block  */
  int matz[] = {1,1,1, 2,1,2,1,2,1, 2,1,2,1,2, 1,1,1};
  int nz3 = sizeof(matz)/sizeof(int);
  assert(nz3==nz1);
  geom->setHexMat(1,matz);
  geom->setHexMat(3,matz);
  geom->setHexMat(4,matz);
  geom->setHexMat(5,matz);
  geom->setHexMat(6,matz);

  Cond3D *cond = new Cond3D();
  cond->setGeom(geom);
  double Mdot = 8.7;	// kg/s
  double six = 6;

  THpro *thpro = new THpro();
  double fkweit[] = {1.0,0.9,0.8};      // multipliers for conductivity
  double fcweit[] = {1.0,0.9,0.9};      // multipliers for heat capacity
  thpro->setFactors(3,fkweit,fcweit);
  double Psys = 3.0e+6;
  double Dcool=1.40e-02; double wetden=9.31;
  int ncool= 570;	// number of coolant holes
  thpro->setChannel(Dcool,wetden,ncool);
  thpro->setMdot(Mdot,Psys);

  cond->setTHpro(thpro);

  cond->setMdot(Mdot/six);
  double corevol = geom->getPowVol();
  double Tin = 600.0;
  cond->setTin(Tin);
  cond->setup();
  int npows = geom->getNpows();
  Vec powden;
  ierr = VecCreate(PETSC_COMM_WORLD,&powden);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetSizes(powden,PETSC_DECIDE,npows);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetUp(powden);      CHKERRABORT(PETSC_COMM_SELF,ierr);
  ierr = VecSetOption(powden,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
        CHKERRABORT(PETSC_COMM_SELF,ierr);
  double totpow = 15.0e+6/six;
  double pden = totpow/corevol;
  VecSet(powden,pden);
  double Mdotcp = Mdot/six * 5195.0;
  if(mpid==0) printf(" Total vol=%.2le power=%.3le (W)  dT=%.2lf K\n",
	corevol,totpow,totpow/Mdotcp);
  
  double errth = 1000.0;
  for(int iter=0; iter < 4; iter++) {
    cond->steady(powden);
    errth = cond->updateTemp();
    if(mpid==0) printf("== errth=%.2le\n",errth);
    if(errth < 1.0) break;
  }
  if(errth > 1.0) {
    printf("*** thermal iteration not converged. errth=%.2le\n",errth);
    exit(0);
  }
  cond->printAx(3);
  cond->printAx(4);
  double Tpeak = cond->getTpeak();
  double TPavg = cond->getTPavg();
  double Texit = cond->getTexit();
  if(mpid==0) printf("== Tpeak=%.3lf Texit=%.3lf\n",Tpeak,Texit);

  /*  SegCond  */
  Vec hcoefvec = cond->getPhcoef();
  Vec Tbulkvec = cond->getPTbulk();

  SegCond *seg = new SegCond();
  seg->setGeom();
  seg->setup(npows);
  seg->setPfactor(4.5);     // multiplier to convert core power density to compact
  seg->setHcoef(hcoefvec);
  seg->steady(Tbulkvec,powden);
  Vec Tcomp = seg->getTpinVec();
  if(prlev > 1) {
    printf("== Tcomp:\n");
    VecView(Tcomp,PETSC_VIEWER_STDOUT_WORLD);
  }
  double Tpin = seg->getTpin();
  double TCpeak = seg->getTpeak();
  printf("== Tpin=%.2lf TCpeak=%.2lf\n",Tpin,TCpeak);
  

  /*  Triso  */
  TriCond *triso = new TriCond();
  triso->setGeom();
  triso->setup(npows);
  triso->setPfactor(4.74e-9);	// convert core power density to kernel power
  triso->steady(powden);
  double Tker = triso->getTker();
  printf(" Tker=%.5lf\n",Tker);
  Vec DTker = triso->getDTker();
  if(prlev > 1) {
    printf("== DTker:\n");
    VecView(DTker,PETSC_VIEWER_STDOUT_WORLD);
  }
  Vec Tfuel;
  VecDuplicate(powden,&Tfuel);
  VecCopy(Tcomp,Tfuel);
  VecAXPY(Tfuel,1.0,DTker);
  if(prlev > 1) {
    printf("== Tfuel:\n");
    VecView(Tfuel,PETSC_VIEWER_STDOUT_WORLD);
  }

  /*  start transient calculation  */
  cond->start();
  seg->start();
  triso->start();
  // cond->setTin(800.0);
  // cond->setMdot(0.0);
  VecSet(powden,0.0);
  // VecScale(powden,2.0);
  double endtime = 1000.0;
  double dtmin = 1.0e-6;
  double dtmax = 10000.0;
  double dt = 1.0e-6;
  double time  = 0.0;
  if(mpid==0) {
     printf("# time   dt    cLTE  Tpeak  TPavg  Texit");
     printf("   sLTE TCpeak Tpin");
     printf("   tLTE Tker\n");
  }
  for(int nt=0; (time < endtime) && (nt < 10000); nt++) {
    double cLTE = cond->step(dt,powden);
    cond->march();
    /* we advance time without worrying to recalculation */
    time += dt;
    if(mpid==0) printf("%.3le %.1le",time,dt);
    Tpeak = cond->getTpeak();
    TPavg = cond->getTPavg();
    Texit = cond->getTexit();
    printf("  %.1le %.2lf %.2lf %.2lf",cLTE,Tpeak,TPavg,Texit);

    double sLTE = seg->step(dt,Tbulkvec,powden);
    seg->march();
    double TCpeak = seg->getTpeak();
    double Tpin = seg->getTpin();
    printf("  %.1le %.2lf %.2lf",sLTE,TCpeak,Tpin);

    double tLTE = triso->step(dt,powden);
    triso->march();
    double Tker = triso->getTker();
    printf("  %.1le %.3lf",tLTE,Tker);

    if(mpid==0) printf("\n");
    /*  control time step  */
    /*  select maximum LTE to control time step size  */
    double rLTE = tLTE;
    if(rLTE < sLTE) rLTE = sLTE;
    if(rLTE < cLTE) rLTE = cLTE;
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

