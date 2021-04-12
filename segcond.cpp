/*
 * 	segcond.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
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

  SegCond *cond = new SegCond();
  cond->setGeom();
  double fqcomp = 4.5/10.56; // core power to compact power density
  
  // cond->setGeom1D(20,10,0.1);
  cond->setup(1);
  double Tbulk = 800.0;		// deg-K
  double hcoef = 1700.0;
  double powden = 15.0e+6*fqcomp;
  powden = 5.0e+6;
  if(mpid==0) printf(" compact powden=%.2le W/m^3\n",powden);
  cond->setHcoef(hcoef);
  cond->steady(Tbulk,powden);
  double Tpin = cond->getTpin();
  if(mpid==0) printf("== Tpin=%.2lf\n",Tpin);
  /*  start transient */
  cond->start();
  double endtime = 100.0;
  double dtmin = 1.0e-6;
  double dtmax = 100.0;
  double dt = 0.001;
  double time  = 0.0;
  double pownew = powden;
  // Tbulk = 0.0;
 if(mpid==0) {
    printf("#  Tbulk=%.2lf pownew=%.1le\n",Tbulk,pownew);
    printf("# time   dt     rLTE Tpeak\n");
  }
  for(int nt=0; (time < endtime) && (nt < 10000); nt++) {
    cond->setHcoef(hcoef);
    double rLTE = cond->step(dt,Tbulk,pownew);
    cond->march();
    time += dt;
    if(mpid==0) printf("%.3le %.1le %.1le",time,dt,rLTE);
    double Tpeak = cond->getTpeak();
    Tpin = cond->getTpin();
   if(mpid==0) {
      printf(" %.2lf %.2lf",Tpeak,Tpin);
      printf("\n");
    }
    /*  control time step  */
    double hstar = dt*pow(fabs(rLTE),-0.333333);
    if(hstar < 0.5*dt) dt = hstar;
    if(hstar > 2.0*dt) dt = hstar;
    if(dt < dtmin) dt = dtmin;
    if(dt > dtmax) dt = dtmax;
    double remtime = endtime-time;
    if(dt > remtime) dt = remtime;
  }
  cond->dumpT();
  Tpin = cond->getTpin();
  if(mpid==0) printf("== Tpin=%.2lf\n",Tpin);

  delete cond;

  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

  PetscFinalize();
};
