/*
 * 	triso conduction problem
 * 	using spherical FVM
 */
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
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

  TriCond *cond = new TriCond();
  cond->setGeom();
  double fpker = 4.74e-9;	// power density to kernel power multiplier
  cond->setup(1);
  double pker = 1.0e-2;	// power per kernel
  pker = 15.0e+7*fpker;
  if(mpid==0) printf("== pker=%.2le mW\n",pker*1000);
  cond->steady(pker);
  double Tker = cond->getTker();
  if(mpid==0) printf("Tker=%.5lf\n",Tker);
  cond->dumpT();
 
  /*  start transient */
  cond->start();
  double endtime = 100.0;
  double dtmin = 1.0e-6;
  double dtmax = 10.0;
  double dt = 1.0e-5;
  double time  = 0.0;
  // pker = 2.0e-1;
  if(mpid==0) {
    printf("#  pker=%.1le\n",pker);
    printf("# time   dt     rLTE Tpeak\n");
  }
  for(int nt=0; (time < endtime) && (nt < 20000); nt++) {
    double rLTE = cond->step(dt,pker);
    cond->march();
    time += dt;
    if(mpid==0) printf("%.3le %.1le %.1le",time,dt,rLTE);
    Tker = cond->getTker();
    if(mpid==0) {
      printf(" %.5lf",Tker);
      printf("\n");
    }
    if(Tker < 1.0e-5) break;
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

  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

  delete cond;
  PetscFinalize();
};
