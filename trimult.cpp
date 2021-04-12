/*
 * 	trimult.cpp
 * 	tricond multiple run
 *
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

  TriCond *triso = new TriCond();
  triso->setGeom();
  /*  prepare multiple pker vector */
  int npows = 4;
  Vec pkervec;
  VecCreate(PETSC_COMM_WORLD,&pkervec);
  VecSetSizes(pkervec,PETSC_DECIDE,npows);
  VecSetUp(pkervec);
  VecSetOption(pkervec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  Vec pkernew;
  VecDuplicate(pkervec,&pkernew);

  /*  create multiplier */
  Vec fvec;
  VecDuplicate(pkervec,&fvec);
  VecSet(fvec,0.0);
  int *indr =  new int[npows];
  double *val = new double[npows];
  double dpow = 1.0/npows;
  for(int kp=0; kp < npows; kp++) {
    indr[kp] = kp;
    val[kp] = kp*dpow;
  }
  VecSetValues(fvec,npows,indr,val,INSERT_VALUES);
  VecAssemblyBegin(fvec);
  VecAssemblyEnd(fvec);


  triso->setup(npows);
  double pker = 1.0e-2;	// power per kernel
  double fpker = 4.74e-9;       // power density to kernel power multiplier
  pker = 15.0e+6*fpker;
  printf("== pker=%.2le mW\n",pker*1000);

  VecSet(pkervec,pker);
  triso->steady(pkervec);
  double Tker = triso->getTker();
  if(mpid==0) printf(" Tker=%.5lf\n",Tker);
  // triso->dumpT();
  Vec dTker = triso->getDTker();
  VecView(dTker,PETSC_VIEWER_STDOUT_WORLD);

  /*  start transient */
  triso->start();
  double endtime = 100.0;
  double dtmin = 1.0e-6;
  double dtmax = 10.0;
  double dt = 1.0e-5;
  double time  = 0.0;
  
  VecSet(pkervec,pker);
  if(mpid==0) printf("#  pker=%.1le\n",pker);
  if(mpid==0) printf("# time    dt      rLTE    Tpeak\n");
  for(int nt=0; (time < endtime) && (nt < 2000); nt++) {
    double rLTE = triso->step(dt,pkervec);
    triso->march();
    time += dt;
    if(mpid==0) printf("%.3le %.1le %.1le",time,dt,rLTE);
    Tker = triso->getTker();
    if(mpid==0) printf(" %.5lf",Tker);
    if(mpid==0) printf("\n");
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
  // triso->dumpT();
  dTker = triso->getDTker();
  VecView(dTker,PETSC_VIEWER_STDOUT_WORLD);

  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
