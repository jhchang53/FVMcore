/*
 * 	segmult.cpp
 * 	segcond multiple run
 */
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include "SegCond.h"
#include "MPwork.h"

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

  SegCond *cond = new SegCond();
  cond->setGeom();
  /*  prepare multiple powden vector  */
  int npows = 5;
  Vec powdenvec,Tbulkvec;
  VecCreate(PETSC_COMM_WORLD,&powdenvec);
  VecSetSizes(powdenvec,PETSC_DECIDE,npows);
  VecSetUp(powdenvec); 
  VecSetOption(powdenvec,VEC_IGNORE_NEGATIVE_INDICES,PETSC_TRUE);
  VecDuplicate(powdenvec,&Tbulkvec);
  Vec powdennew;
  VecDuplicate(powdenvec,&powdennew);

  /*  create multiplier */
  Vec fvec;
  VecDuplicate(powdenvec,&fvec);
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

  // cond->setGeom1D(20,10,0.1);
  cond->setup(npows);
  double Tbulk = 800.0;		// deg-K
  double powden = 5.0e+6;	// power density of compact
  double hcoef = 1700.0;
  VecSet(Tbulkvec,Tbulk);
  VecSet(powdenvec,powden);
  cond->setHcoef(hcoef);
  /*  prepare Thermal property  */
  MPwork *mpwork = new MPwork();
  /* system pressure  */
  mpwork->set(3.0e+6);
  /*  triso dimension (in m) */
  double dKer = 500.0e-6;
  double tlay[] = {100.0,40.0,35.0,40.0};
  double PF = 0.4;
  mpwork->setTriso(dKer,4,tlay,PF);
  double flu = 1.0;
  double fima = 5.0;
  printf("# flu=%2lf fima=%.2lf\n",flu,fima);
  /*  prepare flu and fima in vector form  */
  Vec fluvec,fimavec;
  VecDuplicate(powdenvec,&fluvec);
  VecSet(fluvec,flu);
  VecDuplicate(powdenvec,&fimavec);
  VecSet(fimavec,fima);
  /*  prepare array of vectors  */
  Vec Tempvec[3],condvec[3],rhocpvec[3];
  for(int kv=0; kv < 3; kv++) {
    VecDuplicate(powdenvec,Tempvec+kv);
    VecDuplicate(powdenvec,condvec+kv);
    VecDuplicate(powdenvec,rhocpvec+kv);
  }

  Vec Tcomp_old,Terr;
  VecDuplicate(powdenvec,&Tcomp_old);
  VecDuplicate(powdenvec,&Terr);
  VecSet(Tcomp_old,0.0);
  Vec Tcomp;
  double err;
  for(int iter=0; iter < 8; iter++ ) {
    cond->steady(Tbulkvec,powdenvec);
    cond->dumpT();
    Tcomp = cond->getTpinVec();
    if(prlev) {
      if(mpid==0) printf("Tcomp:\n");
      VecView(Tcomp,PETSC_VIEWER_STDOUT_WORLD);
    }
    Vec Tgra = cond->getTgraVec();
    if(prlev) {
      if(mpid==0) printf("Tgra:\n");
      VecView(Tgra,PETSC_VIEWER_STDOUT_WORLD);
    }
    /*  update thermal property  */
    VecCopy(Tgra,Tempvec[0]);
    VecCopy(Tcomp,Tempvec[2]);
    VecWAXPY(Tempvec[1],1.0,Tgra,Tcomp);
    VecScale(Tempvec[1],0.5);
    mpwork->makeSegProp(fluvec,fimavec,3,Tempvec,condvec,rhocpvec);
    if(prlev) {
      for(int m=0; m < 3; m++) {
        if(mpid==0) printf("condvec[%d]:\n",m);
        VecView(condvec[m],PETSC_VIEWER_STDOUT_WORLD);
      }
    }
    cond->setCond(3,condvec,rhocpvec);
    VecCopy(Tcomp,Terr);
    VecAXPY(Terr,-1.0,Tcomp_old);
    VecNorm(Terr,NORM_1,&err);
    if(mpid==0) printf("== iter=%d err=%.1le\n",iter,err);
    VecCopy(Tcomp,Tcomp_old);
    if(err < 0.1) break;
  }
  /*  for TRISO  */
  
  /*  prepare array of vectors  */
  Vec condTvec[6],rhocpTvec[6];
  for(int kv=0; kv < 6; kv++) {
    VecDuplicate(powdenvec,condTvec+kv);
    VecDuplicate(powdenvec,rhocpTvec+kv);
  }
  mpwork->makeTriProp(fluvec,fimavec,Tcomp_old, 6,condTvec,rhocpTvec);
    for(int m=0; m < 6; m++) {
        if(mpid==0) printf("condTvec[%d]:\n",m);
        VecView(condTvec[m],PETSC_VIEWER_STDOUT_WORLD);
      }

  exit(0);

  /*  start transient */
  cond->start();
  double endtime = 100.0;
  double dtmin = 1.0e-6;
  double dtmax = 100.0;
  double dt = 0.001;
  double time  = 0.0;
  VecCopy(powdenvec,powdennew);
  // powden = 0.0;
  // VecPointwiseMult(powdennew,powdenvec,fvec);
  // printf("fvec:");
  // VecView(fvec,PETSC_VIEWER_STDOUT_WORLD);

  if(mpid==0) {
    printf("#  Tbulk=%.2lf powden=%.1le\n",Tbulk,powden);
    printf("# time   dt     rLTE Tpeak  Tpin\n");
  }
  for(int nt=0; (time < endtime) && (nt < 10000); nt++) {
    double rLTE = cond->step(dt,Tbulkvec,powdennew);
    cond->march();
    time += dt;
    if(mpid==0) printf("%.3le %.1le %.1le",time,dt,rLTE);
    double Tpeak = cond->getTpeak();
    double Tpin = cond->getTpin();
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
  Tcomp = cond->getTpinVec();
  VecView(Tcomp,PETSC_VIEWER_STDOUT_WORLD);

  auto end = std::chrono::high_resolution_clock::now();
  double xtime =std::chrono::duration_cast<std::chrono::milliseconds>(end-start).count();
  if(mpid==0) printf("== execution %.3lf sec\n",xtime/1000);

};
