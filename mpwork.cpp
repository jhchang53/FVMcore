/*
 * 	mpwork.cpp
 * 	prepare material properties for SegCond and TriCond
 */
#include <stdio.h>
#include <stdlib.h>
#include "MPwork.h"

int main()
{
  MPwork *work = new MPwork();
  /* system pressure  */
  work->set(3.0e+6);
  /*  triso dimension (in m) */
  double dKer = 500.0e-6;
  double tlay[] = {100.0,40.0,35.0,40.0};
  double PF = 0.4;
  work->setTriso(dKer,4,tlay,PF);
  double flu = 1.0;
  double fima = 5.0;
  printf("# flu=%2lf fima=%.2lf\n",flu,fima);
  double TK[] = {600.0,600.0,600.0};
  double cond[6],rhocp[6];
  printf("#  TK  Gra  k  rhocp     He   k  rhocp     Pin  k  rhocp\n");
  for(int nt=0; nt < 10; nt++) {
    double tk = 600.0 + 50*nt;
    TK[0] = TK[1] = TK[2] = tk;
    work->makeSegProp(flu,fima,3,TK,cond,rhocp);
    printf("%.1lf",TK[0]);
    for(int m=0; m < 3; m++) printf(" %.2le %.2le",cond[m],rhocp[m]);
    printf("\n");
  }
  /*  TRISO properties  */
  printf("#  TRISO properties ker/buf/iPyC/SiC/oPyC/mat  \n");
  for(int nt=0; nt < 10; nt++) {
    double tk = 600.0 + 50*nt;
    work->makeTriProp(flu,fima,tk,6,cond,rhocp);
    printf("%.1lf",TK[0]);
    for(int m=0; m < 6; m++) printf(" %.2le %.2le",cond[m],rhocp[m]);
    printf("\n");
  }

  delete work;
}
