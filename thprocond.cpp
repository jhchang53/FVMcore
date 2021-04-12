/*
 *  thprocond.cpp
 *  test driver for TH properties : conductivity
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "THpro.h"

int main()
{
  THpro *thpro = new THpro();
  double fkweit[] = {1.0,0.9,0.8};      // multipliers for conductivity
  double fcweit[] = {1.0,0.9,0.9};      // multipliers for heat capacity
  thpro->setFactors(3,fkweit,fcweit);
  double flu = 1.0;
  for(int nt=0; nt < 5; nt++) {
    double TK = 573.15 + 100*nt;
    double cond = thpro->cond(0,TK,flu);
    double rhocp = thpro->rhocp(0,TK);
    printf("T=%.1lf C k=%.3lf W/m/K rhocp=%.3le J/m^3/K\n",TK-273.15,cond,rhocp);
  }
  delete thpro;
};
