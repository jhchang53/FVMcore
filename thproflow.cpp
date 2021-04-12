/*
 * 	thproflow.cpp
 * 	test driver for TH properties heat transfer coefficient
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "THpro.h"

int main()
{
  THpro *thpro = new THpro();
  double Psys = 3.0e+6;	// system pressure (Pa)
  double mdot = 8.691;  // kg/sec
  int nasy = 30;
  int nhole_block = 19;
  double Dcool = 0.014; // m
  double mdot_tri = mdot/(nasy*6);
  printf("mdot_tri=%.2le\n",mdot_tri);
  double pitch = 0.322; // m
  double side = pitch/sqrt(3.0);
  double PI = 2*acos(0.0);
  double area_blk = 3*sqrt(3.0)/2*side*side;
  double wetden = nhole_block*PI*Dcool/area_blk;
  double ncool = nasy*nhole_block;
  printf(" D=%.2le wetden=%.2le ncool=%.2le\n",
	Dcool,wetden,ncool);
  thpro->setChannel(Dcool,wetden,ncool);
  thpro->setMdot(mdot,Psys);
  double totpow = 15.0e+6;
  double powvol = 6*1.682;
  double powden = totpow/powvol;
  printf(" powden=%.3le (W/m3)\n",powden);
  for(int nt=0; nt < 5; nt++) {
    double TK = 573.15 + 100*nt;
    double hcoef = thpro->hcoef(TK);
    double hAcoef = thpro->hAcoef(TK);
    printf("%.1lfC h=%.2le W/m^2/K hA=%.2le W/m^3/K ", 
      TK-273.14,hcoef, hAcoef);
    printf("\n");
  }
  delete thpro;
};

