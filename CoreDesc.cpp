/*
 * 	CoreDesc.cpp
 * 	core description
 * 	calculator to derive thsuper parameters
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CoreDesc.h"

CoreDesc::CoreDesc()
{
  prlev = 0;
  rt3 = sqrt(3.0);
  PI = 2*acos(0.0);
};

CoreDesc::~CoreDesc()
{
};

void CoreDesc::setPrlev(int prl)
{
  prlev = prl;
};	// CoreDesc::setPrlev

void CoreDesc::setMMR(double Totpow_i, double Psys_i, double Tin_i, double Tout)
{
  Totpow = Totpow_i;
  Psys = Psys_i;
  Tin = Tin_i;
  /*  describe MMR core */
  /*  fuel block  */
  block_pitch = 0.322;	// pitch of graphite blocks
  block_height = 0.7;
  block_nfuel = 30;	// number of fuel blocks in a plane
  block_nlayer = 6;
  block_Dfuel = 0.02336;
  block_Dcool = 0.014;	// coolant hole diameter
  /*  fuel pellet  */
  pellet_OD = 0.023;
  pellet_height = 0.0252;
  double gap_compact = (block_Dfuel - pellet_OD)/2;
  pellet_stack_height = 0.65528;
  packing_fraction = 0.4;
  pellet_over = 0.0005;	// SiC over coating
  /*  triso  */
  double triso_enr = 19.75;
  double triso_c = 0.4;
  double triso_o = 1.5;
  triso_Dkernel = 500.0;	// we use um for kernel
  double triso_kernel_den = 10.40;	// (gr/cc)
  triso_OD = 930.0;	// including over coating
  triso_ncoat = 4;
  triso_tbuf = 100.0;	double triso_buf_den = 1.05;
  triso_tipyc = 40.0;	double triso_ipyc_den = 1.90;
  triso_tsic = 35.0;	double triso_sic_den = 3.19;
  triso_topyc = 40.0;	double triso_opyc_den = 1.90;
  /* fuel block  radial direction */
  no_pellets = 48;
  int no_bps = 6;
  no_holes = 21;	// no. of coolant holes per fuel block
  cell_pitch = 0.033333;

  double Cp_he = 5195.0;	// (J/kg/K)
  Mdot = Totpow/((Tout-Tin)*Cp_he);
  if(prlev) printf("== Mdot=%.3lf (kg/s)\n",Mdot);	// no. of coolant holes per fuel block	// no. of coolant holes per fuel block
  makeTH3D();
  makeMeso();
  makeSuper();
  makeTriso();
  int pel_in_block = (pellet_stack_height+1.0e-5)/pellet_height;
  if(prlev) printf(" pel_in_block=%d\n",pel_in_block);
  int pellets = block_nlayer*block_nfuel*pel_in_block*no_pellets;
  double ntriso = pellets*triso_compact;  // 8000;
  double pker = Totpow/ntriso;
  if(prlev) printf(" ntriso=%.2le pker=%.2le\n",ntriso,pker);
  if(prlev) printf(" %d compacts\n",pellets);
};	// CoreDesc::setMMR

void CoreDesc::makeTH3D()
{
  /*  global model parameters  */
  side = block_pitch/rt3;
  fblk = pellet_stack_height;
  fgra = block_height - pellet_stack_height;
  Dcool = block_Dcool;
  block_area = 3*rt3/2*side*side; 	// top area of a hexagon
  double bfheit = block_nlayer*pellet_stack_height;	// height of active fuel
  double volpow = bfheit*block_nfuel*block_area;
  if(prlev) printf(" Power volume=%.3lf m^3\n",volpow);
  powden = Totpow/volpow;
  /*  flow paramters  */
  nchan = no_holes*block_nfuel;
  /* wetted parameter density  */
  wet_den = no_holes*PI*block_Dcool/block_area;
  /*  compute effective conduction area  */
  double farea = no_pellets*PI*pellet_OD*pellet_OD/4;
  double carea = no_holes*PI*Dcool*Dcool/4;
  ffuel = farea/block_area;
  fgraf = (block_area-farea-carea)/block_area;
};	// CoreDesc::makeTH3D

void CoreDesc::makeMeso()
{
  rF = pellet_OD/2;
  rgap = block_Dfuel/2;
  double p = cell_pitch;
  cell_area = 3*rt3/2*p*p;
  double hole_area = PI*block_Dcool*block_Dcool/4;
  rD = sqrt((cell_area-hole_area)/PI);
  /*double f_qcomp,qcomp; // compact power density  */
  compact_area = PI*pellet_OD*pellet_OD/4;
  double bfheit = block_nlayer*pellet_stack_height;
  double volcomp = bfheit*no_pellets*block_nfuel*compact_area;
  qcomp = Totpow/volcomp;
  f_qcomp = block_nfuel*block_area/(no_pellets*block_nfuel*compact_area);
  if(prlev) printf("qcomp=%.3le calc=%.3le\n",qcomp,f_qcomp*powden);
};	//  CoreDesc::makeMeso

void CoreDesc::makeSuper()
{
  rC = block_Dcool/2;
  rA = cell_pitch - (block_Dcool+block_Dfuel)/2;
  rB = sqrt(cell_area/PI);
  /*  fuel fraction in mixture is averaged  */
  double mixarea = PI*(rB*rB-rA*rA);
  double farea = PI*pellet_OD*pellet_OD/4;
  ffrac = no_pellets*farea/(no_holes*mixarea);
  /*  power density in mixture region  */
  qmix = no_pellets*qcomp*farea/(no_holes*mixarea);
  f_qmix = f_qcomp*no_pellets*farea/(no_holes*mixarea);
};	// CoreDesc::makeSuper

void CoreDesc::makeTriso()
{
  dKer = triso_Dkernel;
  ncoat = triso_ncoat;
  tcoat[0] = triso_tbuf;
  tcoat[1] = triso_tipyc;
  tcoat[2] = triso_tsic;
  tcoat[3] = triso_topyc;
  tcoat[4] = triso_OD - triso_Dkernel
	 -2*(triso_tbuf+triso_tipyc+triso_tsic+triso_topyc);
  PF = packing_fraction;
  /*  power per kernel  */
  double tr = triso_OD/2*1.0e-6;
  double voltriso = 4*PI/3*tr*tr*tr;
  double cr = pellet_OD/2-pellet_over;
  double ch = pellet_height - 2*pellet_over;
  double volcompact =  PI*cr*cr*ch;
  triso_compact = PF*volcompact/voltriso;
  pker = volcompact*qcomp/triso_compact;
  f_pker = f_qcomp * volcompact/triso_compact;
};	// CoreDesc::makeTriso

void CoreDesc::print()
{
  printf("=== Global model parameters ===\n");
  printf(" Totpow=%.2le (W)\n",Totpow);
  printf(" Mdot=%.2le (kg/s)\n",Mdot);
  printf(" side=%.2le (m)\n",side);
  printf(" fblk=%.2le (m)\n",fblk);
  printf(" fgra=%.2le (m)\n",fgra);
  printf(" Dcool=%.2le (m)\n",Dcool);
  printf(" powden=%.2le (W/m^3)\n",powden);
  printf(" nchan=%.1lf\n",nchan);	// no. of coolant holes in plane
  printf(" wet_den=%.2le (1/m)\n",wet_den);	// density of wetted perimeter
  printf(" ffuel=%.3lf\n",ffuel);
  printf(" fgraf=%.3lf\n",fgraf);

  printf("=== meso model ===\n");
  printf(" rF  =%.5le (m)\n",rF);
  printf(" rgap=%.5le (m)\n",rgap);
  printf(" rD  =%.5le (m)\n",rD);
  printf(" f_qcomp=%.3le\n",f_qcomp);
  printf(" qcomp=%.3le (W/m^3) %.3le\n",qcomp,f_qcomp*powden);
  printf("=== super model ===\n");
  printf(" rC  =%.5le (m)\n",rC);
  printf(" rA=%.5le (m)\n",rA);
  printf(" rB  =%.5le (m)\n",rB);
  printf(" f_qmix=%.3le\n",f_qmix);
  printf(" qmix=%.3le (W/m^3)  %.3le\n",qmix,f_qmix*powden);
  printf(" ffrac=%.3lf\n",ffrac);
  printf("=== triso model ===\n");
  printf(" dKer  =%.5lf (um)\n",dKer);
  printf(" tlay:");
  for(int j=0; j < ncoat; j++) printf(" %.2lf",tcoat[j]);
  printf(" (um)\n");
  printf(" PF=%.3lf\n",PF);
  printf(" f_pker=%.3le\n",f_pker);
  printf(" pker=%.3le (W) %.3le\n",pker,f_pker*powden);
  printf(" # no. of triso in a compact = %.1lf\n",triso_compact);
};	// CoreDesc::print
