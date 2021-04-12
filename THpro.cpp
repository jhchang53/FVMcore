/*
 * 	THpro.cpp
 * 	TH property for Cond3D
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "THpro.h"

THpro::THpro()
{
  prlev = 0;
  nmats = 0;
  PI = 2*acos(0.0);
  T0 = 273.15;  // STP
  P0 = 1.0e+5;  // STP

  fkweit = NULL;
  fcweit = NULL;
};	// THpro::THpro

THpro::~THpro()
{
  delete [] fkweit;
  delete [] fcweit;
};	// THpro::~THpro

void THpro::setPrlev(int prl)
{
  prlev = 0;
};	// THpro::setPrlev

void THpro::setChannel(double D_i, double wetden_i,
   double ncool_i)
{
  /*  Dcool : coolant hole diameter (m)                 */
  /*  wetden : wetted perimeter density (1/m)           */
  /*    multiplier to heat transfer coefficient         */
  /*  ncool : no. of coolant holes hannels in core plane             */
  /*    used to obtain flow speed in a channel          */
  D = D_i;        // Dcool is give in meter
  wetden = wetden_i;
  ncool = ncool_i;	// used to normalize mdot
};      // THprop::setChan

void THpro::setMdot(double Mdot_i, double Psys_i)
{
  Mdot = Mdot_i;
  Psys = Psys_i;
};	// THpro::setMdot

double THpro::hAcoef(double Tbulk)
{
  return wetden*hcoef(Tbulk);
};	// THpro::getHcoefA

double THpro::hcoef(double Tsys)
{
  /*  heat transfer coefficient (W/m^2/K)   */
  if(prlev) printf("Psys=%.2le Tsys=%.2lf\n",Psys,Tsys);
  assert(Tsys > 200.0);
  assert(Psys > 1.0e+5);
  /*  compute density */
  double rho = rho_He(Tsys,Psys);
  double vis = dynviscosity(Tsys);
  if(prlev) printf(" rho=%.3le vis=%.3le\n",rho,vis);
  double mdot = Mdot/ncool;
  if(prlev) printf(" ncool=%.1lf mdot=%.2le\n",ncool,mdot);
  double Darea = PI*D*D/4;
  double speed = mdot/rho/Darea;
  if(prlev) printf("speed=%.2lf (m/s)\n",speed);
  double Re = rho*speed*D/vis;
  double Pr = Prandtl(Tsys,Psys);
  double Nu = 0.021*pow(Re,0.8)*pow(Pr,0.3);
  double cond = cond_He(Tsys,Psys);
  double hCoef = cond/D*Nu;
  if(prlev) printf(" cond=%.2le (W/m/K) hcoef=%.2le (W/m^2/K)\n",cond,hCoef);
  if(prlev > 1) exit(0);
  return hCoef; // (W/m^2/K)
};      // THpro::hcoef
#ifdef XX
double THpro::density(double T, double P)
{
  /*  density (kg/m^3)                                  */
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double den = 0.17623*(P/P0)/(T/T0) / (1+0.53e-3*(P/P0)/pow(T/T0,1.2));
  return den;
};      // THpro::density
#endif

double THpro::dynviscosity(double T)
{
  /*  dynamic viscosity (kg/m/s)  */
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double dvis = 1.865e-5*pow(T/T0,0.7); // (kg/m/s)
  return dvis;
};      // THpro::dvis
#ifdef XX
double THpro::conductivity(double T, double P)
{
  /*  conductivity (W/m/K)  */
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double k = 0.144*(1+2.7e-4*(P/P0))*pow(T/T0,0.71*(1-2.0e-4*(P/P0)));
  return k;
};      // THpro::conductivity
#endif

double THpro::Prandtl(double T, double P)
{
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double pr = 0.6728/(1+2.7e-4*(P/P0))*pow(T/T0,-(0.01-1.42e-4*(P/P0)));
  return pr;
};      // THpro::Prandtl


/*	thermal conductivity (W/m/K) and heat capacity (J/m^3/K) */

void THpro::setFactors(int nmats_i, double fkweit_i[], double fcweit_i[])
{
  assert(nmats==0);
  nmats = nmats_i;
  fkweit = new double[nmats];	// conductivity multiplier
  fcweit = new double[nmats];	// heat capacity multiplier
  for(int m=0; m < nmats; m++) {
    fkweit[m] = fkweit_i[m];
    fcweit[m] = fcweit_i[m];
  }
};	// THpro::setFactors

double THpro::cond(int m, double TK, double flu)
{
  assert(m < nmats);
  double kcond = cond_IG110(TK,flu);
  return fkweit[m]*kcond;
};	// THpro::cond

double THpro::rhocp(int m, double TK)
{
  assert(m < nmats);
  double cp = cp_IG110(TK);
  double rcp = rho_gra(TK)*cp;  // (J/m^3/K)
  return fcweit[m]*rcp;
};	// THpro::rhocp
