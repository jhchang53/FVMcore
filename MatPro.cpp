/*
 * 	MatPro.cpp
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "MatPro.h"

MatPro::MatPro()
{
  prlev = 0;
  T0 = 273.15;
  P0 = 1.0e+5;
};

MatPro::~MatPro()
{
};

double MatPro::cond0_IG110(double TK)
{
  /*  adopted from Graphite Datasheet, Ben Webster, USNC (2021)  */
  double TC = TK-T0;
  double ku = 130.04+TC*(-1.747e-1+TC*(1.447e-4+TC*(-6.214e-8+TC*1.132e-11)));
  return ku;
};

double MatPro::cond_IG110(double TK, double flu)
{
  double TC = TK-T0;
  double Tmin = 300.0;  double Tmax = 1100.0;
  double fmin = 0.0;    double fmax = 36.44;
  double Tin = (TC-Tmin)/(Tmax-Tmin);   // (4.4)
  double fn = (flu-fmin)/(fmax-fmin);   // (4.5)
  /*  Exh. 4-4  */
  double a1=-0.985;     double a2= -204.14;
  double b1= 0.441;     double b2=-1520.62;
  double c1=-0.261;     double c2=-1197.38;
  double d1= 0.709;     double d2=  291.72;
  double e1=-0.717;     double e2= 1655.56;
  double f1= 0.581;     double f2=  -0.457;
  double g1=-0.170;     double g2= -113.28;
  double h1= 0.156;     double h2= -0.038;
  double A = a1*fn*Tin*Tin + b1*fn*fn + c1*fn*fn*fn + d1*Tin*fn
        + e1*Tin*fn*fn + f1*Tin + g1*fn + h1;   // (4.6)
  double B = a2*fn*Tin*Tin + b2*fn*fn + c2*fn*fn*fn + d2*Tin*fn
        + e2*Tin*fn*fn + f2*Tin + g2*fn + h2;   // (4.7)
  double k0 = cond0_IG110(TK);
  double factor = A+pow(10.0,B);        // (4.3)
  double cond = k0*factor;
  if(cond < 0.1) {
    printf("***\n");
    printf("*** %s:%d cond=%.2le is too low.\n",__FILE__,__LINE__,cond);
    printf(" TK=%.2lf flu=%.2lf\n",TK,flu);
  }
  assert(cond > 0.1);
  return cond;
};      //  MatPro::cond_IG110

double MatPro::cp_IG110(double TK)
{
  /*  adopted from Graphite Datasheet, Ben Webster, USNC (2021)  */
  double TC = TK-T0;
  double cp = 0.638+TC*(3.05e-3+TC*(-3.24e-6+TC*(1.78e-9-3.78e-13*TC)));
  return 1000*cp;
};

double MatPro::rho_gra(double TK)
{
  /*  T.D.Burchell, CNM 4.10 (2012). Table 1 (kg/m^3) */
  double den = 1770; // IG-110 (kg/m^3) at room T
  double rho = den*(1-3*LTE_gra(TK)*(TK-300.0));
  return rho;
};

double MatPro::LTE_gra(double TK)
{
  /* T.D.Burchell, CNM 4.10 (2012). Table 2. */
  double lte = 4.5e-6;  // IG-110
  return lte;
};

double MatPro::cond_He(double T, double P)
{
  /*  He conductivity (W/m/K)  */
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double k = 0.144*(1+2.7e-4*(P/P0))*pow(T/T0,0.71*(1-2.0e-4*(P/P0)));
  return k;
};      // THpro::conductivity

double MatPro::cp_He(double T)
{
  /*  note. specific heat of gas is indept. of temperature  */
  return 5195.0;
};	//  MatPro::cp_He

double MatPro::rho_He(double T, double P)
{
  /*  density (kg/m^3)                                  */
  /*  ref. H.Petersen, The properties of helium, (1970) */
  double den = 0.17623*(P/P0)/(T/T0) / (1+0.53e-3*(P/P0)/pow(T/T0,1.2));
  return den;
};      // MatPro::rho_He

double MatPro::dens_UO2(double TK)
{
  /* density of solid UO2  */
  /* ref. TECDOC-1496, p.115, eqn (2), (3) (2006)       */
  double rho273 = 10960;        // (kg/m3)
  double LT;
  if(TK < 923.0) LT=9.9734e-1+TK*(9.802e-6+TK*(-2.705e-10+TK*4.291e-13));
  else LT=9.9672e-1+TK*(1.179e-5+TK*(-2.429e-9+TK*1.219e-12));
  double rho = rho273/(LT*LT*LT);
  return rho;
};      // MatPro::dens_UO2


double MatPro::cond_UCO(double TK)
{
  /* thermal conductivity of UCO (W/m/K)        */
  /*  adopted from BISON; http://mooseframework.inl.gove/bison/source/  */
  /*  ref. Nabielek et al, Calculation of Particle temperatures in  */
  /*    NSRR tests, JAEA (1992)                                         */
  double TC = TK-T0;
  double k;
  if(TC < 1650) k = 0.0132*exp(0.00188*TC) + 4040.0/(464.0+TC);
  else k = 0.0132*exp(0.00188*TC) + 1.9;
  return k;
};      // MatPro::cond_UCO

double MatPro::cp_UCO(double TK, double Umass, double Cfrac, double Ofrac)
{
  /*  specific heat of UCO (J/kg/K)  */
  /*  adopted from BISON; http://mooseframework.inl.gove/bison/source/  */
  /*  ref. Fink, Thermophysical properties of uranium oxide,    */
  /*    JNM279,1-18 (2000)      */
  double T = TK/1000;
  double MolarMass = (Umass+12*Cfrac+16*Ofrac); // (kg/mol)
  double cp = (52.1743+T*(87.951+T*(-84.2411+T*(31.542-2.6334*T)))
        - 0.71391/(T*T)) / MolarMass;
  return cp;
};      // MatPro::cp_UCO

double MatPro::cond_UCO_Ben(double TK, double fima, double por)
{
  /*  adopted from UCO datasheet 0.3, Ben Webster, USNC (2020)      */
  /*	TK : temperature (K)		*/
  /*	fima: burnup in %fima		*/
  /*    por: kernel porosity            */
  double tc = TK;
  double bu = fima;
  double k0 = 1.0/(0.0375+(2.165e-4)*tc)+(4.715e+9/(tc*tc))*exp(-16361/tc);
  double k1d = ((1.09/pow(bu,3.265))+(0.0642/sqrt(bu))*sqrt(tc)
	* atan(1/((1.09/pow(bu,3.265))+(0.0642/sqrt(bu))*sqrt(tc))));
  double k1p = 1+(((0.019*bu)/(3-(0.019*bu)))*(1/(1+exp(-(tc-1200)/100))));
  double k4r = 1-(0.2/(1+exp(-(tc-900)/80)))*(1-exp(-bu));
  double k2p = (1-por)/(1+(2*por));
  double k = k0 * k1d * k1p * k4r * k2p;
  return k;
};	// MatPro::cond_UCO_Ben

double MatPro::cond_Buf(double rho)
{
  /*  rho : kg/m3  */
  /*  adopted from BISON; http://mooseframework.inl.gove/bison/source/  */
  /*  ref. Miller et al., PARFUME Theory and Model Basis Report.        */
  /*  INL/EXT-08-14497 (Rev.1) (2018)  */
  double k_init = 0.5;  double rho_init = 1000.0;
  double k_theo = 4.0;  double rho_theo = 2250.0;
  double k = k_init*k_theo*rho_theo*(rho_theo-rho_init) /
        (k_theo*rho_theo*(rho_theo - rho)+k_init*rho*(rho_init));
  return k;
};      //  MatPro::cond_Buf

double MatPro::cp_Buf()
{
  /* adopted from BISONi http://mooseframework.inl.gove/bison/source/  */
  /* ref. Barabash et al., The effect of low temperature neutron irradiation */
  /* and annealing on the thermal conductivity of advanced carbon-based */
  /* materials, JNM 307-311:1300-1304, (2002)   */
  double cp = 720.0;    // (J/kg/K)
  return cp;
};      // MatPro::cp_Buf

double MatPro::rho_Buf()
{
  return 1000.0;
};

double MatPro::cond_PyC(double P)
{
  /* conductivity (W/m/K)  */
  /*  P : porosity  */
  /*  adopted Powers, JNM 405 (2010) 74-82  */
  /*  ref. INEEL-EXT-05-02615, Table 1-9, UK  */
  double kpyc = 10.98222*(1-P)/(1+2*P) + 0.00444;
  return kpyc;
};      //  MatPro::cond_PyC

double MatPro::rho_PyC()
{
  return 2000.0;
};	// MatPro::rho_PyC

double MatPro::cp_PyC(double TK)
{
  /*  adopted heat capacity of graphite (J/kg/K)	*/
  /*  ref. T.D.Burchell, CNM 2.10. 286-304 (2012)	*/
  double cp = 1.0/(11.07*pow(TK,-1.644)+0.0003688*pow(TK,0.02191));
  return cp;
};	// MatPro::cp_PyC

double MatPro::cond0_SiC(double TK)
{
  /*  adopted Powers, JNM 405 (2010) 74-82  */
  /*  ref) INEEL-EXT-05-02615, p.31, CEA        */
  double k0 = 42.58+ ((-1.8458e+9/TK + 1.2977e+7)/TK-1.5564e+4)/TK;
  return k0;
};

double MatPro::cond_SiC(double TK, double P)
{
  /*  adopted Powers, JNM 405 (2010) 74-82  */
  /*  ref) INEEL-EXT-05-02615, p.31, CEA        */
  double ksic = cond0_SiC(TK)*(3.91112e-2*exp(2.24732e-3*TK)*(1-P));
  return ksic;
};      // MatPro::cond_SiC

double MatPro::cond_SiC_Ben(double TK, double TI, double flu)
{
  /*  adopted from SiC datasheet 0.3, Ben Webster, USNC (2020)   */
  /*	ti: irradiation temperature				*/
  /*	tc: current material temperature			*/
  double tc = TK-T0;
  double ti = TI-T0;
  double y = flu;
  double ku = 67.6*exp(-0.00142*tc)+24.33;
  double rd = 0.221*exp(-0.00205*ti);
  double kis = 1/((1/ku)+rd);
  double y_mid = 0.00013*exp(0.00801*ti) + 0.00291;
  double y43 = pow(y,1.3333);
  double y_mid43 = pow(y_mid,1.33333);
  double k = (ku+kis)/2-(((ku-kis)/2)*((y43-y_mid43)/(y43+y_mid43)));
  return k;
};	// MatPro::cond_SiC_Ben


double MatPro::rho_SiC()
{
  /*  ref) Snead, JNM 371 (2007) 329-377. Theoretical density  */
  return 3210;
};

double MatPro::cp_SiC(double TK)
{
  /*  ref) Snead, JNM 371 (2007) 329-377  */
  double cp = 925.65 + TK*(0.3772-7.9259e-5*TK) - 3.1946e-7/(TK*TK);
  return cp;
};      //  MatPro::cp_SiC

double MatPro::cond0_NITE(double TK)
{
  /*  unirradiated NITE  */
  /* adopted from NITE SiC datasheet 0.3, Ben Webster, USNC (2020)      */
  double TC = TK-T0;
  double ku = 18.73*exp(-1.64e-3*TC) + 14.09;
  return ku;
};      // MatPro::cond0_NITE

double MatPro::cond_NITE_Ben(double TK, double TI, double flu)
{
  /*  irradiated NITE , flu ( 0.1 MeV)  */
  /*  adopted from NITE SiC datasheet 0.3, Ben Webster, USNC (2020)      */
  double ti = TI-T0;
  double ku = cond0_NITE(TK);
  double rd = 2.31e-1*exp(-1.58e-3*ti);
  double kis = 1.0/(1.0/ku+rd);
  // double flum = 2.74e-4*exp(7.98e-3*ti)+7.66e-3;
  double flum = 2.32e-4*exp(7.98e-3*ti) + 0.00647;
  double fl43 = pow(flu,4.0/3.0);
  double fm43 = pow(flum,4.0/3.0);
  double kirr = (ku+kis)/2 - (ku-kis)/2*(fl43-fm43)/(fl43+fm43);
  return kirr;
};      // MatPro::cond_NITE_Ben

