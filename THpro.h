/*
 * 	THpro.h
 * 	TH property for Cond3D
 * 	
 *    Usage:
 *	 THpro();   - constructor
 *	 ~THpro();  - destructor
 *	 void setChannel(double D, double wetden, double ncool);
 *	   specifies coolant hole dimension to compute heat transfer coef.
 *	   D : coolant hole diameter (m)
 *	   wetden : wetted perimeter density of fuel blocks (m/m^2)
 *	   ncool : number of coolant holes to normalize Mdot
 *	 void setMdot(double Mdot, double Psys);
 *	   Mdot : mass flow rate (kg/s)
 *	   Psys : system pressure (Pa)
 *	 double hAcoef(double TK);
 *	   returns areal heat transfer coef. for porous media model (W/m^3/K)
 *	 double hcoef(double Tsys);
 *	   returts heat transfer coefficient (W/m^2/k)
 *	 void setFactors(int nmats, double fkweit[], double fcweit[]);
 *	   set multiplers to conductivity of block type
 *	 double cond(int m, double TK, double flu);
 *	   returns conductivity (W/m/K)
 *	   m : block type
 *	   TK : temperature (K)
 *	   flu : fast fluence (10^25 n/m^2)
 *	 double rhocp(int m, double TK);
 *	   returns volumetric heat capacity (J/m^3/K)
 *
 */
#include "MatPro.h"

class THpro : public MatPro{
public:
  THpro();
  ~THpro();
  void setPrlev(int prl);
  void setChannel(double D, double wetden, double ncool);
  void setMdot(double Mdot, double Psys);
  double hAcoef(double TK);
  double hcoef(double Tsys);
  void setFactors(int nmats, double fkweit[], double fcweit[]);
  double cond(int m, double TK, double flu);
  double rhocp(int m, double TK);

private:
  /* heat transfer routines */
  // double density(double T, double P);
  double dynviscosity(double T);
  // double conductivity(double T, double P);
  double Prandtl(double T, double P);

  int prlev,mpid;
  double PI,P0,T0;
  double Psys,Mdot;

  int nmats;
  double *fkweit,*fcweit;	// mutipliers to th. conductivity and heat cap.
  double D,wetden,ncool;
};
