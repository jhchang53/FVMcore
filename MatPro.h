#ifndef _inc_MatPro
#define _inc_MatPro
/*
 * 	fundamental class for Material Properties
 */
class MatPro {
public:
  MatPro();
  ~MatPro();
  double cond_IG110(double TK, double flu);
  double cp_IG110(double TK);
  double rho_gra(double TK);
  double cond_He(double TK, double Psys);
  double cp_He(double TK);
  double rho_He(double TK, double Psys);
  double dens_UO2(double TK);
  double cond_UCO(double TK);
  double cond_UCO_Ben(double TK, double fima, double por);
  double cp_UCO(double TK, double Umass, double Cfrac, double Ofrac);
  double cond_Buf(double rho);
  double cp_Buf();
  double rho_Buf();
  double cond_PyC(double por);
  double rho_PyC();
  double cp_PyC(double TK);
  double cond_SiC(double TK, double por);
  double cond_SiC_Ben(double TK, double TI, double flu);
  double rho_SiC();
  double cp_SiC(double TK);
  double cond_NITE_Ben(double TK, double TI, double flu);
private:
  int prlev;
  double T0,P0;
  double cond0_IG110(double TK);
  double LTE_gra(double TK);
  double cond0_SiC(double TK);
  double cond0_NITE(double TK);
};
#endif
