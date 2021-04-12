#ifndef _inc_MPwork
#define _inc_MPwork
/*
 * 	MPwork.h
 * 	material properties for SegCond and TriCond
 *
 *    usage:
 *      void set(double Psys);
 *        Psys: system pressure (Pa) required to compute He gap heat capacity
 *      void setTriso(double dKer, int nlays, double tlay[], double PF);
 *       set Triso dimensions (in meter)
 *        dKer : kernel diameter
 *        nlays : number of coating layers (currently = 4)
 *        tlay[] : coating thickness
 *        PF :  packing fraction (required to compute equiv. matrix volume)
 *      void makeSegProp(double flu, double fima, int nmats, double TK[],
 *        double cond[], double rhocp[]);
 *       compute material properties required for SegCond
 *        flu : fast neutron fluence (10^25 /m^2)
 *        fima : percent fima. used for kernel property
 *        nmats : number of material zone. current 3 (graphite/he gap/compact)
 *        TK : temperature of each zone (K)
 *       output:
 *        cond : array of conductivity (W/m/K)
 *        rhocp : array of heat capacity (J/m^3/K)
 *      void makeSegProp(Vec fluvec, Vec fimavec, int nmats, Vec Tempvec[],
 *        Vec condvec[], Vec rhocpvec[]);
 *       vector version. each vector is same size (npows)
 *       output:
 *         condvec, rhocpvec
 *      void makeTriProp(double flu, double fima, double TK,
 *        int nmats, double cond[], double rhocp[]);
 *       compute material properties required for TriCond
 *         flu : fast neutron fluence (10^25 /m^2)
 *         fima : percent fima. used for kernel property
 *         TK : temperature (K). compute material property using single temperature.
 *         nmats : number of material zone. current 6 (ker/buf/iPyC/SiC/oPyC/NITE)
 *       output:
 *         cond : array of conductivity (W/m/K)
 *         rhocp : array of heat capacity (J/m^3/K)
 *      void makeTriProp(Vec fluvec, Vec fimavec, Vec TKvec,
 *        int nmats, Vec condvec[], Vec rhocpvec[]);
 *       vector version
 *       output:
 *         condvec, rhocpvec
 *
 *   note) it need lapack and petsc.
 *      MPwork.o: MPwork.cpp
 *          mpicxx -c $^ $(LAINC) $(PETSC_CC_INCLUDES)
 *
 */
#include <petscksp.h>
#include "MatPro.h"

class MPwork : public MatPro {
public:
  MPwork();
  ~MPwork();
  void set(double Psys);
  void setTriso(double dKer, int nlays, double tlay[], double PF);
  void makeSegProp(double flu, double fima, int nmats, double TK[],
	double cond[], double rhocp[]);
  void makeSegProp(Vec fluvec, Vec fimavec, int nmats, Vec Tempvec[],
	Vec condvec[], Vec rhocpvec[]);
  void makeTriProp(double flu, double fima, double TK,
	int nmats, double cond[], double rhocp[]);
  void makeTriProp(Vec fluvec, Vec fimavec, Vec TKvec,
	int nmats, Vec condvec[], Vec rhocpvec[]);
private:
  void setCondL(double flu, double fima, double TK);
  double effCond();
  double effRhoCp();

  int prlev;
  double Psys;
  /* triso quantities */
  int nzone;
  double TK_triso,flu_triso,fima_triso;
  double PF;
  double *zrad,*condL,*rhocpL;
  double *Mat_keff,*b_keff;
  int *ipiv_keff;
};
#endif
