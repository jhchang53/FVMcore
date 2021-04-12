/*
 * 	XScapp.h
 */
#include "XS.h"
#include "GeomCAPP.h"
#include "sqlite3.h"


class XScapp : public XS {
public:
  XScapp();
  ~XScapp();
  void setPrlev(int prl);
  void setGeom(GeomCAPP *geom);
  void open(const char *dbfn, int efpd, double Tm_ref, double Tf_ref);
  void getSig(int tz, double Tm, double Tf, double sig[]);
  double getSigtr(int tz, double Tm, double Tf);
  double getSigab(int tz, double Tm, double Tf);
  double getSignf(int tz, double Tm, double Tf);
  double getFlu(int tz);
  double *getFlu() { return flu; };
  double getFima(int tz);
  double *getFima() { return fima; };

private:
  int getNumIds();
  void setCappIndex();
  int indplane(int i, int j, int h);
  int indaxial(int k, int zz);
  void getXSbase(double Tm, double Tf);
  void getTcoefs(double Tm, double Tf, int kind);
  void getFluence();

  int prlev;
  GeomCAPP *geom;
  int ntri,nz,ndata;
  /*  for coord. translation */
  int nzinp,nrow;
  int *map;	// pointer to map
  int *zindbeg;
  sqlite3 *db;
  int efpd;
  int nids;
  int nwords;	// length of cross section
  int *tzs;
  double Tm_ref,Tf_ref,rt_ref;
  double *sigbase,*Tm_coef,*Tf_coef;
  double *flu,*fima;
};

