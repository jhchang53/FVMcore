/*
 * 	XScapp.h
 * 	derived class for CAPP cross section retrieve fro database
 *     Usage:
 * 	 void setGeom(GeomCAPP *geom); - set geometry for CAPP configuration
 *	   geom : CAPP geometry class
 *	 void open(const char *dbfn, int efpd, double Tm_ref, double Tf_ref);
 *	  - open cross section database
 *	   dbfn : database file name
 *	   efpd : efpd of the cross section sets
 *	   Tm_ref : reference mod T given to cross section
 *	   Tf_ref : reference fuel T given to cross section
 *	 void getSig(int tz, double Tm, double Tf, double sig[]);
 *	  - get cross section for element tz 
 *	   tz : element number (tz=t*nz+z)
 *	   Tm : moderator temp. (K)
 *	   Tf : fuel temp (K)
 *	   sig : returning cross section (sigtr/sigab/signf)
 *	 double getSigtr(int tz, double Tm, double Tf);
 *	   returns sigtr
 *	   tz : element number
 *	   Tm : moderator temp. (K)
 *	   Tf : fuel temp (K)
 *	 double getSigab(int tz, double Tm, double Tf);
 *	   returns sigab
 *	   tz : element number
 *         Tm : moderator temp. (K)
 *         Tf : fuel temp (K)
  double getSignf(int tz, double Tm, double Tf);
  double getFlu(int tz);
  double *getFlu() { return flu; };
  double getFima(int tz);
  double *getFima() { return fima; };

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

