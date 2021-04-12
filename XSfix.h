/*
 * 	XS.fix
 */
#include "XS.h"

class XSfix: public XS {
public:
  XSfix();
  ~XSfix();
  void set(int ndata);
  void getSig(int m, double Tm, double Tf, double sig[]);
  double getSigtr(int m, double Tm, double Tf);
  double getSigab(int m, double Tm, double Tf);
  double getSignf(int m, double Tm, double Tf);
private:
  int ndata;
  double *sigtr,*sigab,*signf;
  double *coefabM,*coefnfM;
  double *coefabF,*coefnfF;
};
