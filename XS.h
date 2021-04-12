#ifndef _inc_XS
#define _inc_XS
/*
 * 	XS.h - generic cross section type
 */
class XS {
public:
  XS();
  virtual ~XS();
  virtual void getSig(int m, double Tm, double Tf, double sig[]);
  virtual double getSigtr(int m, double Tm, double Tf);
  virtual double getSigab(int m, double Tm, double Tf);
  virtual double getSignf(int m, double Tm, double Tf);
};
#endif
