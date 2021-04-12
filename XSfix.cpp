/*
 * 	XSfix.cpp
 */
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "XSfix.h"

XSfix::XSfix()
{
  sigtr = NULL;
  sigab = NULL;
  signf = NULL;
  coefabM = NULL;
  coefnfM = NULL;
  coefabF = NULL;
  coefnfF = NULL;
};

XSfix::~XSfix()
{
  delete [] sigtr;
  delete [] sigab;
  delete [] signf;
  delete [] coefabM;
  delete [] coefnfM;
  delete [] coefabF;
  delete [] coefnfF;
};

void XSfix::set(int ndata_i)
{
  ndata = ndata_i;
  sigtr = new double[ndata];
  sigab = new double[ndata];
  signf = new double[ndata];
  coefabM = new double[2];
  coefnfM = new double[2];
  coefabF = new double[2];
  coefnfF = new double[2];

  for(int n=0; n < ndata; n++) {
    sigtr[n] = 0.3;
    sigab[n] = 4.0e-3;
    signf[n] = 4.5e-3;
  }
  coefabM[0] =  1.0e-6;  coefabM[1] =0.0;
  coefnfM[0] =  0.1e-6;  coefnfM[1] =0.0;
  coefabF[0] =  1.0e-4;  coefabF[1] =0.0;
  coefnfF[0] =  0.1e-4;  coefnfF[1] =0.0;
};

void XSfix::getSig(int m, double Tm, double Tf, double sig[])
{
  double x = Tm-800.0;
  double y = sqrt(Tf)-sqrt(800.0);
  sig[0] = sigtr[m];
  sig[1] = sigab[m] + x*(coefabM[0]+x*coefabM[1])
		    + y*(coefabF[0]+x*coefabF[1]);
  double snf = signf[m] + x*(coefnfM[0]+x*coefnfM[1])
		    + y*(coefnfF[0]+x*coefnfF[1]);
  assert(snf >= 0.0);
  sig[2] = snf;
};

double XSfix::getSigtr(int m, double Tm, double Tf)
{
  return sigtr[m];
};

double XSfix::getSigab(int m, double Tm, double Tf)
{
  double x = Tm-800.0;
  double y = sqrt(Tf)-sqrt(800.0);
  return sigab[m] + x*(coefabM[0]+x*coefabM[1])
		  + y*(coefabF[0]+x*coefabF[1]);
};

double XSfix::getSignf(int m, double Tm, double Tf)
{
  double x = Tm-800.0;
  double y = sqrt(Tf)-sqrt(800.0);
  double snf = signf[m] + x*(coefnfM[0]+x*coefnfM[1])
			+ y*(coefnfF[0]+x*coefnfF[1]);
  assert(snf > 0.0);
  return snf;
};

